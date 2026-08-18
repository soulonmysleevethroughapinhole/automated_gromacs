# this workflow will be used for acquiring the empirical structure, 
# their seq AlphaFold and ColabFold strs are generated 
# 
# the generated structures are screened 
# using RMSD analysis, to prevent MD simulations of structures that are too much alike
# 
# template similarity to empirical structure will also be collected, for further analysis
#TARGETS = ['7UU9']
#SOURCES = ['empirical', 'alphafold', 'colabfold']
TARGETS = config.get("pdbs")#, ["7UU9"])
SOURCES = config.get("prediction_sources")#, ["alphafold", "colabfold"])
MODELS = config.get("model_numbers", [1, 2, 3, 4, 5])

import json
import os
import shutil
import urllib.request
import zipfile 
import numpy as np
import re
from pathlib import Path
import subprocess


def aggregate_md_inputs(wildcards):
	resolved_paths = []
	
	# Loop over your targets inside the function
	for pdb in TARGETS:
		# 1. Force Snakemake to wait by explicitly passing the target name to the checkpoint
		checkpoint_output = checkpoints.screen_structure.get(pdb=pdb).output.report
		
		# 2. Read the evaluation data for this specific PDB
		with open(checkpoint_output) as f:
			data = json.load(f)
		
		# 3. Append the correct path to the master list based on structural compliance
		if data["rmsd"] < config["max_allowed_rmsd"]:
			# Target passed! Route it to the production pipeline
			target_length = config["targets"][pdb]["length_ns"]
			resolved_paths.append(f"results/2_simulations/{pdb}/{target_length}ns/md.gro")
		else:
			# Target failed! Route it to an alternative terminal endpoint
			resolved_paths.append(f"results/1_predicted_models/rejected_{pdb}.txt")
			
	return resolved_paths

#def aggregate_md_inputs(wildcards):
#    # 1. Force Snakemake to wait for the checkpoint to complete and evaluate
#    checkpoint_output = checkpoints.screen_structure.get(**wildcards).output.report
#    
#    # 2. Read the evaluation data
#    with open(checkpoint_output) as f:
#        data = json.load(f)
#    
#    # 3. Decide the DAG path based on structural compliance
#    if data["rmsd"] < config["max_allowed_rmsd"]:
#        # Target passed! Route it to the production pipeline
#        target_length = config["targets"][wildcards.pdb]["length_ns"]
#        return f"results/2_simulations/{wildcards.pdb}/{target_length}ns/md.gro"
#    else:
#        # Target failed! Route it to an alternative terminal endpoint
#        return f"results/1_predicted_models/rejected_{wildcards.pdb}.txt"


#rule all:
#    input:
#        # Dynamically resolved target paths
#        aggregate_md_inputs, pdb=TARGETS

# rule 1: GET the empirical pdb from rcsb
rule get_empirical_structure:
	output:
		cif_file = 'results/structures/{pdb}/empirical/native.cif'
	log:
		'logs/{pdb}/rcsb_download.log'
	run:
		url = f'https://files.rcsb.org/download/{wildcards.pdb}.cif'
		try:
			print(f'fetching {wildcards.pdb} from {url}')
			urllib.request.urlretrieve(url, output.cif_file)
		except Exception as e:
			with open(log[0], 'w') as log_file:
				log_file.write(f"Error fetching {wildcards.pdb}: {e}\n")
			raise e
# rule 1.25: fetch the sequence from rcsb
rule fetch_pdb_fasta:
	output:
		seq_file = "data/sequences/{pdb}_rcsb.fasta"
	log:
		"logs/{pdb}/fetch_fasta.log"
	run:
		url = f"https://www.rcsb.org/fasta/entry/{wildcards.pdb}"
		data = urllib.request.urlopen(url).read().decode("utf-8")
		with open(output.seq_file, "w") as f:
			f.write(data)
		with open(log[0], "w") as logf:
			logf.write(f"Fetched {wildcards.pdb} FASTA from RCSB\n")

# rule 1.5: extract data from cif file to get the sequence and other metadata
rule extract_cif_data:
	input:
		cif_file = 'results/structures/{pdb}/empirical/native.cif'
	output:
		#seq_file = 'data/sequences/{pdb}.fasta',
		metadata = 'results/structures/{pdb}/empirical/metadata.json'
	log:
		'logs/{pdb}/extract_cif_data.log'
	shell:
		#"python scripts/extract_cif_data.py -i {input.cif_file} -s {output.seq_file} -m {output.metadata} > {log} 2>&1"
		"python scripts/extract_cif_data.py -i {input.cif_file} -m {output.metadata} > {log} 2>&1"

#rule reconstruct_sequence:
#	input:
#		metadata = 'results/structures/{pdb}/empirical/metadata.json',
#		#seq_file = 'data/sequences/{pdb}.fasta'
#	output:
#		reconstructed_seq = 'data/sequences/{pdb}_reconstructed.fasta'
#	log:
#		'logs/{pdb}/reconstruct_sequence.log'
#	shell:
#		"python scripts/reconstruct_sequence.py -i {input.seq_file} -o {output.reconstructed_seq} > {log} 2>&1"

# rule 1.75: convert cif to pdb for downstream processing
rule convert_cif_to_pdb:
	input:
		cif_file = 'results/structures/{pdb}/empirical/native.cif'
	output:
		pdb_file = 'results/structures/{pdb}/empirical/native.pdb'
	log:
		'logs/{pdb}/convert_cif_to_pdb.log'
	shell:
	   "python scripts/convert_cif_to_pdb.py --input {input.cif_file} --output {output.pdb_file} > {log} 2>&1"

# rule 1.875: reconstruct sequence from cif and metadata
rule reconstruct_sequence:
	input:
		cif = 'results/structures/{pdb}/empirical/native.cif',
		metadata = 'results/structures/{pdb}/empirical/metadata.json'
	output:
		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
		residue_table = 'results/structures/{pdb}/empirical/{pdb}_canonical_residues.json'
	log:
		'logs/{pdb}/reconstruct_sequence.log'
	shell:
		"python scripts/reconstruct_sequence.py "
		"--input {input.cif} "
		"--metadata {input.metadata} "
		"--output {output.canonical_seq} "
		"--table-out {output.residue_table} > {log} 2>&1"
# rule 2: insert missing AAs 
rule clean_pdb_structure:
	input:
		cif_file = 'results/structures/{pdb}/empirical/native.cif',
		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
		residue_table = 'results/structures/{pdb}/empirical/{pdb}_canonical_residues.json'
	output:
		clean_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb'
	log:
		'logs/{pdb}/clean_structure.log'
	shell:
		"python scripts/clean_structure.py "
		"--input {input.cif_file} "
		"--sequence {input.canonical_seq} "
		"--residue-table {input.residue_table} "
		"--output {output.clean_pdb} > {log} 2>&1"

rule validate_str_seq_match:
	input:
		fixed_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
		fasta = "data/sequences/{pdb}_rcsb.fasta"
	output:
		flag = "results/structures/{pdb}/empirical/sequence_validation.ok"
	run:

		AA_3TO1 = {
			'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H',
			'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q',
			'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
		}

		# 1. Read Canonical FASTA
		fasta_seq = ""
		with open(input.fasta, "r") as f:
			for line in f:
				if not line.startswith(">"):
					fasta_seq += line.strip()

		# 2. Extract Sequence from Fixed PDB (using C-alpha backbone atoms)
		pdb_seq = []
		last_residue_id = None

		with open(input.fixed_pdb, "r") as f:
			for line in f:
				if line.startswith("ATOM") or line.startswith("HETATM"):
					atom_name = line[12:16].strip()
					if atom_name == "CA":
						res_name = line[17:20].strip()
						chain_id = line[21]
						res_seq = line[22:26].strip()
						
						# Unique residue identifier (chain + residue number)
						residue_id = f"{chain_id}_{res_seq}"
						
						if residue_id != last_residue_id:
							# Convert 3-letter code to 1-letter code (default to 'X' if non-standard)
							pdb_seq.append(AA_3TO1.get(res_name, "X"))
							last_residue_id = residue_id

		pdb_seq_str = "".join(pdb_seq)

		# 3. Validation Check
		if fasta_seq == pdb_seq_str:
			print(f"✅ [SUCCESS] {wildcards.pdb} structure sequence matches canonical FASTA exactly ({len(fasta_seq)} residues).")
			with open(output.flag, "w") as out:
				out.write(f"VALIDATED: Length={len(fasta_seq)}\n")
				out.write(f"FASTA: {fasta_seq}\n")
				out.write(f"PDB:   {pdb_seq_str}\n")
		else:
			# Report rich diagnostic info before failing
			print("\n❌ [ERROR] Sequence mismatch detected!")
			print(f"   Canonical FASTA length: {len(fasta_seq)}")
			print(f"   Fixed PDB seq length:   {len(pdb_seq_str)}")
			
			# Highlight mismatch position if lengths match but letters differ
			if len(fasta_seq) == len(pdb_seq_str):
				for i, (f_aa, p_aa) in enumerate(zip(fasta_seq, pdb_seq_str)):
					if f_aa != p_aa:
						print(f"   Mismatch at index {i+1} (Residue {i+1}): FASTA='{f_aa}' vs PDB='{p_aa}'")
			
			raise ValueError(
				f"Sequence of fixed structure {input.fixed_pdb} does not match {input.fasta}. "
				"Halting before querying AlphaFold."
			)

# rule 3: extract sequence of fixed empirical structure
#rule extract_sequence:
#    input:
#        clean_pdb = 'results/structures/{pdb}/empirical/clean_structure.pdb'
#    output:
#        seq_file = 'data/sequences/{pdb}.fasta'
#    log:
#        'logs/{pdb}/extract_sequence.log'
#    shell:
#        "python scripts/extract_sequence.py {input.clean_pdb} {output.seq_file} > {log} 2>&1"


# rule 4: generate protein predictions using alphafold or cf
#rule query_alphafold_prediction:
#	input:
#		clean_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
#		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
#	output:
#		zip_archive = 'results/predictions_raw/{pdb}/alphafold_raw.zip'
#	log:
#		'logs/{pdb}/api_query_alphafold.log'
#	run:
#		shell(
#			f"python scripts/query_alphafold_client.py "
#			f"--input {input.canonical_seq} --output {output.zip_archive}"
#		)

rule prepare_prediction_data:
	input:
		clean_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
		flag = "results/structures/{pdb}/empirical/sequence_validation.ok"
	output:
		metadata = 'results/predictions_raw/{pdb}/{source}_raw_settings.txt'
	log:
		'logs/{pdb}/api_query_{source}.log'
	run:
		manual_path = config.get("manual_prediction_dirs", {}).get(wildcards.source)

		with open(input.canonical_seq, "r") as handle:
			sequence = "".join(line.strip() for line in handle if not line.startswith(">"))

		if wildcards.source not in ["alphafold", "colabfold"]:
			raise ValueError(f"Unknown prediction source: {wildcards.source}")
		
		if wildcards.source == 'alphafold':
			settings = {
				"source": 'alphafold',
				"pdb": wildcards.pdb,
				"manual_prediction_dir": manual_path,
				"sequence_length": len(sequence),
				"template": 'No template',
			}
		elif wildcards.source == 'colabfold':
			settings = {
				"source": "colabfold",
				"pdb": wildcards.pdb,
				"manual_prediction_dir": manual_path,
				"sequence_length": len(sequence),
				'msa_mode': 'single_sequence',
				'max_msa': '16:32',
				'pair_mode': 'unpaired',
				'num_recycles': 1,
				'use_amber': False,
				'templates': 'None'
			}

		os.makedirs(os.path.dirname(output.metadata), exist_ok=True)
		with open(output.metadata, "w") as handle:
			handle.write(f"sequence:\n{sequence}\n\nsettings:\n")
			for key, value in settings.items():
				handle.write(f"{key}: {value}\n")

rule check_pred_archive:
	input:
		metadata = 'results/predictions_raw/{pdb}/{source}_raw_settings.txt'
	output:
		zip_archive = 'results/predictions_raw/{pdb}/{source}_raw.zip'
	log:
		'logs/{pdb}/check_{source}_archive.log'
	run:
		manual_path = config.get("manual_prediction_dirs", {}).get(wildcards.source) + wildcards.pdb + f'/{wildcards.source}_{wildcards.pdb}_raw.zip'
		if not manual_path:
			raise WorkflowError(
				f"No manual {wildcards.source} archive path configured. \n"
				"Run the notebook, export the zip, and place it in the expected location. \n"
			)
		if not os.path.exists(manual_path):
			raise WorkflowError(
				f"{wildcards.source} archive not found at {manual_path}. \n"
				"Please generate it and rerun the workflow. \n"
			)

		os.makedirs(os.path.dirname(output.zip_archive), exist_ok=True)
		shutil.copy(manual_path, output.zip_archive)

# NOTE: maybe considering getting alphafold database structures
#rule prepare_alphafold_prediction:
#	input:
#		clean_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
#		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
#	output:
#		metadata = 'results/predictions_raw/{pdb}/alphafold_raw_settings.txt'
#	log:
#		'logs/{pdb}/api_query_alphafold.log'
#	run:
#		manual_path = config.get("manual_prediction_dirs", {}).get("alphafold")
#
#		with open(input.canonical_seq, "r") as handle:
#			sequence = "".join(line.strip() for line in handle if not line.startswith(">"))
#
#
#		settings = {
#			"source": "alphafold",
#			"pdb": wildcards.pdb,
#			"manual_prediction_dir": manual_path,
#			"sequence_length": len(sequence),
#			"template": 'No template',
#		}
#
#
#
#		os.makedirs(os.path.dirname(output.metadata), exist_ok=True)
#		with open(output.metadata, "w") as handle:
#			handle.write(f"sequence:\n{sequence}\n\nsettings:\n")
#			for key, value in settings.items():
#				handle.write(f"{key}: {value}\n")
#
#rule prepare_colabfold_prediction:
#	input:
#		clean_pdb = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
#		canonical_seq = 'results/structures/{pdb}/empirical/{pdb}_canonical.fasta',
#	output:
#		metadata = 'results/predictions_raw/{pdb}/colabfold_raw_settings.txt'
#	log:
#		'logs/{pdb}/api_query_colabfold.log'
#	run:
#		manual_path = config.get("manual_prediction_dirs", {}).get("colabfold")
#
#		with open(input.canonical_seq, "r") as handle:
#			sequence = "".join(line.strip() for line in handle if not line.startswith(">"))
#
#		settings = {
#			"source": "colabfold",
#			"pdb": wildcards.pdb,
#			"manual_prediction_dir": manual_path,
#			"sequence_length": len(sequence),
#			'msa_mode': 'single_sequence',
#			'max_msa': '16:32',
#			'pair_mode': 'unpaired',
#			'num_recycles': 1,
#			'use_amber': False,
#			'templates': 'None'
#		}
#
#		os.makedirs(os.path.dirname(output.metadata), exist_ok=True)
#		with open(output.metadata, "w") as handle:
#			handle.write(f"sequence:\n{sequence}\n\nsettings:\n")
#			for key, value in settings.items():
#				handle.write(f"{key}: {value}\n")

#rule check_alphafold_archive:
#	input:
#		metadata = 'results/predictions_raw/{pdb}/alphafold_raw_settings.txt'
#	output:
#		zip_archive = 'results/predictions_raw/{pdb}/alphafold_raw.zip'
#	log:
#		'logs/{pdb}/check_alphafold_archive.log'
#	run:
#		manual_path = config.get("manual_prediction_dirs", {}).get("alphafold") + wildcards.pdb
#		if not manual_path:
#			raise WorkflowError(
#				"No manual alphafold archive path configured. "
#				"Run the notebook, export the zip, and place it in the expected location."
#			)
#		if not os.path.exists(manual_path):
#			raise WorkflowError(
#				f"alphafold archive not found at {manual_path}. "
#				"Please generate it and rerun the workflow."
#			)
#
#		os.makedirs(os.path.dirname(output.zip_archive), exist_ok=True)
#		shutil.copy(manual_path, output.zip_archive)


# rule 5: reorganize files after extraction 
rule extract_organize_models:
	input:
		#zip_archive = 'data/manual_predictions/{source}/{pdb}/{source}_{pdb}_raw.zip'
		zip_archive = 'results/predictions_raw/{pdb}/{source}_raw.zip',
	output:
		models = expand("results/structures/{{pdb}}/{{source}}/model_{model_num}.pdb", model_num=MODELS),
		metadata = 'results/structures/{pdb}/{source}/model_metadata.json'
	run:
		out_dir = f'results/structures/{wildcards.pdb}/{wildcards.source}'
		os.makedirs(out_dir, exist_ok=True)

		# temporary unpacking directory
		tmp_extract = Path(f'results/predictions_raw/{wildcards.pdb}/{wildcards.source}_tmp')
		tmp_extract.mkdir(parents=True, exist_ok=True)

		with zipfile.ZipFile(input.zip_archive, 'r') as zip_ref:
			zip_ref.extractall(tmp_extract)

		extracted_contents = list(tmp_extract.iterdir())
		if len(extracted_contents) == 1 and extracted_contents[0].is_dir():
			# ColabFold nested case: Point directly to the inner child folder
			data_dir = extracted_contents[0]
		else:
			# AlphaFold flat case: Use the base temp directory itself
			data_dir = tmp_extract

		# locate structures, jsons and map them to the correspponding model numbers
		# find model files and json files recursively (handle nested folders inside the zip)
		found_models = []
		json_files = []
		# NOTE: OLD
		#for root, dirs, files in os.walk(tmp_extract):
		#	for fname in files:
		#		lower = fname.lower()
		#		if lower.endswith(('.pdb', '.cif')):
		#			found_models.append(os.path.join(root, fname))
		#		elif lower.endswith('.json'):
		#			json_files.append(os.path.join(root, fname))
		
		# NOTE: NEW
		# in case of alphafold result, the extracted folder has .cif as children which correspond to the found_models
		# in case of colabfold result, the extracted tmp fld will have a .result folder as child, which contains the .pdb 
		# -> the numbering of the models will be determined by the .pdb / .cif files,
		# --> by regexing the string of the model, the full_data and scores can be extracted
		# colabfold models are 1 indexed, while alphafold models are 0 indexed  
		# from what I've seen, but ensure num of models checks out, and that if the numbering is shifted (in case of af, the model numbers should just be shifted by 1 to match the config numbering)

		
		for file in os.listdir(data_dir):
			if wildcards.source == "alphafold":
				str_ext = '.cif'
			elif wildcards.source == "colabfold":
				str_ext = '.pdb'
			else:
				raise WorkflowError(f"Unknown source {wildcards.source} for model extraction")

			if file.lower().endswith(str_ext):
				found_models.append(os.path.join(data_dir, file))
			if file.lower().endswith('.json'):
				json_files.append(os.path.join(data_dir, file))
		found_models = sorted(found_models, key = lambda x: x.split('_model_')[1])
		json_files = sorted(json_files)
		

		# map found model files to requested MODELS order; convert CIF -> PDB when needed
		metrics_summary = {}
		# placeholder for any top-level job/request info found in the archive
		
		job_request = None
		for json_pth in json_files:
			name = json_pth.lower()
			if name.endswith('config.json') or name.endswith('job_request.json'):
				job_request = json_pth
				json_files.remove(json_pth)
				break
		#match wildcards.source:
		#case 'alphafold':
		#	job_request = '<>_job_request.json'

		if job_request != None:
			with open(job_request, 'r') as f:
				metrics_summary['job_request'] = json.load(f)
		else:
			raise WorkflowError(f"No job_request.json or config.json found in archive {input.zip_archive}")	

		if len(found_models) == 0:
			raise WorkflowError(f"No model files (.pdb or .cif) found in archive {input.zip_archive}")
		if len(found_models) != len(MODELS):
			raise WorkflowError(f"Warning: Number of model files found ({len(found_models)}) does not match expected number ({len(MODELS)}) for {wildcards.source}/{wildcards.pdb}. Not proceeding with available models.")

		for idx, model_num in enumerate(MODELS):
			src_path = found_models[idx]

			orig_name = os.path.basename(src_path)
			target_path = os.path.join(out_dir, f"model_{model_num}.pdb")
			# ensure output directory exists
			os.makedirs(os.path.dirname(target_path), exist_ok=True)

			if orig_name.lower().endswith('.pdb'):
				shutil.move(src_path, target_path)
			elif orig_name.lower().endswith('.cif'):
				# convert CIF -> PDB using existing conversion script
				subprocess.check_call([
					"python", "scripts/convert_cif_to_pdb.py",
					"--input", src_path,
					"--output", target_path
				])

			if wildcards.source == 'alphafold':
				# Alphafold models are 0-indexed, so they should match the index directly
				# Extract the model number from the filename to verify
				#extracted_model_num = None

				match = re.search(r"^(.*?)(?:_model_(\d+))$", Path(orig_name).stem)
				#base_name = None
				if match:
					extracted_model_num = int(match.group(2))
					base_name = match.group(1)
					if extracted_model_num != idx:
						raise WorkflowError(
							f"Alphafold model numbering mismatch: expected model_{model_num} at index {idx}, "
							f"but found {os.path.basename(src_path)}. Please check the archive."
						)
				else:
					raise WorkflowError(f"Could not extract model number from filename {os.path.basename(src_path)} for Alphafold source.")

				model_key = f"model_{model_num}"
				metrics_summary[model_key] = {"original_file": orig_name}


				# select the metadata file corresponding to this model
				fd_target_stem = f"{base_name}_full_data_{extracted_model_num}"
				matching_json = next(
					(p for p in json_files if Path(p).stem == fd_target_stem), 
					None
				)
				
				if matching_json:
					print(f"Found match: {matching_json}")
					with open(matching_json, 'r') as f:
						metrics_summary[model_key]['metadata_full_data'] = json.load(f)
				else:
					raise WorkflowError(f"No JSON found matching stem: {fd_target_stem} for Alphafold model {extracted_model_num}. Please check the archive.")
					#print(f"No JSON found matching stem: {target_stem}")

				sc_target_stem = f'{base_name}_summary_confidences_{extracted_model_num}' 
				matching_json = next(
					(p for p in json_files if Path(p).stem == sc_target_stem), 
					None
				)

				if matching_json:
					print(f"Found match: {matching_json}")
					with open(matching_json, 'r') as f:
						metrics_summary[model_key]['metadata_summary_confidences'] = json.load(f)
				else:
					raise WorkflowError(f"No JSON found matching stem: {sc_target_stem} for Alphafold model {extracted_model_num}. Please check the archive.")
			elif wildcards.source == 'colabfold':
				model_match = re.search(r"_model_(\d+)", Path(orig_name).stem)
				extracted_model_num = int(model_match.group(1)) if model_match else None

				match = re.search(r"^(.*?)_unrelaxed_(rank_.*)$", Path(orig_name).stem)

				if match:
					#extracted_model_num = int(match.group(2))
					#base_name = match.group(1)
					if extracted_model_num != model_num:
						raise WorkflowError(
							f"Colabfold model numbering mismatch: expected model_{model_num} at index {idx}, "
							f"but found {os.path.basename(src_path)}. Please check the archive."
						)

					base_prefix = match.group(1)
					shared_suffix = match.group(2)

					expected_json_stem = f"{base_prefix}_scores_{shared_suffix}" 
				else:
					raise WorkflowError(f"Could not extract model number from filename {os.path.basename(src_path)} for Colabfold source.")

				model_key = f"model_{model_num}"
				metrics_summary[model_key] = {"original_file": orig_name}

				# select the metadata file corresponding to this model
				matching_json = next(
					(p for p in json_files if Path(p).stem == expected_json_stem),
					None
				)

				if matching_json:
					print(f"Found match: {matching_json}")
					with open(matching_json, 'r') as f:
						metrics_summary[model_key]['metadata_scores'] = json.load(f)
				else:
					raise WorkflowError(f"No JSON found matching stem: {expected_json_stem} for Colabfold model {src_path}. Please check the archive.")
			else:
				raise WorkflowError(f"Unknown source {wildcards.source} for model extraction")

		shutil.rmtree(tmp_extract)

		with open(output.metadata, 'w') as f:
			json.dump(metrics_summary, 
				f, 
				indent=4,
				sort_keys=True
			)



# rule 6: rmsd screening
# Register for your get_final_targets function
checkpoint evaluate_predictions:
	input:
		clean_empirical = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
		models = expand("results/structures/{{pdb}}/{{source}}/model_{model_num}.pdb", model_num=MODELS)
	output:
		eval_dir = directory("results/structures/{pdb}/{source}/evaluations")
	run:
		# delegate batch RMSD evaluation to scripts/screen_rmsd.py
		os.makedirs(output.eval_dir, exist_ok=True)

		cmd = [
			"python", "scripts/screen_rmsd.py",
			"--input-emp", input.clean_empirical,
			"--output-dir", output.eval_dir,
			"--threshold", str(config.get("rmsd_min_threshold", 1.0))
		]

		# append model paths
		cmd.extend(list(input.models))

		subprocess.check_call(cmd)

''' TO EVALUATE PREDICTIONS, this gets the B factor col of af models to get the mean pLDDT score for each model
plddts = []
with open(input.pdb_file, "r") as pdb:
	for line in pdb:
		if line.startswith("ATOM") or line.startswith("HETATM"):
			# B-factor is located from index position 60 to 66 in standard PDB format
			b_factor = float(line[60:66].strip())
			plddts.append(b_factor)

mean_plddt = sum(plddts) / len(plddts) if plddts else 0.0
'''