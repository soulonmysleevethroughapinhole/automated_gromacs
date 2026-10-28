# this file will be used for regenerating .mdp with energygrps
# set to Protein SOL, then mdrun -rerun ing the trajectories
# so that intramoleculer energy can be calculated from the 
# extracted energy components
# This should allow the identification native-like states

import os
import shutil
import pandas as pd

ruleorder: schedule_std100_rerun > schedule_custom_rerun

def get_rerun_energy_targets(config):
	targets = []

	for pdb in config.get("pdbs", []):
		# handle empiricals
		empirical_base = f'results/gromacs/{pdb}/empirical/canonical_structure'
		# Gather custom protocols and standard_100ns
		custom_entries = config.get("custom_simulations", {}).get(pdb, {}).get("empirical_canonical_structure", [])
		protocols = [e["protocol"] for e in custom_entries if "protocol" in e]
		if not protocols:
			protocols = ["standard_100ns"]
			
		for proto in protocols:
			done_file = f"{empirical_base}/{proto}/md_results/md_completed.txt"
			if os.path.exists(done_file):
				targets.append(f"data/gmx_rerun/{pdb}/empirical/canonical_structure/{proto}/gmx_rerun_energies_parsed.csv")
		#parsed_energies = "data/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/"

		# 2. Prediction sources (alphafold, colabfold, etc.)
		for source in config.get("prediction_sources", []):
			# Scan filesystem directly for completed prediction runs
			pred_dir = f"results/gromacs/{pdb}/{source}"
			if os.path.exists(pred_dir):
				for model_id in os.listdir(pred_dir):
					model_path = os.path.join(pred_dir, model_id)
					if os.path.isdir(model_path):
						for proto in os.listdir(model_path):
							done_file = os.path.join(model_path, proto, "md_results", "md_completed.txt")
							if os.path.exists(done_file):
								targets.append(f"data/gmx_rerun/{pdb}/{source}/{model_id}/{proto}/gmx_rerun_energies_parsed.csv")

	return targets

def parse_mm_energy_xvg(xvg_file):
	headers = []
	data = []
	
	with open(xvg_file, 'r') as f:
		for line in f:
			line = line.strip()
			# Capture column legend names dynamically from @ s0 legend "..."
			if line.startswith('@ s') and 'legend' in line:
				label = line.split('legend')[-1].replace('"', '').strip()
				headers.append(label)
			elif line.startswith(('@', '#')):
				continue
			elif line:
				parts = [float(x) for x in line.split()]
				if parts:
					data.append(parts)
					
	if not data:
		raise ValueError(f"No numeric data found in {xvg_file}")
		
	num_cols = len(data[0])
	# Assign header names or default if legends were missing
	if headers and len(headers) == num_cols - 1:
		columns = ['Time_ps'] + headers
	else:
		columns = ['Time_ps'] + [f"col_{i}" for i in range(1, num_cols)]
		
	return pd.DataFrame(data, columns=columns)


rule schedule_custom_rerun:
	input:
		custom_mdp = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/generated.mdp",
		gro = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_{protocol}_md.gro",
		#top = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/topol.top",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		xtc = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_{protocol}_md.xtc",
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/md_completed.txt",
	output:
		modified_mdp = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/protsolenergygrps.mdp",
		gro = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.gro",
		top = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.top",
		xtc = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.xtc",
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/logs/schedule_custom_rerun.log"
	run:
		os.makedirs(os.path.dirname(output.modified_mdp), exist_ok=True)
		os.makedirs(os.path.dirname(log[0]), exist_ok=True)

		with open(log[0], 'w') as log_file:
			log_file.write(f"Pre-processing MDP file for rerun: {input.custom_mdp}\n")
			
			with open(input.custom_mdp, 'r') as f:
				mdp_lines = f.readlines()

			cleaned_lines = []
			annealing_keys = ("annealing", "annealing_npoints", "annealing_time", "annealing_temp")
			
			for line in mdp_lines:
				stripped_line = line.strip()
				if not stripped_line or stripped_line.startswith(";"):
					cleaned_lines.append(line)
					continue

				key = stripped_line.split("=")[0].strip().replace("-", "_")
				
				if key in annealing_keys:
					log_file.write(f"Stripped simulated annealing directive: {line}")
					continue
					
				cleaned_lines.append(line)

			with open(output.modified_mdp, 'w') as f:
				f.write("; --- RERUN OVERRIDES ---\n")
				f.write("energygrps = Protein SOL\n")
				f.write("annealing  = no\n")
				f.write("; -----------------------\n\n")
				f.writelines(cleaned_lines)
				
			# KEEP INSIDE THE WITH BLOCK:
			log_file.write(f"Successfully generated {output.modified_mdp}\n")

		# Operations outside log_file scope:
		shutil.copy(input.gro, output.gro)
		shutil.copy(input.top, output.top)
		shutil.copy(input.xtc, output.xtc)

rule schedule_std100_rerun:
	input:
		mdp = 'config/gromacs_settings/protSOLenergygrp/md.mdp',
		gro = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.gro",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		xtc = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.xtc",
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/md_completed.txt",
	output:
		mdp = "results/gmx_rerun/{pdb}/{source}/{model_id}/standard_100ns/protsolenergygrps.mdp",
		gro = "results/gmx_rerun/{pdb}/{source}/{model_id}/standard_100ns/original.gro",
		top = "results/gmx_rerun/{pdb}/{source}/{model_id}/standard_100ns/original.top",
		xtc = "results/gmx_rerun/{pdb}/{source}/{model_id}/standard_100ns/original.xtc",
	run:
		os.makedirs(os.path.dirname(output.mdp), exist_ok=True)

		shutil.copy(input.mdp, output.mdp)
		shutil.copy(input.gro, output.gro)
		shutil.copy(input.top, output.top)
		shutil.copy(input.xtc, output.xtc)

rule perform_rerun:
	input:
		mdp = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/protsolenergygrps.mdp",
		gro = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.gro",
		top = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.top",
		xtc = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/original.xtc",
	output:
		tpr = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy.tpr",
		edr = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy.edr",
		log = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy.log",
	log:
		"results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/logs/perform_rerun.log"
	params:
		out_prefix = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy"
	threads: 4
	shell:
		"""
			mkdir -p $(dirname {log})
	
			gmx grompp -f {input.mdp} -c {input.gro}  -p {input.top} -o {output.tpr} -maxwarn 1 > {log} 2>&1
			gmx mdrun -s {output.tpr} -rerun {input.xtc} -deffnm {params.out_prefix} -nt {threads} > {log} 2>&1
		"""

rule extract_rerun_energies:
	input:
		tpr = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy.tpr",
		edr = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy.edr",
	output:
		energies = "data/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/gmx_rerun_energies.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/extract_gmx_rerun_energies.log"
	shell:
		"""
		mkdir -p $(dirname {output.energies})
		echo "1 2 3 4 5 9 15 16 17 18 19 20 0" | gmx energy \
			-s {input.tpr} \
			-f {input.edr} \
			-o {output.energies} > {log} 2>&1
		"""


# parse rerun xvg
rule parse_rerun_MM_energies:
	input:
		energies = "data/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/gmx_rerun_energies.xvg"
	output:
		parsed_energies = "data/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/gmx_rerun_energies_parsed.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/parse_gmx_rerun_energies.log"
	run:
		try:
			df_energies = parse_mm_energy_xvg(input.energies)
			df_energies.to_csv(output.parsed_energies, index=False)

			with open(log[0], "w") as l:
				l.write(f"Successfully parsed energies into {output.parsed_energies}\n")

		except Exception as e:
			with open(log[0], "w") as l:
				l.write(f"Error parsing energies: {str(e)}\n")
			raise e
