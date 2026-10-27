# this file will be used for regenerating .mdp with energygrps
# set to Protein SOL, then mdrun -rerun ing the trajectories
# so that intramoleculer energy can be calculated from the 
# extracted energy components
# This should allow the identification native-like states

import os
import shutil

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
				targets.append(f"results/gmx_rerun/{pdb}/empirical/canonical_structure/{proto}/rerun_energy.edr")

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
								targets.append(f"results/gmx_rerun/{pdb}/{source}/{model_id}/{proto}/rerun_energy.edr")

	return targets

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
	run:
		os.makedirs(os.path.dirname(output.modified_mdp), exist_ok=True)

		newline = "energygrps = Protein SOL\n"
		with open(input.custom_mdp, 'r') as f:
			mdp_lines = f.read()

		with open(output.modified_mdp, 'w') as f:
			f.write(newline)
			f.write(mdp_lines)

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
	params:
		out_prefix = "results/gmx_rerun/{pdb}/{source}/{model_id}/{protocol}/rerun_energy"
	threads: 4
	shell:
		"""
			gmx grompp -f {input.mdp} -c {input.gro}  -p {input.top} -o {output.tpr} -maxwarn 1
			gmx mdrun -s {output.tpr} -rerun {input.xtc} -deffnm {params.out_prefix} -nt {threads}
		"""