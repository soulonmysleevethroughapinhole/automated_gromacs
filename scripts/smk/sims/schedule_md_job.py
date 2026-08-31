#!/usr/bin/env python3
import subprocess
import argparse
import os
import shutil
import yaml

with open("config/config.yml") as f:
	try:
		config = yaml.safe_load(f)
	except:
		raise NameError('Config not found for schedule_md_job.py')	

def parse_args():
	parser = argparse.ArgumentParser(
		description="Decide the scheduling for MD job."
	)
	parser.add_argument(
		"--in_top",
		"-t",
		required=True,
		help="Path to the input top file.",
	)
	parser.add_argument(
		'--in_gro',
		'-g',
		required=True,
		help='Path to input gro'
	)
	parser.add_argument(
		'--in_mdp',
		'-m',
		required=True,
		help='Path to input mdp'
	)
	parser.add_argument(
		"--out_scheduling",
		"-sch",
		required=True,
		help="Path to the output scheduling file.",
	)
	parser.add_argument(
		'--out_tpr',
		'-o',
		required=True,
		help='Path to the output tpr'
	)
	parser.add_argument(
		'--pdb',
		'-p',
		required=True,
		help='PDB of model'
	)
	parser.add_argument(
		'--source',
		'-src',
		required=True,
		help='Source of model'
	)
	parser.add_argument(
		'--model_id',
		'-mid',
		required=True,
		help='Model ID'
	)
	parser.add_argument(
		'--abs_log',
		'-l',
		required=True,
		help='Absolute path of log path'
	)
	parser.add_argument(
		'--protocol',
		'-prot',
		required=True,
		help='Which protocol for MD to use'
	)
	return parser.parse_args()

def main():
	args = parse_args()

	pdb = args.pdb
	source = args.source
	model_id = args.model_id
	#current_protocol = "standard"
	current_protocol = args.protocol

	in_top = args.in_top
	in_gro = args.in_gro
	in_mdp = args.in_mdp

	out_tpr = args.out_tpr
	out_scheduling = args.out_scheduling
	abs_log = args.abs_log

	mdp_abs = os.path.abspath(in_mdp)

	def_compute_target = config.get("default_compute_target", "local")
	# get target from config for this PDB, or use default
	pdb_config_entries = config.get("custom_simulations", {})\
					.get(pdb, {})\
					.get(f'{source}_{model_id}', [])
	pdb_config = next((entry for entry in pdb_config_entries if entry.get("protocol") == current_protocol), {})
	compute_target = pdb_config.get("compute_target", def_compute_target)
	#compute_target = pdb_config.get("compute_target", None)
	work_dir = os.path.dirname(out_scheduling)
	os.makedirs(work_dir, exist_ok=True)
	# Copy the necessary files into the work directory
	shutil.copy(in_gro, os.path.join(work_dir, os.path.basename(in_gro)))
	shutil.copy(in_top, os.path.join(work_dir, os.path.basename(in_top)))
	shutil.copy(in_mdp, os.path.join(work_dir, os.path.basename(in_mdp)))
	exec_dir = os.path.dirname(out_scheduling)
	# TPR path will be created inside exec_dir with the basename of output.tpr
	tpr_path = os.path.join(exec_dir, os.path.basename(out_tpr))
	# absolute path for grompp logs
	#log_abs = os.path.join(exec_dir, "grompp.log")
	# Compile the .tpr ONLY if it doesn't exist yet
	if not os.path.exists(tpr_path):
		subprocess.run("""
			cd {exec_dir}
			gmx grompp -f {mdp_abs} \
				-o $(basename {out_tpr}) \
				-c $(basename {in_gro}) \
				-r $(basename {in_gro}) \
				-p $(basename {in_top}) -maxwarn 1 > {abs_log} 2>&1
		""")
	# Clean up the temporary copies we made in exec_dir
	os.remove(os.path.join(exec_dir, os.path.basename(in_gro)))
	os.remove(os.path.join(exec_dir, os.path.basename(in_top)))
	os.remove(os.path.join(exec_dir, os.path.basename(in_mdp)))
	
	scheduling_info = {}
	if compute_target in ["local", 'HPC']: # write local scheduling into job description
		scheduling_info['COMPUTE'] = compute_target
		with open (out_scheduling, 'w') as f:
			yaml.safe_dump(scheduling_info, f, sort_keys=False)
	else:
		raise ValueError(f"{compute_target} Compute type not permitted ")

if __name__ == '__main__':
	main()