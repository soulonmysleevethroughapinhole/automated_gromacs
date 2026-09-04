
import os
import shutil
import tempfile
import subprocess
import logging

rule create_gmx_em_job:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	output:
		job_description = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_em_job.job"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/energy_minimization/create_job.log"
	shell:
		"python scripts/smk/gmx_em/create_gmx_em_job.py {input.tar} {output.job_description} > {log} 2>&1"


rule run_HPC_gmx_em:
	input:
		job_description = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_em_job.job"
	output:
		tar = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz",
		em_completed = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/em_completed.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/energy_minimization/run_HPC_gmx_em.log"
	shell:
		"python scripts/smk/gmx_em/run_HPC_gmx_em.py {input.job_description} {output.em_completed} > {log} 2>&1"

rule energy_minimize_structure_single:
	input:
		structure = "results/structures/{pdb}/{source}/{model_id}.pdb",
		mdp_st = "config/gmx_em/emw_steep_hydr.mdp",
		mdp_cg = "config/gmx_em/emw_cg_hydr.mdp"
	output:
		em_structure = "results/gmx_em/{pdb}/{source}/em_{model_id}.pdb"
	log:
		"logs/{pdb}/{source}/{model_id}/energy_minimize_structure.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/energy_minimize_structure.log"
		)
	threads: 12
	run:
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		# Setup Logger
		logger = logging.getLogger(f"em_single_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}")
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
		fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
		ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}"
		logger.info("Starting local energy minimization for structure: %s", target_id)

		input_pdb_abs = os.path.abspath(input.structure)
		mdp_st_abs = os.path.abspath(input.mdp_st)
		mdp_cg_abs = os.path.abspath(input.mdp_cg)
		output_pdb_abs = os.path.abspath(output.em_structure)

		os.makedirs(os.path.dirname(output_pdb_abs), exist_ok=True)

		# Execute inside a isolated temporary folder
		with tempfile.TemporaryDirectory() as work_dir:
			logger.info("Created temporary execution directory: %s", work_dir)

			def run_cmd(cmd_str, input_text=None):
				logger.info("Running: %s", cmd_str)
				res = subprocess.run(
					cmd_str,
					shell=True,
					cwd=work_dir,
					input=input_text,
					text=True,
					capture_output=True
				)
				with open(log_path, "a") as lf:
					lf.write(f"=== Command: {cmd_str} ===\n")
					lf.write(f"STDOUT:\n{res.stdout}\n")
					lf.write(f"STDERR:\n{res.stderr}\n\n")

				if res.returncode != 0:
					err_msg = f"Command failed with exit code {res.returncode}: {cmd_str}\n{res.stderr}"
					logger.error(err_msg)
					raise RuntimeError(err_msg)

			# 1. Topology & Hydrogen Generation (pdb2gmx)
			run_cmd(f"gmx pdb2gmx -f '{input_pdb_abs}' -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh")

			# 2. Box Definition & Solvation
			run_cmd("gmx editconf -f processed.gro -o newbox.gro -d 1.0 -bt cubic")
			run_cmd("gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top")

			# 3. Add Ions (Neutralization)
			run_cmd(f"gmx grompp -f '{mdp_st_abs}' -c solv.gro -r solv.gro -o ion_prep.tpr -p topol.top -maxwarn 2")
			run_cmd("gmx genion -s ion_prep.tpr -o ion_b4em.gro -p topol.top -pname NA -nname CL -neutral", input_text="SOL\n")

			# 4. Stage 1: Steepest Descent Minimization
			#run_cmd(f"gmx grompp -f '{mdp_st_abs}' -c ion_b4em.gro -r ion_b4em.gro -o st.tpr -p topol.top -maxwarn 2")
			#run_cmd(f"gmx mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log -ntomp {threads}")
			#run_cmd(f"gmx mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log -ntmpi 1 -ntomp {threads} -nb gpu -pme gpu -bonded gpu")
			run_cmd(f"gmx grompp -f '{mdp_st_abs}' -c ion_b4em.gro -r ion_b4em.gro -o st.tpr -p topol.top -maxwarn 2")
			run_cmd(f"gmx mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log -ntmpi 1 -ntomp {threads} -nb gpu -pme cpu -bonded cpu")

			# 5. Stage 2: Conjugate Gradient Minimization
			#run_cmd(f"gmx grompp -f '{mdp_cg_abs}' -c after_st.gro -r after_st.gro -o cg.tpr -p topol.top -maxwarn 2")
			#run_cmd(f"gmx mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log -ntomp {threads}")
			#run_cmd(f"gmx mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log -ntmpi 1 -ntomp {threads} -nb gpu -pme gpu -bonded gpu")
			run_cmd(f"gmx grompp -f '{mdp_cg_abs}' -c after_st.gro -r after_st.gro -o cg.tpr -p topol.top -maxwarn 2")
			run_cmd(f"gmx mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log -ntmpi 1 -ntomp {threads} -nb gpu -pme cpu -bonded cpu")
			
			# 6. Extract Final Cleaned Protein Structure
			run_cmd("gmx trjconv -s cg.tpr -f after_cg.gro -o output_ec.pdb -pbc mol -ur compact", input_text="Protein\n")

			produced_pdb = os.path.join(work_dir, "output_ec.pdb")
			if not os.path.exists(produced_pdb):
				raise FileNotFoundError(f"Expected output structure not generated: {produced_pdb}")

			shutil.copy2(produced_pdb, output_pdb_abs)
			logger.info("✅ Successfully generated energy-minimized structure: %s", output_pdb_abs)