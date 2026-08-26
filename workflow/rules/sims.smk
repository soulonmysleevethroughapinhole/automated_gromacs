import os
import shutil
import tarfile
import yaml
from dotenv import load_dotenv
load_dotenv() # Load environment variables from .env file

import logging
import sys
import traceback
import re
import subprocess
import time

HPC_HOST = os.getenv("HPC_HOST")
HPC_USER = os.getenv("HPC_USER")
HPC_ACCOUNT = os.getenv("HPC_ACCOUNT")
HPC_REMOTE_BASE = os.getenv("HPC_REMOTE_BASE")
HPC_PARTITION = os.getenv("HPC_PARTITION")
HPC_GMX_HOME = os.getenv("HPC_GMX_HOME")
HPC_TIME = os.getenv("HPC_TIME")
SUBMIT_HPC = os.getenv("SUBMIT_HPC", "0")
SSH_DAEMON_PORT = os.getenv("SSH_DAEMON_PORT", "22")  # Default to port 22 if not set


# Rule 1.75: Prepare protein and  generate topology
rule prepare_system:
	input:
		#pdb_clean = "results/proteins/{pdb}/cleaned_protein.pdb",
		pdb_file = "results/structures/{pdb}/{source}/{model_id}.pdb",
	output:
		gro_processed = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/{pdb}_{source}_{model_id}_processed.gro",
		gro_box       = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/{pdb}_{source}_{model_id}_newbox.gro",
		gro_solv      = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/{pdb}_{source}_{model_id}_solv.gro",
		top           = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",

	params:
		water = "tip3p",
		ff = "amber99sb-ildn",
		pdb_abs = lambda wildcards: os.path.abspath(
			f"results/structures/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}.pdb"
		),
		# NOTE: maybe lambda basenames for output files
		#log_abs = lambda log: os.path.abspath(str(log))
		#log_abs = lambda wildcards, log: os.path.abspath(str(log[0]))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/pdb2gmx.log"
		)
	log:
		"logs/{pdb}/{source}/{model_id}/pdb2gmx.log"
	shell:
		"""
		# absolute path of PDB
		# PDB_ABS_path = $(readlink -f {input.pdb_file})
		# absolute path of log file
		#LOG_ABS=$(readlink -f {log})

		mkdir -p results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns
		cd results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns

		gmx pdb2gmx -water {params.water} -ff {params.ff} -ignh \
			-f {params.pdb_abs} \
			-o $(basename {output.gro_processed}) \
			-p $(basename {output.top}) >> {params.log_abs} 2>&1

		gmx editconf -f $(basename {output.gro_processed}) -o $(basename {output.gro_box}) -d 1.0 -bt cubic >> {params.log_abs} 2>&1
		gmx solvate -cp $(basename {output.gro_box}) -cs spc216.gro -o $(basename {output.gro_solv}) -p $(basename {output.top}) >> {params.log_abs} 2>&1


		"""
# Rule 2.1 add ions to neutralize 
rule add_ions:
	input:
		gro_solv = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/{pdb}_{source}_{model_id}_solv.gro",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		mdp = 'config/gromacs_settings/interruptable_config_ultimate/emw_steep.mdp'
	output:
		tpr_em = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/em_setup.tpr",
		gro_ions = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/ion_b4em.gro",
	log:
		"logs/{pdb}/{source}/{model_id}/genion.log"
	params:
		mdp_abs = lambda wildcards, input: os.path.abspath(input.mdp),
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/genion.log"
		)

	#LOG_ABS=$(readlink -f {log})
	shell:
		"""
		#cd results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns
		cd $(dirname {output.gro_ions})

		gmx grompp -v -f {params.mdp_abs} \
			-c $(basename {input.gro_solv}) \
			-o $(basename {output.tpr_em}) \
			-p $(basename {input.top}) -maxwarn 1 >> {params.log_abs} 2>&1
			
		echo "SOL" | gmx genion -s $(basename {output.tpr_em}) \
			-o $(basename {output.gro_ions}) \
			-p $(basename {input.top}) -pname NA -nname CL -neutral >> {params.log_abs} 2>&1
		"""

# 3 st minimization
rule minimize_steepest:
	input:
		gro_ions = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/ion_b4em.gro",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		mdp = 'config/gromacs_settings/interruptable_config_ultimate/emw_steep.mdp',

	output:
		gro_st = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_st.gro",
		tpr_st = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/st.tpr",
		trr_st = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/st.trr",
		log_st = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/st.log",
	log:
		"logs/{pdb}/{source}/{model_id}/gmx_mdrun_steep.log"
	params:
		mdp_abs = lambda wildcards, input: os.path.abspath(input.mdp),
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/gmx_mdrun_steep.log"
		)
	shell:
		"""
		#LOG_ABS=$(readlink -f {log})

		cd results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns

		gmx grompp -v -f {params.mdp_abs} \
			-c $(basename {input.gro_ions}) \
			-o $(basename {output.tpr_st}) \
			-p $(basename {input.top}) -maxwarn 1 >> {params.log_abs} 2>&1
		
		gmx mdrun -v -ntmpi 1 -s $(basename {output.tpr_st}) \
			-o $(basename {output.trr_st}) \
			-c $(basename {output.gro_st}) \
			-g $(basename {output.log_st}) >> {params.log_abs} 2>&1
		"""
#rule 4: cg minimization
rule minimize_conjugate:
	input:
		gro_st = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_st.gro",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		mdp = 'config/gromacs_settings/interruptable_config_ultimate/emw_cg.mdp',
	output:
		tpr_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/cg.tpr",
		gro_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		trr_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/cg.trr",
		log_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/cg.log",
	log:
		"logs/{pdb}/{source}/{model_id}/gmx_mdrun_cg.log"
	params:
		mdp_abs = lambda wildcards, input: os.path.abspath(input.mdp),
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/gmx_mdrun_cg.log"
		)

	shell:
		"""
		#LOG_ABS=$(readlink -f {log})

		# cd results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns
		cd $(dirname {output.tpr_cg})

		gmx grompp -v -f {params.mdp_abs} \
			-c $(basename {input.gro_st}) \
			-o $(basename {output.tpr_cg}) \
			-p $(basename {input.top}) -maxwarn 1 >> {params.log_abs} 2>&1
		
		gmx mdrun -v -ntmpi 1 -s $(basename {output.tpr_cg}) \
			-o $(basename {output.trr_cg}) \
			-c $(basename {output.gro_cg}) \
			-g $(basename {output.log_cg}) >> {params.log_abs} 2>&1
		"""

# NOTE:
# NOTE: Maybe scheduling is not necessary!!!
# NOTE: it's possible to just ssh into the server when md is ran and thats it
# NOTE: there'd just need to be a step to check if it's on the server first, if not 
# NOTE: then package the JOB folder, and tf it to server and extract to run the mdrun .jobSS
# NOTE:
# decide whether to submit job to HPC or run locally, based on environment variable
rule schedule_md_job:
	input:
		gro_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		mdp = "config/gromacs_settings/interruptable_config_ultimate/md.mdp",
	output:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/scheduling.yml",
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
	log:
		"logs/{pdb}/{source}/{model_id}/schedule_md_job.log"
	params:
		mdp_abs = lambda wildcards, input: os.path.abspath(input.mdp),
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/schedule_md_job.log"
		)

	run:
		def_compute_target = config.get("default_compute_target", "local")
		# get target from config for this PDB, or use default
		pdb_config_entries = config.get("custom_simulations", {})\
						.get(wildcards.pdb, {})\
						.get(f'{wildcards.source}_{wildcards.model_id}', [])

		current_protocol = "standard"
		pdb_config = next((entry for entry in pdb_config_entries if entry.get("protocol") == current_protocol), {})

		compute_target = pdb_config.get("compute_target", def_compute_target)
		#compute_target = pdb_config.get("compute_target", None)
		work_dir = os.path.dirname(output.scheduling)
		os.makedirs(work_dir, exist_ok=True)

		# Copy the necessary files into the work directory
		shutil.copy(input.gro_cg, os.path.join(work_dir, os.path.basename(input.gro_cg)))
		shutil.copy(input.top, os.path.join(work_dir, os.path.basename(input.top)))
		shutil.copy(input.mdp, os.path.join(work_dir, os.path.basename(input.mdp)))

		exec_dir = os.path.dirname(output.scheduling)
		# TPR path will be created inside exec_dir with the basename of output.tpr
		tpr_path = os.path.join(exec_dir, os.path.basename(output.tpr_file))
		# absolute path for grompp logs
		#log_abs = os.path.join(exec_dir, "grompp.log")

		# Compile the .tpr ONLY if it doesn't exist yet
		if not os.path.exists(tpr_path):
			shell("""
				cd {exec_dir}
				gmx grompp -f {params.mdp_abs} \
					-o $(basename {output.tpr_file}) \
					-c $(basename {input.gro_cg}) \
					-r $(basename {input.gro_cg}) \
					-p $(basename {input.top}) -maxwarn 1 > {params.log_abs} 2>&1
			""")

		# Clean up the temporary copies we made in exec_dir
		os.remove(os.path.join(exec_dir, os.path.basename(input.gro_cg)))
		os.remove(os.path.join(exec_dir, os.path.basename(input.top)))
		os.remove(os.path.join(exec_dir, os.path.basename(input.mdp)))
		
		if compute_target == "local": # write local scheduling into job description
			with open(output.scheduling, 'w') as f:
				f.write("COMPUTE: local")
		elif compute_target == "HPC":
			with open(output.scheduling, 'w') as f:
				f.write("COMPUTE: HPC")
		else:
			raise ValueError(f"{compute_target} Compute type not permitted ")

# create molecular dynamics job, which can be submitted komondor HPC later, or run locally. 
rule create_md_job:
	input:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/scheduling.yml"
	output:
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
	log:
		"logs/{pdb}/{source}/{model_id}/create_md_job.log"
	params:
		#log_abs = lambda wildcards, log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/create_md_job.log"
		)
	run:
		# 1. Setup Isolated Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger = logging.getLogger(f"create_md_job_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}")
		logger.setLevel(logging.INFO)
		logger.handlers.clear()  # Prevent duplicate handlers on re-runs

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

		# File Handler
		fh = logging.FileHandler(log_path, mode="w")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		# Console Handler
		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}"
		logger.info("Initializing MD job creation for target: %s", target_id)

		try:
			# 2. Parse Scheduling Data
			if not os.path.exists(input.scheduling):
				raise FileNotFoundError(f"Scheduling file not found at {input.scheduling}")

			with open(input.scheduling, "r") as f:
				scheduling_info = yaml.load(f, Loader=yaml.FullLoader) or {}

			compute_target = scheduling_info.get("COMPUTE")
			logger.info("Parsed compute target: '%s'", compute_target)

			# 3. Handle Local Target
			if compute_target == "local":
				logger.info("Local execution target confirmed for %s. Writing placeholder job script.", target_id)
				os.makedirs(os.path.dirname(output.job_description), exist_ok=True)
				with open(output.job_description, "w") as f:
					f.write("BLANK\n")
				logger.info("Placeholder job script created: %s", output.job_description)

			# 4. Handle HPC Target
			elif compute_target == "HPC":
				submit_hpc = str(SUBMIT_HPC).strip() == "1"

				if not submit_hpc:
					logger.error("HPC submission is disabled in environment (SUBMIT_HPC=%s).", SUBMIT_HPC)
					raise ValueError(f"HPC submission disabled for {target_id}. Set SUBMIT_HPC=1 in .env to enable.")

				job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
				script_path = os.path.abspath(output.job_description)
				deffnm = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"

				# Generate Slurm Batch Script
				slurm_script = f"""#! /bin/bash
# Komondor Slurm Template for Gromacs 2025.4 (Compiled: MPI noCUDA)
# Run: CPU, MPI
# KV

#SBATCH --job-name={job_name}
#SBATCH --output={deffnm}_%j.out
#SBATCH --error={deffnm}_%j.err
#SBATCH --account={HPC_ACCOUNT}
#SBATCH --partition={HPC_PARTITION}
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time={HPC_TIME}
#SBATCH --no-requeue
#SBATCH --exclusive

# Module load, export
export OMP_NUM_THREADS=16
#export GMX_GPU_DD_COMMS=true
#export GMX_GPU_PME_PP_COMMS=true
#export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# Module setup
module --force purge
module load PrgEnv-gnu/8.6.0
module load gcc/12.2.0
module load cray-mpich
module load cray-fftw

# set network and cray mpich env
export OMP_NUM_THREADS=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_PROVIDER=cxi
export MPICH_OFI_STARTUP_CONNECT=1

# fix shared  library pathing
# export CRAY_FFTW_DIR="${{CRAY_FFTW_DIR:-$FFTW_DIR}}"
# export LD_LIBRARY_PATH="${{CRAY_FFTW_DIR}}/lib:${{LD_LIBRARY_PATH}}"

# --- Source GROMACS Binaries ---
export gmxhome={HPC_GMX_HOME}
export PATH="${{gmxhome}}/bin:${{PATH}}"
export LD_LIBRARY_PATH="${{gmxhome}}/lib64:${{gmxhome}}/lib:${{LD_LIBRARY_PATH}}"




# --- Helper Function for MDRun ---
# This checks if a checkpoint exists for the specific step to resume it

# Helper Function for MDRun
run_md() {{
	local name="$1"
	if [ -f "${{name}}.cpt" ]; then
		echo "--> Resuming ${{name}} from checkpoint..."
		srun gmx_mpi mdrun -v -deffnm "${{name}}" -dlb yes -cpi "${{name}}.cpt"
	else
		echo "--> Starting ${{name}}..."
		srun gmx_mpi mdrun -v -deffnm "${{name}}" -dlb yes
	fi
}}

# Production MD Execution
if [ ! -f "{deffnm}.gro" ]; then
	run_md "{deffnm}"
fi
"""
				os.makedirs(os.path.dirname(script_path), exist_ok=True)
				with open(script_path, "w") as fh:
					fh.write(slurm_script)
				logger.info("Generated Slurm batch script: %s", script_path)

				# Create Tarball Archive
				job_dir = os.path.dirname(output.job_description)
				job_targz = f"{job_dir}.tar.gz"
				logger.info("Archiving directory '%s' into '%s'...", job_dir, job_targz)

				with tarfile.open(job_targz, "w:gz") as tar:
					for fn in os.listdir(job_dir):
						p = os.path.join(job_dir, fn)
						tar.add(p, arcname=os.path.basename(fn))
				logger.info("Archive created successfully (Size: %.2f KB)", os.path.getsize(job_targz) / 1024.0)

				# Remote Sync via SSH/Rsync
				remote_dir = os.path.join(HPC_REMOTE_BASE, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns")
				ssh_target = f"{HPC_USER}@{HPC_HOST}"
				logger.info("Deploying archive to remote target %s:%s", ssh_target, remote_dir)

				shell("""
					ssh {ssh_target} "mkdir -p '{remote_dir}'"
					rsync -avz "{job_targz}" "{ssh_target}:{remote_dir}/JOB.tar.gz"
					ssh {ssh_target} "test -f '{remote_dir}/JOB.tar.gz' && echo 'Remote payload verified at {remote_dir}/JOB.tar.gz'"
				""")
				logger.info("Remote transfer and payload verification completed.")

				# Append Submission Metadata
				with open(input.scheduling, "a") as f:
					f.write("JOB_STATUS: Submitted\n")
					f.write(f"REMOTE_DIR: {remote_dir}\n")
				logger.info("Updated scheduling file '%s' with submission metadata.", input.scheduling)

			else:
				logger.error("Invalid COMPUTE target '%s' in %s", compute_target, input.scheduling)
				raise ValueError(f"Unknown compute target in scheduling info: {compute_target}")

			logger.info("Rule create_md_job completed successfully for %s.", target_id)

		except Exception as err:
			logger.exception("Execution failed in create_md_job for %s: %s", target_id, str(err))
			raise


# Rule 5: Run 100n Molecular Dynamics

rule run_molecular_dynamics:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/production_mdrun.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/production_mdrun.log"
		)
	resources:
		gpu = 1
	run:
		# Setup Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger = logging.getLogger(f"run_md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}")
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
		fh = logging.FileHandler(log_path, mode="a")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
		job_dir = os.path.dirname(input.job_description)
		local_dir = os.path.abspath(os.path.dirname(output.done))
		step_log_file = os.path.join(local_dir, "simulation_steps.log")

		def log_sim_step(step_name, details=""):
			"""Writes timestamped step tracking events to simulation_steps.log."""
			timestamp = subprocess.check_output('date +"%Y-%m-%d %H:%M:%S"', shell=True).decode().strip()
			line = f"[{timestamp}] [STEP: {step_name}] {details}\n"
			with open(step_log_file, "a") as f:
				f.write(line)
		
		def get_remote_gromacs_progress(
			ssh_target: str, remote_dir: str, prefix: str
		) -> str | None:
			"""Fetches real-time ETA or step count directly from Slurm .err or gmx log."""
			# First, look for live stdout/stderr dumps matching job error log patterns
			remote_cmd = (
				f"cd {remote_dir} && "
				f"if [ -f *.err ]; then tail -n 20 *.err | grep -i 'will finish'; "
				f"elif [ -f {prefix}.log ]; then tail -n 100 {prefix}.log | grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+\\.'; "
				f"fi"
			)
		
			res = subprocess.run(
				["ssh", ssh_target, remote_cmd],
				capture_output=True,
				text=True,
				check=False,
			)
		
			if res.returncode != 0 or not res.stdout.strip():
				return None
		
			last_line = res.stdout.strip().splitlines()[-1]
		
			# If captured from .err file (e.g. "step 6144100, will finish Wed Aug 26 04:04:18 2026")
			if "will finish" in last_line:
				# Extract the step and ETA cleanly
				return last_line.strip()
		
			# Fallback for parsing step/ps from .log
			parts = last_line.split()
			try:
				ps_val = float(parts[1])
				return f"{ps_val / 1000.0:.2f} ns"
			except (IndexError, ValueError):
				return None

		logger.info("Starting MD execution phase for target: %s", target_id)
		
		# Preserve first start time if log exists, or mark INIT
		if not os.path.exists(step_log_file):
			log_sim_step("INIT", f"Target: {target_id}")
		else:
			log_sim_step("RESUME_WORKFLOW", f"Workflow checked target: {target_id}")

		try:
			with open(input.job_description, "r") as f:
				job_description = f.read().strip()

			# -------------------------------------------------------------
			# LOCAL EXECUTION PATHWAY
			# -------------------------------------------------------------
			if job_description == "BLANK":
				logger.info("Target configured for LOCAL execution.")
				log_sim_step("EXEC_MODE", "Local execution requested")
				cpt_path = os.path.join(job_dir, f"{prefix}.cpt")

				if os.path.exists(cpt_path):
					logger.info("Active checkpoint detected. Resuming local MD run...")
					log_sim_step("MD_RESUME", f"Resuming from checkpoint {prefix}.cpt")
					shell("""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-cpi $(basename {cpt_path}) \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")
				else:
					logger.info("Launching fresh local MD run...")
					log_sim_step("MD_START", f"Starting fresh mdrun for {prefix}")
					shell("""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")

				shell("""
					cd "{job_dir}"
					[ -f "{prefix}.xtc" ] && ln -sf "{prefix}.xtc" ../md.xtc
					[ -f "{prefix}.tpr" ] && ln -sf "{prefix}.tpr" ../md.tpr
				""")
				log_sim_step("MD_COMPLETE", "Local run finished successfully")

			# HPC EXECUTION PATHWAY
			# -------------------------------------------------------------
			else:
				ssh_target = f"{HPC_USER}@{HPC_HOST}"
				clean_base = HPC_REMOTE_BASE.lstrip("~/")
				remote_dir = os.path.join(
					clean_base,
					f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns",
				)
				job_name = (
					f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
				)
				job_script = os.path.basename(input.job_description)

				log_sim_step("EXEC_MODE", f"HPC execution target ({ssh_target})")

				# 1. Check remote completion upfront
				check_cmd = [
					"ssh",
					ssh_target,
					f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES || echo NO",
				]
				remote_done = (
					subprocess.check_output(check_cmd, text=True).strip()
				)

				if remote_done == "YES":
					logger.info(
						"🎉 Remote simulation already finished! (Found %s.gro"
						" on HPC)",
						prefix,
					)
					log_sim_step(
						"CHECK_REMOTE",
						f"Simulation finished ({prefix}.gro present)",
					)
				else:
					# 2. Check if job is active in queue (Returns Job ID and State)
					squeue_cmd = [
						"ssh",
						ssh_target,
						f"squeue -u {HPC_USER} -n {job_name} -h -o '%i %t'",
					]
					q_res = subprocess.run(
						squeue_cmd, capture_output=True, text=True, check=False
					)
					q_out = q_res.stdout.strip()

					job_id = None
					if q_out:
						# Extract Job ID if already running/pending
						job_id = q_out.split()[0]
						logger.info(
							"⏳ Slurm Job '%s' already active in queue (Job"
							" ID: %s).",
							job_name,
							job_id,
						)
						log_sim_step(
							"SLURM_ACTIVE",
							f"Job ID {job_id} running/queued on HPC",
						)
					else:
						logger.info(
							"🚀 Submitting new Slurm job to Komondor (%s)...",
							ssh_target,
						)
						submit_cmd = [
							"ssh",
							ssh_target,
							f"cd '{remote_dir}' && [ -f JOB.tar.gz ] && tar"
							f" -xzf JOB.tar.gz && sbatch {job_script}",
						]
						sub_res = subprocess.run(
							submit_cmd,
							capture_output=True,
							text=True,
							check=False,
						)

						if (
							sub_res.returncode == 0
							and "Submitted batch job" in sub_res.stdout
						):
							job_id = sub_res.stdout.strip().split()[-1]
							logger.info(
								"✅ Slurm job successfully submitted!"
								" Assigned Job ID: %s",
								job_id,
							)
							log_sim_step(
								"SLURM_SUBMIT",
								f"Submitted batch job ID: {job_id}",
							)
						else:
							logger.error(
								"❌ Failed to submit Slurm job: %s",
								sub_res.stderr,
							)
							log_sim_step(
								"SLURM_FAILED",
								f"Submission failed: {sub_res.stderr.strip()}",
							)
							raise RuntimeError(
								f"Slurm sbatch submission failed: {sub_res.stderr}"
							)

					# 3. Monitor Slurm Execution Loop
					if job_id:
						logger.info("Monitoring Slurm Job %s...", job_id)
						last_progress_msg = ""

						while True:
							# Query job state safely (check=False avoids CalledProcessError when job finishes)
							if job_id and str(job_id).isdigit():
								check_q = [
									"ssh",
									ssh_target,
									f"squeue -j {job_id} -h -o '%t'",
								]
								q_check = subprocess.run(
									check_q,
									capture_output=True,
									text=True,
									check=False,
								)
							
								# Clean stdout: take only the first token (e.g. "R" or "PD")
								stdout_clean = q_check.stdout.strip()
								job_state = stdout_clean.split()[0] if stdout_clean else ""
							
								if not job_state:
									logger.info(
										"Job %s left the queue. Verifying completion...",
										job_id,
									)
									break  # Job completed or died; exit polling loop
							
								if job_state == "R":
									progress = get_remote_gromacs_progress(ssh_target, remote_dir, prefix)
									if progress and progress != last_progress_msg:
										log_sim_step(
											"HEARTBEAT",
											f"Job {job_id} running - {progress}",
										)
										last_progress_msg = progress
									else:
										log_sim_step(
											"HEARTBEAT",
											f"Job {job_id} actively executing on HPC",
										)
								else:
									log_sim_step(
										"HEARTBEAT",
										f"Job {job_id} queued (State: {job_state})",
									)

							time.sleep(180)  # Poll every 3 minutes

						# 4. Final Verification: Confirm simulation produced output
						post_check = [
							"ssh",
							ssh_target,
							f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES"
							" || echo NO",
						]
						final_done = (
							subprocess.check_output(post_check, text=True)
							.strip()
						)

						if final_done == "YES":
							log_sim_step(
								"SLURM_FINISHED",
								f"Slurm Job {job_id} completed successfully",
							)
						else:
							log_sim_step(
								"SLURM_FAILED",
								f"Slurm Job {job_id} ended without producing"
								f" {prefix}.gro",
							)
							raise RuntimeError(
								f"HPC Job {job_id} terminated unexpectedly"
								f" ({prefix}.gro missing)."
							)

				# Sync and retrieve results locally
				logger.info("📦 Archiving and retrieving simulation artifacts from HPC...")
				log_sim_step("RETRIEVE_START", "Fetching remote output files via rsync")



				shell(f"""
					SSH_TARGET="{ssh_target}"
					REMOTE_DIR="{remote_dir}"
					PREFIX="{prefix}"
					LOCAL_DIR="{local_dir}"
					LOG_FILE="{log_path}"

					ssh "$SSH_TARGET" "cd '$REMOTE_DIR' && tar -czf md_results.tar.gz $PREFIX.*" >> "$LOG_FILE" 2>&1
					rsync -avz "$SSH_TARGET:$REMOTE_DIR/md_results.tar.gz" "$LOCAL_DIR/" >> "$LOG_FILE" 2>&1
					mkdir "$LOCAL_DIR/md_results"
					tar -xzf "$LOCAL_DIR/md_results.tar.gz" -C "$LOCAL_DIR/md_results"

					# I don't want to bother with this yet.
					#[ -f "$LOCAL_DIR/$PREFIX.xtc" ] && ln -sf "$PREFIX.xtc" "$LOCAL_DIR/md.xtc"
					#[ -f "$LOCAL_DIR/$PREFIX.tpr" ] && ln -sf "$PREFIX.tpr" "$LOCAL_DIR/md.tpr"

					rm -f "$LOCAL_DIR/md_results.tar.gz"
					ssh "$SSH_TARGET" "rm -f '$REMOTE_DIR/md_results.tar.gz'" >> "$LOG_FILE" 2>&1
				""")

				log_sim_step("RETRIEVE_COMPLETE", "Downloaded, extracted, and cleaned up md_results.tar.gz")

			# Sentinel output
			with open(output.done, "w") as f:
				f.write(f"MD simulation completed for {target_id}.\n")
			
			log_sim_step("DONE", "Pipeline rule finished successfully.")

		except Exception as err:
			log_sim_step("ERROR", str(err))
			logger.exception("Execution failed in run_molecular_dynamics for %s: %s", target_id, str(err))
			raise

#rule md_repackaging:



# Rule 6: Process Trajectory and Extract Low-Energy/Representative Snapshots
# rule process_trajectory:
#     input:
#         xtc = "results/gromacs/{pdb}/md/trajectory.xtc",
#         tpr = "results/gromacs/{pdb}/md/sim.tpr",
#     output:
#         snapshots = directory("results/{pdb}/snapshots/")
#     script:
#         "../scripts/cluster_trajectory.py"

# --- STEP 6: TRAJECTORY CLEANING & PBC WRAPPING CORRECTION ---
rule pbc_correction_and_extract:
	input:
		md_checkpoint_guard = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt",
		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md.xtc",
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md.tpr",
	output:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/FRAMES_compressed.tar.gz"
	log:
		"logs/{pdb}/{source}/{model_id}/trjconv_pbc.log"
	shell:
		"""
		# 1. Recenter molecular unity boundaries
		echo "Protein" | gmx trjconv -f {input.xtc_md} -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc -pbc mol -ur compact > {log} 2>&1
		
		# 2. Fit rotational and translational structural drift
		echo "Protein Protein Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc -center -fit rot+trans >> {log} 2>&1
		
		# 3. Chop trajectory into individual PDB frames
		mkdir -p results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/FRAMES
		echo "Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/FRAMES/frame.pdb -s {input.tpr_md} -sep >> {log} 2>&1
		
		# Compress and archive individual coordinate files
		tar -czf {output.tar} -C results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id} FRAMES >> {log} 2>&1
		
		# Clean up massive intermediate trajectories to preserve space
		# rm -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc
		"""


# --- STEP 7: AUTOMATED PYMOL MOVIE GENERATION ---
rule generate_pymol_movie:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/FRAMES_compressed.tar.gz"
	output:
		movie = "results/gromacs/{pdb}/{source}/{model_id}/movies/100ns_md_trajectory.mov"
	log:
		"logs/{pdb}/{source}/{model_id}/pymol_render.log"
	shell:
		"python scripts/render_pymol_movie.py {input.tar} {output.movie} > {log} 2>&1"
		
# Rule 4: Run Quantum Mechanical / Excited-State Calculations on Snapshots
# rule run_quantum_mechanics:
	# input:
		# snapshots = "results/{pdb}/snapshots/"
	# output:
		# qm_out = "results/{pdb}/qm_results.dat"
	# shell:
		# # Loops through extracted snapshots and runs MOPAC or ORCA
		# """
		# for f in {input.snapshots}/*.inp; do
			# {config[mopac_command]} $f
		# done
		# touch {output.qm_out}
		# """

# # Rule 5: Compile calculations and plot final theoretical UV-Vis or Enthalpy graph
# rule plot_results:
	# input:
		# qm_out = "results/{pdb}/qm_results.dat"
	# output:
		# plot = "results/{pdb}/final_spectra.png"
	# script:
		# "../scripts/parse_qm_spectra.py"
