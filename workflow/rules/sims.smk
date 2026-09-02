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
from pathlib import Path

wildcard_constraints:
	pdb = "[^/]+",
	source = "[^/]+",
	model_id = "[^/]+",
	protocol = "[^/]+"
	
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

def resolve_slurm_job_script(wildcards, protocol):

	deffnm = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
	#job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{protocol}"
	job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
	
	if protocol=='standard_100ns':
		return f"""#! /bin/bash
# Komondor Slurm Template for Gromacs 2025.4 (Compiled: MPI noCUDA)
# Run: CPU, MPI
# KV

#SBATCH --job-name={job_name}
#SBATCH --output={deffnm}_%j.out
#SBATCH --error={deffnm}_%j.err
#SBATCH --account={HPC_ACCOUNT}
#SBATCH --partition={HPC_PARTITION}
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --time={HPC_TIME}
#SBATCH --no-requeue
#SBATCH --exclusive

# Module setup
module --force purge
module load PrgEnv-gnu/8.6.0
module load gcc/12.2.0
module load cray-mpich
module load cray-fftw

export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# set network and cray mpich env
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
export MPICH_SMP_SINGLE_COPY_MODE=CMA
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
		srun --nodes=2 \
			--ntasks-per-node=16 \
			--cpus-per-task=8 \
			--cpu-bind=cores \ 
			gmx_mpi mdrun -v -deffnm "${{name}}" \
				-dds 0.8 \
				-ntomp 8 \
				-rcon 0 \
				-dlb yes \
				-cpi "${{name}}.cpt"
	else
		echo "--> Starting ${{name}}..."
		srun --nodes=2 \
			--ntasks-per-node=16 \
			--cpus-per-task=8 \
			--cpu-bind=cores \
			gmx_mpi mdrun -v -deffnm "${{name}}" \
				-ntomp 8 \
				-dds 0.8 \
				-rcon 0 \
				-dlb yes 
	fi
}}

# Production MD Execution
if [ ! -f "{deffnm}.gro" ]; then
	run_md "{deffnm}"
fi
"""
	elif protocol=='extended_1000ns':

		return f"""#! /bin/bash
		ns=900
		ps_to_ext=$(( ns * 1000 ))

		gmx convert-tpr -s md.tpr -extend $ps_to_ext -o md_7KPH_extended.tpr

		rm md.tpr

		gmx mdrun -s md_7KPH_extended.tpr -cpi state.cpt -deffnm md # -noappend
		"""
	elif protocol=='simulated_annealation':
		
		return f"""#! /bin/bash
		#!/bin/bash

		# --- CONFIGURATION VARIABLES ---
		PDBCODE="7KPH"
		SOURCE_DIR="$HOME/compchem_research/workspace/data/gromacs/runs/100ns/staged/v2/7KPH/empirical/7KPH"

		echo "================================================================="
		echo " GROMACS Simulated Annealing Execution Wrapper"
		echo "================================================================="

		# --- CHECK IF THIS IS A RESUMPTION OR A NEW RUN ---
		# If md.cpt exists, an active annealing run was interrupted and we resume it.
		if [ -f "md.cpt" ]; then
			echo "[!] Detected existing progress checkpoint (md.cpt)."
			echo "[->] Continuing the interrupted Simulated Annealing simulation..."

			gmx mdrun -s md_${PDBCODE}_anneal.tpr -cpi md.cpt -deffnm md

		else
			echo "[+] No active annealing run found. Initializing fresh protocol..."
			echo "[+] Copying baseline structure and topologies from archive..."

			# Copy the baseline structure template and the pristine 100ns checkpoint mark
			cp "${SOURCE_DIR}/md.tpr" ./
			cp "${SOURCE_DIR}/md.cpt" ./state.cpt

			# CRUCIAL: grompp needs your system topology layout to build a new physics matrix.
			# We copy the topol.top and any inclusion files (.itp) from your source directory.
			cp "${SOURCE_DIR}/topol.top" ./ 2>/dev/null || echo "[!] Check: Verify topol.top is in your source directory."
			cp "${SOURCE_DIR}"/*.itp ./ 2>/dev/null

			echo "[+] Compiling new Simulated Annealing topology via gmx grompp..."
			# We use your custom 'anneal.mdp'. By passing '-t state.cpt', grompp extracts
			# the exact coordinates and velocities from the end of your 100ns run.
			gmx grompp -f anneal.mdp -c md.tpr -t state.cpt -p topol.top -o md_${PDBCODE}_anneal.tpr

			# Clean up the original reference .tpr to keep the directory tidy
			rm md.tpr

			echo "[->] Launching Simulated Annealing production run..."
			gmx mdrun -s md_${PDBCODE}_anneal.tpr -deffnm md
		fi

		echo "================================================================="
		echo " Execution block completed."
		echo "================================================================="
		"""
	else:
		raise ValueError(f"{protocol} not recognized as an approved protocol for slurm job")


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
#	shell:
#		'python scripts/sims/schedule_md_job.py \
#			--in_top=input.top \
#			--in_gro=input.gro_cg  \
#			--in_mdp=input.mdp \
#			--out_scheduling=output.scheduling' \
#			--out_tpr=output.tpr_file \
#			--pdb=wildcards.pdb \
#			--source=wildcards.source \
#			--model_id=wildcards.model_id \
#			--abs_log=params.log_abs
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
		
		scheduling_info = {}

		if compute_target in ["local", 'HPC']: # write local scheduling into job description
			scheduling_info['COMPUTE'] = compute_target
			with open (output.scheduling, 'w') as f:
				yaml.safe_dump(scheduling_info, f, sort_keys=False)
			#with open(output.scheduling, 'w') as f:
			#	f.write("COMPUTE: local")
		#elif compute_target == "HPC":
			#with open(output.scheduling, 'w') as f:
			#	f.write("COMPUTE: HPC")
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

			#with open(input.scheduling, "r") as f:
			#	scheduling_info = yaml.load(f, Loader=yaml.FullLoader) or {}
			with open(input.scheduling, "r") as f:
				scheduling_info = yaml.safe_load(f) or {}
			

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

				script_path = os.path.abspath(output.job_description)
				#deffnm = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
				#job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"

				# Generate Slurm Batch Script

##SBATCH --nodes=4
##SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
				slurm_script = resolve_slurm_job_script(wildcards, protocol='standard_100ns')

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
				#ssh_target = f"{HPC_USER}@{HPC_HOST}"
				ssh_target = "komondor"
				logger.info("Deploying archive to remote target %s:%s", ssh_target, remote_dir)

				shell("""
					ssh -o BatchMode=yes {ssh_target} "mkdir -p '{remote_dir}'"
					rsync -e "ssh -o BatchMode=yes" -avz "{job_targz}" "{ssh_target}:{remote_dir}/JOB.tar.gz"
					ssh -o BatchMode=yes {ssh_target} "test -f '{remote_dir}/JOB.tar.gz' && echo 'Remote payload verified at {remote_dir}/JOB.tar.gz'"
				""")
				logger.info("Remote transfer and payload verification completed.")

				# Append Submission Metadata
				#with open(input.scheduling, "a") as f:
				#	f.write("JOB_STATUS: Submitted\n")
				#	f.write(f"REMOTE_DIR: {remote_dir}\n")
				scheduling_info['JOB_STATUS'] = 'Submitted'
				scheduling_info['REMOTE_DIR'] = str(remote_dir)

				with open (input.scheduling, 'w') as f:
					yaml.safe_dump(scheduling_info, f, sort_keys=False)
				
				logger.info("Updated scheduling file '%s' with submission metadata.", input.scheduling)

			else:
				logger.error("Invalid COMPUTE target '%s' in %s", compute_target, input.scheduling)
				raise ValueError(f"Unknown compute target in scheduling info: {compute_target}")

			logger.info("Rule create_md_job completed successfully for %s.", target_id)

		except Exception as err:
			logger.exception("Execution failed in create_md_job for %s: %s", target_id, str(err))
			raise

# NOTE: Rule BRANCHING
# split between LOCAL and HPC compute
# Rule 5: Run 100n Molecular Dynamics

def get_compute_target(wildcards):
	protocol = getattr(wildcards, "protocol", "standard_100ns")

	def_compute_target = config.get("default_compute_target", "local")
	# get target from config for this PDB, or use default
	pdb_config_entries = config.get("custom_simulations", {})\
					.get(wildcards.pdb, {})\
					.get(f'{wildcards.source}_{wildcards.model_id}', [])

	current_protocol = "standard"
	pdb_config = next((entry for entry in pdb_config_entries if entry.get("protocol") == current_protocol), {})

	compute_target = pdb_config.get("compute_target", def_compute_target)
	if compute_target in ["HPC", "local"]:
		return compute_target

	scheduling_file = Path(
		f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{protocol}/JOB/scheduling.yml"
	)

	if scheduling_file.exists():
		try:
			with open(scheduling_file, 'r') as f:
				data = yaml.safe_load(f) or {}
				# Extract 'COMPUTE' key directly from scheduling.yml
				target = data.get("COMPUTE")
				if target in ["HPC", "local"]:
					return target
		except Exception as e:
			print(f"Error reading {scheduling_file}: {e}")

	# Fallback to config default if scheduling.yml is missing or unreadable
	return config.get("default_compute_target", "local")


def det_compute_scheduling(wildcards):
	protocol = getattr(wildcards, "protocol", "standard_100ns")

	mode = get_compute_target(wildcards)
	if mode in ['local', 'HPC']:
		return f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{protocol}/{mode}_md_completed.txt"
	else:
		raise ValueError(f"Compute type '{mode}' not permitted.")

def log_sim_step(step_name, step_log_file, details=""):
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

rule run_local_md:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/local_md_completed.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/local_production_mdrun.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/local_production_mdrun.log"
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

		target_mode = get_compute_target(wildcards)

		if target_mode != "local":
			logger.info("Target %s is configured for HPC mode. Skipping local rule.", log)



		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
		fh = logging.FileHandler(log_path, mode="a")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
		job_dir = os.path.abspath(os.path.dirname(input.job_description))
		local_dir = os.path.abspath(os.path.dirname(output.done))
		step_log_file = os.path.join(job_dir, "simulation_steps.log")

		logger.info("Starting MD execution phase for target: %s", target_id)

		if not os.path.exists(step_log_file):
			log_sim_step("INIT", step_log_file, f"Target: {target_id}")
		else:
			log_sim_step("RESUME_WORKFLOW", step_log_file, f"Workflow checked target: {target_id}")

		try:
			with open(input.job_description, "r") as f:
				job_description = f.read().strip()

			if job_description == "BLANK":
				logger.info("Target configured for LOCAL execution.")
				log_sim_step("EXEC_MODE", step_log_file, "Local execution requested")
				cpt_path = os.path.join(job_dir, f"{prefix}.cpt")

				if os.path.exists(cpt_path):
					logger.info("Active checkpoint detected. Resuming local MD run...")
					log_sim_step("MD_RESUME", step_log_file, f"Resuming from checkpoint {prefix}.cpt")
					shell(f"""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-cpi {prefix}.cpt \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")
				else:
					logger.info("Launching fresh local MD run...")
					log_sim_step("MD_START", step_log_file, f"Starting fresh mdrun for {prefix}")
					shell(f"""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")

				# Ensure md_results directory exists before linking/copying
				#md_results_dir = os.path.abspath(os.path.join(local_dir, "..", "md_results"))
				#md_results_dir = local_dir
				#os.makedirs(md_results_dir, exist_ok=True)

				shell(f"""
					cd "{job_dir}"
					[ -f "{prefix}.xtc" ] && cp -f "{prefix}.xtc" "{local_dir}/{prefix}.xtc" || true
					[ -f "{prefix}.tpr" ] && cp -f "{prefix}.tpr" "{local_dir}/{prefix}.tpr" || true
				""")

				# Write sentinel file
				with open(output.done, "w") as f:
					f.write(f"Local MD completed for {target_id}.\n")

				log_sim_step("MD_COMPLETE", step_log_file, "Local run finished successfully")
			else:
				raise ValueError(f"Job description is '{job_description}', expected 'BLANK' for local runs")

		except Exception as err:
			log_sim_step("ERROR", step_log_file, str(err))
			logger.exception("Execution failed in run_molecular_dynamics for %s: %s", target_id, str(err))
			raise

rule run_HPC_md:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
	output:
		#md_dir = directory("results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results"),
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/HPC_md_completed.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/HPC_production_mdrun.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/HPC_production_mdrun.log"
		)
	resources:
		gpu = 0,     # Zero local GPUs used! Passive SSH / monitoring thread only
		mem_mb = 500
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
		step_log_file = os.path.join(job_dir, "simulation_steps.log")

		logger.info("Starting MD execution phase for target: %s", target_id)

		if not os.path.exists(step_log_file):
			log_sim_step("INIT", step_log_file, f"Target: {target_id}")
		else:
			log_sim_step("RESUME_WORKFLOW", step_log_file, f"Workflow checked target: {target_id}")

		try:
			with open(input.job_description, "r") as f:
				job_description = f.read().strip()

			ssh_target = "komondor"
			clean_base = HPC_REMOTE_BASE.lstrip("~/")
			remote_dir = os.path.join(
				clean_base,
				f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns",
			)
			job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
			job_script = os.path.basename(input.job_description)

			log_sim_step("EXEC_MODE", step_log_file, f"HPC execution target ({ssh_target})")

			# Guarantee remote directory exists before testing files
			subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, f"mkdir -p '{remote_dir}'"], check=False)

			# 1. Check remote completion upfront
			check_cmd = [
				"ssh",
				"-o", "BatchMode=yes",
				ssh_target,
				f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES || echo NO",
			]
			check_res = subprocess.run(check_cmd, capture_output=True, text=True, check=False)

			if check_res.returncode == 0 and "YES" in check_res.stdout:
				remote_done = "YES"
			else:
				remote_done = "NO"

			if remote_done == "YES":
				logger.info("🎉 Remote simulation already finished! (Found %s.gro on HPC)", prefix)
				log_sim_step("CHECK_REMOTE", step_log_file, f"Simulation finished ({prefix}.gro present)")
			else:
				# 2. Check if job is active in queue (Returns Job ID if running or pending)
				squeue_cmd = [
					"ssh",
					"-o", "BatchMode=yes",
					ssh_target,
					f"squeue -u {HPC_USER} -n {job_name} -h -o '%i %t'",
				]
				q_res = subprocess.run(squeue_cmd, capture_output=True, text=True, check=False)
				q_out = q_res.stdout.strip()

				job_id = None
				if q_out:
					# Job is ALREADY active! Extract Job ID and attach monitor loop to it.
					job_id = q_out.split()[0]
					logger.info("⏳ Slurm Job '%s' already active in queue (Job ID: %s). Attaching monitor loop...", job_name, job_id)
					log_sim_step("SLURM_ACTIVE", step_log_file, f"Job ID {job_id} running/queued on HPC")
				else:
					# Job is not in queue; submit a new Slurm batch job
					logger.info("🚀 Submitting new Slurm job to Komondor %s @(%s)...", job_name, ssh_target)
					submit_cmd = [
						"ssh",
						"-o", "BatchMode=yes",
						ssh_target,
						f"cd '{remote_dir}' && [ -f JOB.tar.gz ] && tar -xzf JOB.tar.gz; sbatch {job_script}",
					]
					sub_res = subprocess.run(submit_cmd, capture_output=True, text=True, check=False)

					if sub_res.returncode == 0 and "Submitted batch job" in sub_res.stdout:
						job_id = sub_res.stdout.strip().split()[-1]
						logger.info("✅ Slurm job successfully submitted! Assigned Job ID: %s", job_id)
						log_sim_step("SLURM_SUBMIT", step_log_file, f"Submitted batch job ID: {job_id}")
					else:
						logger.error("❌ Failed to submit Slurm job: %s", sub_res.stderr)
						log_sim_step("SLURM_FAILED", step_log_file, f"Submission failed: {sub_res.stderr.strip()}")
						raise RuntimeError(f"Slurm sbatch submission failed: {sub_res.stderr}")

				# 3. Monitor Slurm Execution Loop (EVALUATED FOR BOTH ACTIVE & NEW JOBS)
				if job_id:
					logger.info("Monitoring Slurm Job %s...", job_id)
					last_progress_msg = ""

					while True:
						if job_id and str(job_id).isdigit():
							check_q = [
								"ssh",
								"-o", "BatchMode=yes",
								ssh_target,
								f"squeue -j {job_id} -h -o '%t'",
							]
							q_check = subprocess.run(check_q, capture_output=True, text=True, check=False)
							stdout_clean = q_check.stdout.strip()
							job_state = stdout_clean.split()[0] if stdout_clean else ""

							if not job_state:
								logger.info("Job %s left the queue. Verifying completion...", job_id)
								break  # Job completed or left queue; exit polling loop

							if job_state == "R":
								progress = get_remote_gromacs_progress(ssh_target, remote_dir, prefix)
								if progress and progress != last_progress_msg:
									log_sim_step("HEARTBEAT", step_log_file, f"Job {job_id} running - {progress}")
									last_progress_msg = progress
								else:
									log_sim_step("HEARTBEAT", step_log_file, f"Job {job_id} actively executing on HPC")
							else:
								log_sim_step("HEARTBEAT", step_log_file, f"Job {job_id} queued (State: {job_state})")

						time.sleep(900)  # Poll every 15 minutes

					# 4. Final Verification: Confirm simulation produced output
					post_check = [
						"ssh",
						"-o", "BatchMode=yes",
						ssh_target,
						f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES || echo NO",
					]
					post_res = subprocess.run(post_check, capture_output=True, text=True, check=False)
					final_done = post_res.stdout.strip() if post_res.returncode == 0 else "NO"

					if final_done == "YES":
						log_sim_step("SLURM_FINISHED", step_log_file, f"Slurm Job {job_id} completed successfully")
					else:
						log_sim_step("SLURM_FAILED", step_log_file, f"Slurm Job {job_id} ended without producing {prefix}.gro")
						raise RuntimeError(f"HPC Job {job_id} terminated unexpectedly ({prefix}.gro missing).")

			# Sync and retrieve results locally
			logger.info("📦 Archiving and retrieving simulation artifacts from HPC...")
			log_sim_step("RETRIEVE_START", step_log_file, "Fetching remote output files via rsync")

			shell(f"""
				SSH_TARGET="{ssh_target}"
				REMOTE_DIR="{remote_dir}"
				LOCAL_DIR="{local_dir}"
				LOG_FILE="{log_path}"

				mkdir -p "$LOCAL_DIR/md_results"

				# Directly pull all simulation artifacts, excluding job scripts/archives
				rsync -e "ssh -o BatchMode=yes" -avz \
					--exclude="JOB" \
					--exclude="JOB.tar.gz" \
					"$SSH_TARGET:$REMOTE_DIR/" \
					"$LOCAL_DIR/" >> "$LOG_FILE" 2>&1
			""")

#			shell(f"""
#				SSH_TARGET="{ssh_target}"
#				REMOTE_DIR="{remote_dir}"
#				PREFIX="{prefix}"
#				LOCAL_DIR="{local_dir}"
#				LOG_FILE="{log_path}"
#
#				# 1. Create remote archive using Python's {prefix} variable directly
#				ssh -o BatchMode=yes "$SSH_TARGET" "cd '$REMOTE_DIR' && tar -czf md_results.tar.gz {prefix}* *.tpr *.xtc 2>/dev/null || tar -czf md_results.tar.gz {prefix}* 2>/dev/null || true" >> "$LOG_FILE" 2>&1
#
#				# 2. Fetch archive locally
#				rsync -avz "$SSH_TARGET:$REMOTE_DIR/md_results.tar.gz" "$LOCAL_DIR/" >> "$LOG_FILE" 2>&1
#
#				# 3. Extract locally
#				mkdir -p "$LOCAL_DIR/md_results"
#				tar -xzf "$LOCAL_DIR/md_results.tar.gz" -C "$LOCAL_DIR/md_results"
#
#				# 4. Clean up temporary archives
#				rm -f "$LOCAL_DIR/md_results.tar.gz"
#				ssh -o BatchMode=yes "$SSH_TARGET" "rm -f '$REMOTE_DIR/md_results.tar.gz'" >> "$LOG_FILE" 2>&1
#			""")

			logger.info(f"MD DATA ACQUIRED CHECK @ {local_dir}/md_results")
			log_sim_step("RETRIEVE_COMPLETE", step_log_file, "Downloaded, extracted, and cleaned up md_results.tar.gz")

			# Sentinel output
			with open(output.done, "w") as f:
				f.write(f"MD simulation completed for {target_id}.\n")

			log_sim_step("DONE", step_log_file, "Pipeline rule finished successfully.")

		except Exception as err:
			log_sim_step("ERROR", step_log_file, str(err))
			logger.exception("Execution failed in run_molecular_dynamics for %s: %s", target_id, str(err))
			raise

rule finalize_md:
	input:
		det_compute_scheduling = det_compute_scheduling
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/md_completed.txt",
		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.xtc",
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.tpr",
	run:
		compute_target = get_compute_target(wildcards)
		out_dir = os.path.abspath(os.path.dirname(output.done))
		os.makedirs(out_dir, exist_ok=True)

		job_dir = os.path.abspath(f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/JOB")
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"

		# ------------------------------------------------------------------
		# COMMON FALLBACK: Ensure .tpr exists in md_results
		# (.tpr is created locally during grompp in JOB/, copy if remote tar missed it)
		# ------------------------------------------------------------------
		local_tpr_src = os.path.join(job_dir, f"{prefix}.tpr")
		if not os.path.exists(output.tpr_md) and os.path.exists(local_tpr_src):
			shutil.copy2(local_tpr_src, output.tpr_md)

		# ------------------------------------------------------------------
		# FALLBACK 2: Check for generic 'md.xtc' or misplaced files
		# ------------------------------------------------------------------
		if not os.path.exists(output.xtc_md):
			# Check if extracted as md.xtc inside md_results/
			generic_xtc = os.path.join(out_dir, "md.xtc")
			if os.path.exists(generic_xtc):
				os.rename(generic_xtc, output.xtc_md)

			# Check if extracted into parent standard_100ns directory
			parent_dir = os.path.abspath(os.path.join(out_dir, ".."))
			parent_xtc = os.path.join(parent_dir, f"{prefix}.xtc")
			if os.path.exists(parent_xtc):
				shutil.move(parent_xtc, output.xtc_md)
		# ------------------------------------------------------------------
		# HPC PATH
		# ------------------------------------------------------------------
		if compute_target == "HPC":
			# Check that required outputs exist before writing sentinel
			if not os.path.exists(output.xtc_md):
				raise FileNotFoundError(f"HPC execution finished, but {output.xtc_md} was not retrieved into md_results!")
			if not os.path.exists(output.tpr_md):
				raise FileNotFoundError(f"HPC execution finished, but {output.tpr_md} is missing from md_results!")

			with open(output.done, "w") as f:
				f.write(f"MD simulation (HPC) finalized for {wildcards.pdb}/{wildcards.source}/{wildcards.model_id}\n")
			return

		# ------------------------------------------------------------------
		# LOCAL PATH
		# ------------------------------------------------------------------
		if compute_target == "local":
			job_desc = os.path.join(job_dir, f"{prefix}_job.job")

			if not os.path.exists(job_desc):
				raise FileNotFoundError(f"Expected job description not found: {job_desc}")

			with open(job_desc, "r") as fh:
				content = fh.read().strip()

			if content != "BLANK":
				raise ValueError(f"Job description must be 'BLANK' for local finalize, got: '{content}'")

			# Collect all matching MD outputs from JOB into md_results
			if os.path.isdir(job_dir):
				for fn in os.listdir(job_dir):
					if fn.startswith(prefix):
						src = os.path.join(job_dir, fn)
						dst = os.path.join(out_dir, fn)
						try:
							shutil.copy2(src, dst)
						except Exception:
							pass

			# Validate that output files are present
			if not os.path.exists(output.xtc_md) or not os.path.exists(output.tpr_md):
				raise FileNotFoundError(f"Local MD completed, but missing expected outputs (.xtc / .tpr) in {out_dir}")

			with open(output.done, "w") as f:
				f.write(f"MD simulation (local) finalized and artifacts collected for {wildcards.pdb}/{wildcards.source}/{wildcards.model_id}\n")
			return

		raise ValueError(f"Unknown compute target when finalizing MD: {compute_target}")
		
#rule finalize_md:
#	input:
##		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/local_md_completed.txt"
##		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/HPC_md_completed.txt"
#		det_compute_scheduling = det_compute_scheduling
#	output:
#		#done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt"
#		#md_cpt = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.cpt"
#		#md_edr = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.edr"
#		#md_log = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.log"
#		#md_xtc = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.xtc"
#		#md_results = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results"
#		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/md_completed.txt",
#		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.xtc",
#		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.tpr",
#
#
#	run:
#		# Finalize MD: if HPC just write sentinel; if local, validate BLANK job and
#		# collect MD outputs into md_results then write sentinel file.
#		scheduling_path = f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/JOB/scheduling.yml"
#		compute_target = get_compute_target(wildcards)
#
#		# Ensure output directory exists
#		out_dir = os.path.abspath(os.path.dirname(output.done))
#		os.makedirs(out_dir, exist_ok=True)
#
#		# HPC path: nothing to fetch here (already retrieved in run_HPC_md), just write sentinel
#		if compute_target == "HPC":
#			with open(output.done, "w") as f:
#				f.write(f"MD simulation (HPC) finalized for {wildcards.pdb}/{wildcards.source}/{wildcards.model_id}\n")
#			return
#
#		# Local path: verify the JOB descriptor indicates local execution then collect files
#		if compute_target == "local":
#			job_dir = os.path.abspath(f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/JOB")
#			job_desc = os.path.join(job_dir, f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md_job.job")
#
#			if not os.path.exists(job_desc):
#				raise FileNotFoundError(f"Expected job description not found: {job_desc}")
#
#			with open(job_desc, "r") as fh:
#				content = fh.read().strip()
#
#			if content != "BLANK":
#				raise ValueError(f"Job description must be 'BLANK' for local finalize, got: '{content}'")
#
#			# Collect MD outputs into md_results
#			#md_results = os.path.join(out_dir, "md_results")
#			#os.makedirs(md_results, exist_ok=True)
#			md_results = out_dir
#
#			prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
#
#			# Copy files from JOB directory matching prefix
#			if os.path.isdir(job_dir):
#				for fn in os.listdir(job_dir):
#					if fn.startswith(prefix):
#						src = os.path.join(job_dir, fn)
#						dst = os.path.join(md_results, fn)
#						try:
#							shutil.copy2(src, dst)
#						except Exception:
#							# best-effort copy; continue on failure
#							pass
#
#			# Also copy common md artifacts in the parent standard_100ns directory
#			# dont do that that'd copy data from preparation
#			#parent_dir = os.path.abspath(os.path.dirname(job_dir))
#			#extra_candidates = [f"{prefix}.xtc", f"{prefix}.tpr", f"{prefix}.cpt", "md.xtc", "md.tpr", "md.cpt", "md.edr", "md.log"]
#			#for cand in extra_candidates:
#			#	src = os.path.join(parent_dir, cand)
#			#	if os.path.exists(src):
#			#		try:
#			#			shutil.copy2(src, os.path.join(md_results, os.path.basename(src)))
#			#		except Exception:
#			#			pass
#
#			# Write sentinel
#			with open(output.done, "w") as f:
#				f.write(f"MD simulation (local) finalized and artifacts collected for {wildcards.pdb}/{wildcards.source}/{wildcards.model_id}\n")
#			return
#
#		# Unknown compute target
#		raise ValueError(f"Unknown compute target when finalizing MD: {compute_target}")


# --- STEP 6: TRAJECTORY CLEANING & PBC WRAPPING CORRECTION ---
rule pbc_correction_and_extract:
	input:
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/md_completed.txt",
		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.xtc",
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.tpr",
	output:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/trjconv_pbc.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/trjconv_pbc.log"
		),
		tar_abs = lambda wildcards: os.path.abspath(
			f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/frames/FRAMES_compressed.tar.gz"
		)
	shell: 
		"""
		LOG_ABS="{params.log_abs}"
		TAR_ABS="{params.tar_abs}"
		
		mkdir -p $(dirname "$LOG_ABS")
		mkdir -p $(dirname "$TAR_ABS")

		# Change execution directory to md_results
		cd $(dirname {input.xtc_md})

		# 1. Recenter molecular unity boundaries
		echo "Protein" | gmx trjconv \
			-f $(basename {input.xtc_md}) \
			-s $(basename {input.tpr_md}) \
			-o md_whole.xtc \
			-pbc mol -ur compact > "$LOG_ABS" 2>&1

		# 2. Fit rotational/translational drift (printf with double backslash guarantees two distinct lines)
		printf "Protein Protein Protein" | gmx trjconv \
			-f md_whole.xtc \
			-s $(basename {input.tpr_md}) \
			-o md_clean.xtc \
			-center -fit rot+trans >> "$LOG_ABS" 2>&1

		# 3. Chop trajectory into individual PDB frames
		mkdir -p FRAMES
		echo "Protein" | gmx trjconv \
			-f md_clean.xtc \
			-s $(basename {input.tpr_md}) \
			-o FRAMES/frame.pdb \
			-sep >> "$LOG_ABS" 2>&1

		# 4. Compress and archive coordinate frames into absolute output destination
		tar -czf "$TAR_ABS" FRAMES >> "$LOG_ABS" 2>&1

		# 5. Clean up temporary intermediate trajectory files and uncompressed frames
		rm -f md_whole.xtc md_clean.xtc
		rm -rf FRAMES
		"""

# --- STEP 7: AUTOMATED PYMOL MOVIE GENERATION ---
rule generate_pymol_movie:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	output:
		movie = "results/movies/{pdb}/{source}/{model_id}/{protocol}_md_trajectory.mov"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/pymol_render.log"
	shell:
		"python scripts/render_pymol_movie.py {input.tar} {output.movie} > {log} 2>&1"


# should automatically execute 
#rule pbc_correction_and_extract:
#	input:
#		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/md_completed.txt",
#		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.xtc",
#		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.tpr",
#	output:
#		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
#	log:
#		"logs/{pdb}/{source}/{model_id}/{protocol}/trjconv_pbc.log"
#	params:
#		log_abs = lambda wildcards: os.path.abspath(
#			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/trjconv_pbc.log"
#		)
#	shell: 
#		"""
#		LOG_ABS="{params.log_abs}"
#		cd $(dirname {input.xtc_md})
#
#		# 1. Recenter molecular unity boundaries
#		echo "Protein" | gmx trjconv -f $(basename {input.xtc_md}) -s $(basename {input.tpr_md}) -o md_whole.xtc -pbc mol -ur compact > "$LOG_ABS" 2>&1
#		
#		# 2. Fit rotational and translational structural drift
#		echo "Protein Protein Protein" | gmx trjconv -f md_whole.xtc -s $(basename {input.tpr_md}) -o md_clean.xtc -center -fit rot+trans >> "$LOG_ABS" 2>&1
#		
#		# 3. Chop trajectory into individual PDB frames
#		mkdir -p frames
#		echo "Protein" | gmx trjconv -f md_clean.xtc -o FRAMES/frame.pdb -s $(basename {input.tpr_md}) -sep >> "$LOG_ABS" 2>&1
#		
#		# Compress and archive individual coordinate files
#		mkdir -p $(dirname {output.tar})
#		#mkdir FRAMES
#
#		tar -czf {output.tar} -C . frames >> "$LOG_ABS" 2>&1
#		
#		# Clean up massive intermediate trajectories to preserve space
#		# rm -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc
		"""
#		"""
#		cd $(dirname {output.xtc_md})
#		# 1. Recenter molecular unity boundaries
#		echo "Protein" | gmx trjconv -f {input.xtc_md} -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/md_whole.xtc -pbc mol -ur compact > {log} 2>&1
#		
#		# 2. Fit rotational and translational structural drift
#		echo "Protein Protein Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_whole.xtc -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/md_clean.xtc -center -fit rot+trans >> {log} 2>&1
#		
#		# 3. Chop trajectory into individual PDB frames
#		mkdir -p results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/FRAMES
#		echo "Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_clean.xtc -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/FRAMES/frame.pdb -s {input.tpr_md} -sep >> {log} 2>&1
#		
#		# Compress and archive individual coordinate files
#		mkdir -p $(dirname {output.tar})
#		tar -czf {output.tar} -C results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results FRAMES >> {log} 2>&1
#		
#		# Clean up massive intermediate trajectories to preserve space
#		# rm -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc
#		"""
###############################################################
### NOTE: CUSTOM SIMULATIONS BUILT ON TOP OF STANDARD 100ns SIM
### NOTE: CUSTOM SIMULATIONS BUILT ON TOP OF STANDARD 100ns SIM
### NOTE: CUSTOM SIMULATIONS BUILT ON TOP OF STANDARD 100ns SIM
### TODO: separate .py for reused rules
###############################################################

ruleorder: schedule_md_job > schedule_custom_md
ruleorder: create_md_job > create_custom_md
ruleorder: run_HPC_md > run_custom_md_HPC
ruleorder: run_local_md > run_custom_md_local
ruleorder: finalize_md > finalize_custom_md
# resource allocation
rule schedule_custom_md:
	input:
		#gro_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		#top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		#mdp = "config/gromacs_settings/interruptable_config_ultimate/md.mdp",
		std_done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/md_results/md_completed.txt",
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/md_results/{pdb}_{source}_{model_id}_md.tpr",
	output:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/scheduling.yml",
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md.tpr",
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/schedule_md_custom_job.log"
	params:
		#mdp_abs = lambda wildcards, input: os.path.abspath(input.mdp),
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/schedule_md_custom_job.log"
		)

	shell:
		'python scripts/sims/schedule_md_job.py - - - -'




# job preparation
rule create_custom_md:
	input:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/scheduling.yml"
	output:
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/create_custom_md_job.log"
	params:
		#log_abs = lambda wildcards, log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/create_custom_md_job.log"
		)
	shell:
		"python scripts/sims/create_md_job.py - - - -"

# submit job
rule run_custom_md_local:
	input:
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"

		#tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.tpr",
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/local_md_completed.txt"
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/local_production_mdrun.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/local_custom_production_mdrun.log"
		)
		
	shell:
		'python scripts/sims/run_md_job.py - - - -'

rule run_custom_md_HPC:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/HPC_md_completed.txt"
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/HPC_production_mdrun.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/HPC_custom_production_mdrun.log"
		)
	resources:
		gpu = 0,     # Zero local GPUs used! Passive SSH / monitoring thread only
		mem_mb = 500
#	shell:
#		'python scripts/sims/run_md_job.py - - - -'
	run:
		execute_hpc_md(
			wildcards=wildcards,
			job_description_path=input.job_description,
			output_done_path=output.done,
			log_path=params.log_abs,
			hpc_user=HPC_USER,
			hpc_host=HPC_HOST,
			hpc_remote_base=HPC_REMOTE_BASE
		)

rule finalize_custom_md:
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
		#protocol = "[a-zA-Z0-9_]+(?<!standard_100ns)" # Matches anything except 'standard_100ns'
	input:
		def_compute_scheduling = det_compute_scheduling
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/md_completed.txt",
		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.xtc",
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_md.tpr",

	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/finalize_custom_md.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/finalize_custom_md.log"
		)		
	shell:
		'python scripts/sims/finalize_md_job.py - - - -> {log} 2>&1'



# INACTIVATE!!
#rule run_molecular_dynamics:
#	input:
#		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
#		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
#	output:
#		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt"
#	log:
#		"logs/{pdb}/{source}/{model_id}/production_mdrun.log"
#	params:
#		log_abs = lambda wildcards: os.path.abspath(
#			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/production_mdrun.log"
#		)
#	resources:
#		gpu = 1
#	run:
#		# Setup Logger
#		log_path = params.log_abs
#		os.makedirs(os.path.dirname(log_path), exist_ok=True)
#
#		logger = logging.getLogger(f"run_md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}")
#		logger.setLevel(logging.INFO)
#		logger.handlers.clear()
#
#		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
#		fh = logging.FileHandler(log_path, mode="a")
#		fh.setFormatter(fmt)
#		logger.addHandler(fh)
#
#		ch = logging.StreamHandler()
#		ch.setFormatter(fmt)
#		logger.addHandler(ch)
#
#		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}"
#		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
#		job_dir = os.path.dirname(input.job_description)
#		local_dir = os.path.abspath(os.path.dirname(output.done))
#		step_log_file = os.path.join(local_dir, "simulation_steps.log")
#
#
#
#		logger.info("Starting MD execution phase for target: %s", target_id)
#		
#		# Preserve first start time if log exists, or mark INIT
#		if not os.path.exists(step_log_file):
#			log_sim_step("INIT", f"Target: {target_id}")
#		else:
#			log_sim_step("RESUME_WORKFLOW", f"Workflow checked target: {target_id}")
#
#		try:
#			with open(input.job_description, "r") as f:
#				job_description = f.read().strip()
#
#			# -------------------------------------------------------------
#			# LOCAL EXECUTION PATHWAY
#			# -------------------------------------------------------------
#			if job_description == "BLANK":
#				logger.info("Target configured for LOCAL execution.")
#				log_sim_step("EXEC_MODE", "Local execution requested")
#				cpt_path = os.path.join(job_dir, f"{prefix}.cpt")
#
#				if os.path.exists(cpt_path):
#					logger.info("Active checkpoint detected. Resuming local MD run...")
#					log_sim_step("MD_RESUME", f"Resuming from checkpoint {prefix}.cpt")
#					shell("""
#						cd "{job_dir}"
#						gmx mdrun -v -ntmpi 1 \
#							-deffnm {prefix} \
#							-cpi $(basename {cpt_path}) \
#							-nb gpu -pme gpu >> "{log_path}" 2>&1
#					""")
#				else:
#					logger.info("Launching fresh local MD run...")
#					log_sim_step("MD_START", f"Starting fresh mdrun for {prefix}")
#					shell("""
#						cd "{job_dir}"
#						gmx mdrun -v -ntmpi 1 \
#							-deffnm {prefix} \
#							-nb gpu -pme gpu >> "{log_path}" 2>&1
#					""")
#
#				shell("""
#					cd "{job_dir}"
#					[ -f "{prefix}.xtc" ] && ln -sf "{prefix}.xtc" ../md.xtc
#					[ -f "{prefix}.tpr" ] && ln -sf "{prefix}.tpr" ../md.tpr
#				""")
#				log_sim_step("MD_COMPLETE", "Local run finished successfully")
#
#			# HPC EXECUTION PATHWAY
#			# -------------------------------------------------------------
#			else:
#				ssh_target = f"{HPC_USER}@{HPC_HOST}"
#				clean_base = HPC_REMOTE_BASE.lstrip("~/")
#				remote_dir = os.path.join(
#					clean_base,
#					f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns",
#				)
#				job_name = (
#					f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
#				)
#				job_script = os.path.basename(input.job_description)
#
#				log_sim_step("EXEC_MODE", f"HPC execution target ({ssh_target})")
#
#				# 1. Check remote completion upfront
#				check_cmd = [
#					"ssh",
#					ssh_target,
#					f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES || echo NO",
#				]
#				remote_done = (
#					subprocess.check_output(check_cmd, text=True).strip()
#				)
#
#				if remote_done == "YES":
#					logger.info(
#						"🎉 Remote simulation already finished! (Found %s.gro"
#						" on HPC)",
#						prefix,
#					)
#					log_sim_step(
#						"CHECK_REMOTE",
#						f"Simulation finished ({prefix}.gro present)",
#					)
#				else:
#					# 2. Check if job is active in queue (Returns Job ID and State)
#					squeue_cmd = [
#						"ssh",
#						ssh_target,
#						f"squeue -u {HPC_USER} -n {job_name} -h -o '%i %t'",
#					]
#					q_res = subprocess.run(
#						squeue_cmd, capture_output=True, text=True, check=False
#					)
#					q_out = q_res.stdout.strip()
#
#					job_id = None
#					if q_out:
#						# Extract Job ID if already running/pending
#						job_id = q_out.split()[0]
#						logger.info(
#							"⏳ Slurm Job '%s' already active in queue (Job"
#							" ID: %s).",
#							job_name,
#							job_id,
#						)
#						log_sim_step(
#							"SLURM_ACTIVE",
#							f"Job ID {job_id} running/queued on HPC",
#						)
#					else:
#						logger.info(
#							"🚀 Submitting new Slurm job to Komondor (%s)...",
#							ssh_target,
#						)
#						submit_cmd = [
#							"ssh",
#							ssh_target,
#							f"cd '{remote_dir}' && [ -f JOB.tar.gz ] && tar"
#							f" -xzf JOB.tar.gz && sbatch {job_script}",
#						]
#						sub_res = subprocess.run(
#							submit_cmd,
#							capture_output=True,
#							text=True,
#							check=False,
#						)
#
#						if (
#							sub_res.returncode == 0
#							and "Submitted batch job" in sub_res.stdout
#						):
#							job_id = sub_res.stdout.strip().split()[-1]
#							logger.info(
#								"✅ Slurm job successfully submitted!"
#								" Assigned Job ID: %s",
#								job_id,
#							)
#							log_sim_step(
#								"SLURM_SUBMIT",
#								f"Submitted batch job ID: {job_id}",
#							)
#						else:
#							logger.error(
#								"❌ Failed to submit Slurm job: %s",
#								sub_res.stderr,
#							)
#							log_sim_step(
#								"SLURM_FAILED",
#								f"Submission failed: {sub_res.stderr.strip()}",
#							)
#							raise RuntimeError(
#								f"Slurm sbatch submission failed: {sub_res.stderr}"
#							)
#
#					# 3. Monitor Slurm Execution Loop
#					if job_id:
#						logger.info("Monitoring Slurm Job %s...", job_id)
#						last_progress_msg = ""
#
#						while True:
#							# Query job state safely (check=False avoids CalledProcessError when job finishes)
#							if job_id and str(job_id).isdigit():
#								check_q = [
#									"ssh",
#									ssh_target,
#									f"squeue -j {job_id} -h -o '%t'",
#								]
#								q_check = subprocess.run(
#									check_q,
#									capture_output=True,
#									text=True,
#									check=False,
#								)
#							
#								# Clean stdout: take only the first token (e.g. "R" or "PD")
#								stdout_clean = q_check.stdout.strip()
#								job_state = stdout_clean.split()[0] if stdout_clean else ""
#							
#								if not job_state:
#									logger.info(
#										"Job %s left the queue. Verifying completion...",
#										job_id,
#									)
#									break  # Job completed or died; exit polling loop
#							
#								if job_state == "R":
#									progress = get_remote_gromacs_progress(ssh_target, remote_dir, prefix)
#									if progress and progress != last_progress_msg:
#										log_sim_step(
#											"HEARTBEAT",
#											f"Job {job_id} running - {progress}",
#										)
#										last_progress_msg = progress
#									else:
#										log_sim_step(
#											"HEARTBEAT",
#											f"Job {job_id} actively executing on HPC",
#										)
#								else:
#									log_sim_step(
#										"HEARTBEAT",
#										f"Job {job_id} queued (State: {job_state})",
#									)
#
#							time.sleep(180)  # Poll every 3 minutes
#
#						# 4. Final Verification: Confirm simulation produced output
#						post_check = [
#							"ssh",
#							ssh_target,
#							f"[ -f '{remote_dir}/{prefix}.gro' ] && echo YES"
#							" || echo NO",
#						]
#						final_done = (
#							subprocess.check_output(post_check, text=True)
#							.strip()
#						)
#
#						if final_done == "YES":
#							log_sim_step(
#								"SLURM_FINISHED",
#								f"Slurm Job {job_id} completed successfully",
#							)
#						else:
#							log_sim_step(
#								"SLURM_FAILED",
#								f"Slurm Job {job_id} ended without producing"
#								f" {prefix}.gro",
#							)
#							raise RuntimeError(
#								f"HPC Job {job_id} terminated unexpectedly"
#								f" ({prefix}.gro missing)."
#							)
#
#				# Sync and retrieve results locally
#				logger.info("📦 Archiving and retrieving simulation artifacts from HPC...")
#				log_sim_step("RETRIEVE_START", "Fetching remote output files via rsync")
#
#
#				shell(f"""
#					SSH_TARGET="{ssh_target}"
#					REMOTE_DIR="{remote_dir}"
#					PREFIX="{prefix}"
#					LOCAL_DIR="{local_dir}"
#					LOG_FILE="{log_path}"
#
#					ssh "$SSH_TARGET" "cd '$REMOTE_DIR' && tar -czf md_results.tar.gz $PREFIX.*" >> "$LOG_FILE" 2>&1
#					rsync -avz "$SSH_TARGET:$REMOTE_DIR/md_results.tar.gz" "$LOCAL_DIR/" >> "$LOG_FILE" 2>&1
#					mkdir -p "$LOCAL_DIR/md_results"
#					tar -xzf "$LOCAL_DIR/md_results.tar.gz" -C "$LOCAL_DIR/md_results"
#
#					# I don't want to bother with this yet.
#					#[ -f "$LOCAL_DIR/$PREFIX.xtc" ] && ln -sf "$PREFIX.xtc" "$LOCAL_DIR/md.xtc"
#					#[ -f "$LOCAL_DIR/$PREFIX.tpr" ] && ln -sf "$PREFIX.tpr" "$LOCAL_DIR/md.tpr"
#
#					rm -f "$LOCAL_DIR/md_results.tar.gz"
#					ssh "$SSH_TARGET" "rm -f '$REMOTE_DIR/md_results.tar.gz'" >> "$LOG_FILE" 2>&1
#				""")
#
#				logger.info(f"MD DATA ACQUIRED CHECK @ {local_dir}/md_results")
#
#				log_sim_step("RETRIEVE_COMPLETE", "Downloaded, extracted, and cleaned up md_results.tar.gz")
#
#			# Sentinel output
#			with open(output.done, "w") as f:
#				f.write(f"MD simulation completed for {target_id}.\n")
#			
#			log_sim_step("DONE", "Pipeline rule finished successfully.")
#
#		except Exception as err:
#			log_sim_step("ERROR", str(err))
#			logger.exception("Execution failed in run_molecular_dynamics for %s: %s", target_id, str(err))
#			raise

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
#rule pbc_correction_and_extract:
#	input:
#		md_checkpoint_guard = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt",
#		xtc_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md.xtc",
#		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md.tpr",
#	output:
#		tar = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/FRAMES_compressed.tar.gz"
#	log:
#		"logs/{pdb}/{source}/{model_id}/trjconv_pbc.log"
#	shell:
#		"""
#		# 1. Recenter molecular unity boundaries
#		echo "Protein" | gmx trjconv -f {input.xtc_md} -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc -pbc mol -ur compact > {log} 2>&1
#		
#		# 2. Fit rotational and translational structural drift
#		echo "Protein Protein Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc -s {input.tpr_md} -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc -center -fit rot+trans >> {log} 2>&1
#		
#		# 3. Chop trajectory into individual PDB frames
#		mkdir -p results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/FRAMES
#		echo "Protein" | gmx trjconv -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc -o results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/FRAMES/frame.pdb -s {input.tpr_md} -sep >> {log} 2>&1
#		
#		# Compress and archive individual coordinate files
#		tar -czf {output.tar} -C results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id} FRAMES >> {log} 2>&1
#		
#		# Clean up massive intermediate trajectories to preserve space
#		# rm -f results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_whole.xtc results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns/md_clean.xtc
#		"""
#
#
		
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
