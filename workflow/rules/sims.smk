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

#from scripts.smk.sims.hpc_engine import HPCJobManager
sys.path.insert(0, os.path.abspath("scripts/smk/sims"))
from hpc_engine import HPCJobManager

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

	##job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{protocol}"

	if protocol is None:
		protocol = getattr(wildcards, "protocol", "standard_100ns")

	protocol_stem = protocol.split("_")[0]
	job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{protocol}"

	time_limit = HPC_TIME  # Default time limit from environment variable

	if protocol == "standard_100ns":
		deffnm = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
	else:
		#time_limit = "48:00:00"  # Set a longer time limit for extended or other protocols
		time_limit = '7-00:00:00' # 

		deffnm = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{protocol}_md"

	# Calculate absolute target total time in picoseconds (-until)
	length_match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns|ps)', protocol, re.IGNORECASE)
	if length_match:
		val, unit = float(length_match.group(1)), length_match.group(2).lower()
		if unit == "us":
			until_ps = int(val * 1_000_000)
		elif unit == "ns":
			until_ps = int(val * 1_000)
		else:
			until_ps = int(val)
	else:
		until_ps = 1_000_000

	# Include 'standard_100ns' and 'standard' in valid protocols
	valid_protocols = ["standard", "standard_100ns", "simulated-annealing", "extended"]

	if protocol_stem not in valid_protocols:
		raise ValueError(f"{protocol} not recognized as an approved protocol for slurm job")
	
	if protocol_stem=='extended':
		md_execution_block = f"""\
# --- Extended MD Execution Path ---
# --- Extended MD Execution Path (Option B: -until {until_ps} ps) ---
if [ -f "parent.cpt" ] && [ -f "parent.tpr" ] && [ ! -f "{deffnm}.cpt" ]; then
	echo "--> Extending TPR until target total duration ({until_ps} ps)..."
	gmx_mpi convert-tpr -s "parent.tpr" -until {until_ps} -o "{deffnm}.tpr"	
	echo "--> Starting extended run from previous checkpoint..."
	srun --nodes=2 \\
		--ntasks-per-node=16 \\
		--cpus-per-task=8 \\
		--cpu-bind=cores \\
		gmx_mpi mdrun -v -deffnm "{deffnm}" \\
			-ntomp 8 \\
			-dds 0.8 \\
			-rcon 0 \\
			-dlb yes \\
			-cpi "{deffnm}_prev.cpt"
elif [ -f "{deffnm}.cpt" ] && [ ! -f "{deffnm}.gro" ]; then
	echo "--> Resuming active extended run from local checkpoint..."
	srun --nodes=2 \\
		--ntasks-per-node=16 \\
		--cpus-per-task=8 \\
		--cpu-bind=cores \\
		gmx_mpi mdrun -v -deffnm "{deffnm}" \\
			-ntomp 8 \\
			-dds 0.8 \\
			-rcon 0 \\
			-dlb yes \\
			-cpi "{deffnm}.cpt"
else
	if [ -f "{deffnm}.gro" ]; then
		echo "--> Extended MD trajectory already complete ({deffnm}.gro found)."
	else
		echo "❌ ERROR: Required previous checkpoint ({deffnm}_prev.cpt) not found!"
		exit 1
	fi
fi
"""

# --- Helper Function for MDRun ---
# This checks if a checkpoint exists for the specific step to resume it
	else:
# 2. Logic for Standard & Simulated Annealing Runs	else:
		md_execution_block = f"""\
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

	#else:
	#	raise ValueError(f"{protocol} not recognized as an approved protocol for slurm job")

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
#SBATCH --time={time_limit}
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

{md_execution_block}
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

		current_protocol = "standard_100ns"
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

	# Access config entries using combined nested key
	nested_key = f"{wildcards.source}_{wildcards.model_id}"
	pdb_config_entries = (
		config.get("custom_simulations", {})
		.get(wildcards.pdb, {})
		.get(nested_key, [])
	)

	pdb_config = next((entry for entry in pdb_config_entries if entry.get("protocol") == protocol), {})
	compute_target = pdb_config.get("compute_target", def_compute_target)

	if compute_target in ["HPC", "local"]:
		return compute_target

	# Fallback to checking scheduling.yml if present
	scheduling_file = Path(
		f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{protocol}/JOB/scheduling.yml"
	)

	if scheduling_file.exists():
		try:
			with open(scheduling_file, 'r') as f:
				data = yaml.safe_load(f) or {}
				target = data.get("COMPUTE")
				if target in ["HPC", "local"]:
					return target
		except Exception:
			pass

	return def_compute_target


def det_compute_scheduling(wildcards):
	protocol = getattr(wildcards, "protocol", "standard_100ns")
	mode = get_compute_target(wildcards)

	if mode in ['local', 'HPC']:
		# MUST MATCH output in run_custom_md_local / run_custom_md_HPC
		return f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{protocol}/{mode}_md_completed.txt"
	else:
		raise ValueError(f"Compute type '{mode}' not permitted.")

def log_sim_step(step_name, step_log_file, details=""):
	"""Writes timestamped step tracking events to simulation_steps.log."""
	timestamp = subprocess.check_output('date +"%Y-%m-%d %H:%M:%S"', shell=True).decode().strip()
	line = f"[{timestamp}] [STEP: {step_name}] {details}\n"
	with open(step_log_file, "a") as f:
		f.write(line)

def get_simulation_progress(ssh_target: str, remote_dir: str) -> str | None:
	"""Helper to parse the latest step/time progress from remote md.log or slurm .out files."""
	cmd = (
		f"ssh -o BatchMode=yes {ssh_target} "
		f"\"for f in '{remote_dir}'/*.log '{remote_dir}'/*.out '{remote_dir}'/JOB/*.log '{remote_dir}'/JOB/*.out; do "
		f"[ -f \\\"\$f\\\" ] && tail -n 30 \\\"\$f\\\" | grep -E 'Step|Time|Vol|ETA' | tail -n 1; "
		f"done | tail -n 1\""
	)
	try:
		res = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
		out = res.stdout.strip()
		return out if out else None
	except Exception:
		return None
		
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

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_standard_100ns_md"
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
		fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
		ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns"
		
		# FIX: Standard 100ns uses standard prefix without protocol suffix
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
		
		job_dir = os.path.dirname(input.job_description)
		local_dir = os.path.abspath(os.path.dirname(output.done))
		step_log_file = os.path.join(job_dir, "simulation_steps.log")

		clean_base = HPC_REMOTE_BASE.lstrip("~/")
		remote_dir = os.path.join(clean_base, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns")
		#job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
		job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_standard_100ns"

		hpc = HPCJobManager(
			ssh_target="komondor",
			remote_dir=remote_dir,
			local_dir=local_dir,
			job_name=job_name,
			job_script=os.path.basename(input.job_description),
			log_path=log_path,
			logger=logger,
			step_log_file=step_log_file,
			hpc_user=HPC_USER,
		)

		try:
			def md_progress(ssh_target, r_dir):
				return get_remote_gromacs_progress(ssh_target, r_dir, prefix)

			hpc.execute_pipeline(
				completion_check_file=f"{prefix}.gro",
				target_subdir="./md_results_HPC/",
				poll_interval_sec=60,
				get_progress_fn=md_progress,
				unpack_job_archive=True
			)

			with open(output.done, "w") as f:
				f.write(f"HPC MD simulation completed for {target_id}.\n")

			log_sim_step("DONE", step_log_file, "HPC MD rule finished successfully.")

		except Exception as err:
			log_sim_step("ERROR", step_log_file, str(err))
			logger.exception("Execution failed for %s: %s", target_id, str(err))
			raise

rule finalize_md:
	input:
		det_compute_scheduling = det_compute_scheduling
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/md_completed.txt",
		xtc_md = protected("results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.xtc"),
		tpr_md = protected("results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.tpr"),
		cpt_md = protected("results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_results/{pdb}_{source}_{model_id}_md.cpt")  # <-- CRITICAL FOR DAG RESOLUTION
	log:
		"logs/{pdb}/{source}/{model_id}/finalize_md.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/finalize_md.log"
		)
	run:
		# ------------------------------------------------------------------
		# SHORT-CIRCUIT: Skip if already finalized locally
		# ------------------------------------------------------------------
		out_dir = os.path.abspath(os.path.dirname(output.done))

		if os.path.exists(output.done) and os.path.exists(output.xtc_md) and os.path.exists(output.tpr_md):
			print(f"Target {wildcards.pdb}/{wildcards.source}/{wildcards.model_id} already finalized locally. Skipping HPC pull.")
			return

		# Setup Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger = logging.getLogger(f"finalize_md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}")
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
		fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
		ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
		compute_target = get_compute_target(wildcards)

		os.makedirs(out_dir, exist_ok=True)

		protocol_dir = os.path.abspath(f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns")
		job_dir = os.path.join(protocol_dir, "JOB")
		hpc_staging_dir = os.path.join(protocol_dir, "md_results_HPC")

		logger.info("Starting finalize_md phase for %s [Compute Target: %s]", target_id, compute_target)

		try:
			# ------------------------------------------------------------------
			# HPC PATH: Relocate from md_results_HPC -> md_results
			# ------------------------------------------------------------------
			if compute_target == "HPC":
				if os.path.exists(hpc_staging_dir):
					logger.info("Staging directory found at %s. Moving artifacts to %s...", hpc_staging_dir, out_dir)
					moved_count = 0

					for item in os.listdir(hpc_staging_dir):
						src = os.path.join(hpc_staging_dir, item)
						dst = os.path.join(out_dir, item)

						if os.path.isfile(src):
							shutil.move(src, dst)
							logger.info("Moved file: %s -> %s", item, dst)
							moved_count += 1
						elif os.path.isdir(src):
							logger.info("Scanning nested directory inside staging: %s", item)
							for nested_fn in os.listdir(src):
								nested_src = os.path.join(src, nested_fn)
								nested_dst = os.path.join(out_dir, nested_fn)
								shutil.move(nested_src, nested_dst)
								logger.info("Moved nested file: %s -> %s", nested_fn, nested_dst)
								moved_count += 1

					logger.info("Successfully relocated %d items from HPC staging.", moved_count)
					shutil.rmtree(hpc_staging_dir, ignore_errors=True)
					logger.info("Cleaned up staging folder: %s", hpc_staging_dir)
				else:
					logger.warning("HPC staging directory %s does not exist. Checking protocol root for fallback files...", hpc_staging_dir)
					# Fallback: Move matching files directly from protocol_dir if they were downloaded at protocol root
					for ext in [".xtc", ".tpr", ".cpt", ".gro", ".edr", ".log"]:
						src = os.path.join(protocol_dir, f"{prefix}{ext}")
						dst = os.path.join(out_dir, f"{prefix}{ext}")
						if os.path.exists(src) and not os.path.exists(dst):
							shutil.move(src, dst)
							logger.info("Moved fallback file from protocol root: %s -> %s", src, dst)

				# Fallback check for generic 'md.xtc' naming
				generic_xtc = os.path.join(out_dir, "md.xtc")
				if not os.path.exists(output.xtc_md) and os.path.exists(generic_xtc):
					os.rename(generic_xtc, output.xtc_md)
					logger.info("Renamed generic 'md.xtc' to expected '%s'", output.xtc_md)

				# Fallback copy for .tpr from local JOB folder if missing
				local_tpr_src = os.path.join(job_dir, f"{prefix}.tpr")
				if not os.path.exists(output.tpr_md) and os.path.exists(local_tpr_src):
					shutil.copy2(local_tpr_src, output.tpr_md)
					logger.info("Copied fallback TPR structure from %s to %s", local_tpr_src, output.tpr_md)

				# Validation checks
				if not os.path.exists(output.xtc_md):
					err_msg = f"HPC execution finished, but {output.xtc_md} was not found in {out_dir}!"
					logger.error(err_msg)
					raise FileNotFoundError(err_msg)

				if not os.path.exists(output.tpr_md):
					err_msg = f"HPC execution finished, but {output.tpr_md} is missing from {out_dir}!"
					logger.error(err_msg)
					raise FileNotFoundError(err_msg)

				with open(output.done, "w") as f:
					f.write(f"MD simulation (HPC) finalized for {target_id}\n")

				logger.info("✅ HPC MD finalization successfully completed for %s", target_id)
				return

			# ------------------------------------------------------------------
			# LOCAL PATH
			# ------------------------------------------------------------------
			if compute_target == "local":
				job_desc = os.path.join(job_dir, f"{prefix}_job.job")

				if not os.path.exists(job_desc):
					err_msg = f"Expected job description not found: {job_desc}"
					logger.error(err_msg)
					raise FileNotFoundError(err_msg)

				with open(job_desc, "r") as fh:
					content = fh.read().strip()

				if content != "BLANK":
					err_msg = f"Job description must be 'BLANK' for local finalize, got: '{content}'"
					logger.error(err_msg)
					raise ValueError(err_msg)

				if os.path.isdir(job_dir):
					logger.info("Collecting local MD outputs from %s into %s...", job_dir, out_dir)
					copied_count = 0
					for fn in os.listdir(job_dir):
						if fn.startswith(prefix):
							src = os.path.join(job_dir, fn)
							dst = os.path.join(out_dir, fn)
							try:
								shutil.copy2(src, dst)
								copied_count += 1
							except Exception as c_err:
								logger.warning("Failed to copy %s to %s: %s", src, dst, str(c_err))
					logger.info("Copied %d local simulation artifacts.", copied_count)

				generic_xtc = os.path.join(out_dir, "md.xtc")
				if not os.path.exists(output.xtc_md) and os.path.exists(generic_xtc):
					os.rename(generic_xtc, output.xtc_md)
					logger.info("Renamed generic 'md.xtc' to '%s'", output.xtc_md)

				if not os.path.exists(output.xtc_md) or not os.path.exists(output.tpr_md):
					err_msg = f"Local MD completed, but missing expected outputs (.xtc / .tpr) in {out_dir}"
					logger.error(err_msg)
					raise FileNotFoundError(err_msg)

				with open(output.done, "w") as f:
					f.write(f"MD simulation (local) finalized and artifacts collected for {target_id}\n")

				logger.info("✅ Local MD finalization successfully completed for %s", target_id)
				return

			raise ValueError(f"Unknown compute target when finalizing MD: {compute_target}")

		except Exception as err:
			logger.exception("Finalization failed for target %s: %s", target_id, str(err))
			raise
				
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

def get_pbc_correction_input_deps(wildcards):
	protocol = wildcards.protocol
	base_path = f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{protocol}/md_results"

	if protocol == "standard_100ns":
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
	else:
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{protocol}_md"

	return {
		"done": f"{base_path}/md_completed.txt",
		"xtc_md": f"{base_path}/{prefix}.xtc",
		"tpr_md": f"{base_path}/{prefix}.tpr"
	}

rule pbc_correction_and_extract:
	input:
		unpack(get_pbc_correction_input_deps)
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

		echo "Protein Protein Protein" | gmx trjconv \
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


def parse_length_ns(full_protocol, cfg_dict):
	"""
	Parses simulation length in nanoseconds.
	1. Checks explicit `length_ns` key in config dictionary.
	2. Parses the duration suffix from full_protocol (e.g., 'extended_1000ns' -> 1000.0, 'extended_1us' -> 1000.0).
	"""
	if "length_ns" in cfg_dict:
		return float(cfg_dict["length_ns"])

	# Match numbers followed by 'us' or 'ns' anywhere in the protocol name
	match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns)', full_protocol, re.IGNORECASE)
	if match:
		val, unit = float(match.group(1)), match.group(2).lower()
		return val * 1000.0 if unit == "us" else val

	raise ValueError(f"Could not parse length from protocol string '{full_protocol}' or config: {cfg_dict}")

def get_md_input_dependencies(wildcards):
	pdb = wildcards.pdb
	source = wildcards.source
	model_id = wildcards.model_id
	protocol = wildcards.protocol

	nested_key = f"{source}_{model_id}"
	entries = (
		config.get("custom_simulations", {})
		.get(pdb, {})
		.get(nested_key, [])
	)

	target_entry = next((e for e in entries if e.get("protocol") == protocol), None)

	# 1. Extension run: require parent checkpoint and TPR from parent's md_results/
	if target_entry and "md_to_extend" in target_entry:
		parent_protocol = target_entry["md_to_extend"]
		parent_dir = f"results/gromacs/{pdb}/{source}/{model_id}/{parent_protocol}/md_results"

		# Baseline uses standard naming; custom protocols include {parent_protocol} in stem
		if parent_protocol == "standard_100ns":
			parent_file_stem = f"{pdb}_{source}_{model_id}_md"
		else:
			parent_file_stem = f"{pdb}_{source}_{model_id}_{parent_protocol}_md"

		return {
			"cpt": f"{parent_dir}/{parent_file_stem}.cpt",
			"tpr": f"{parent_dir}/{parent_file_stem}.tpr"
		}

	# 2. Fresh run / Fallback: require EM outputs
	return {
		"gro_cg": f"results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		"top": f"results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top"
	}



# for extended runs... ahh
# stage resolution map: N -> N-1
#STAGE_PRECEDENCE = {
#	'100ns': None,
#	'1us': '100ns',
#	'2us': '1us',
#	'5us': '2us',
#	'10us': '5us',
#}



ruleorder: schedule_md_job > schedule_custom_md
ruleorder: create_md_job > create_custom_md_job
ruleorder: run_HPC_md > run_custom_md_HPC
ruleorder: run_local_md > run_custom_md_local
ruleorder: finalize_md > finalize_custom_md
# resource allocation

rule create_custom_mdp:
	input:
		config_file = "config/config.yaml"
	output:
		custom_mdp = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/generated.mdp"
	params:
		protocol = lambda wildcards: wildcards.protocol
	run:
		pdb = wildcards.pdb
		source = wildcards.source
		model_id = wildcards.model_id
		protocol = wildcards.protocol

		# Combine source and model_id to match config layout (e.g., 'empirical_canonical_structure')
		nested_key = f"{source}_{model_id}"

		# Retrieve target entries list from config
		entries = (
			config.get("custom_simulations", {})
			.get(pdb, {})
			.get(nested_key, [])
		)

		# Locate specific entry matching this wildcard protocol
		target_entry = next((e for e in entries if e.get("protocol") == protocol), None)

		if not target_entry:
			raise ValueError(
				f"Protocol '{protocol}' not found in custom_simulations for {pdb}/{nested_key}"
			)

		# Write matching protocol configuration to a temporary YAML file
		with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
			yaml.dump(target_entry, tmp, sort_keys=False)
			tmp_yaml_path = tmp.name

		try:
			# Generate target .mdp file
			shell("python scripts/smk/sims/generate_custom_mdp.py '{protocol}' '{tmp_yaml_path}' '{output.custom_mdp}'")
		finally:
			if os.path.exists(tmp_yaml_path):
				os.remove(tmp_yaml_path)



rule schedule_custom_md:
	input:
		unpack(get_md_input_dependencies),
		custom_mdp = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/generated.mdp"
	output:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/scheduling.yml",
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_{protocol}_md.tpr"
	wildcard_constraints:
		#protocol = "(?!standard_100ns$)[^/]+"
		#protocol = "^(?!standard_100ns$).+"   # Strict Regex: Matches everything EXCEPT 'standard_100ns'
		protocol = "(?!standard_100ns$)[^/]+"

	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/schedule_md_custom_job.log"
	params:
		mdp_abs = lambda wildcards, input: os.path.abspath(input.custom_mdp),
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/schedule_md_custom_job.log"
		)
	run:
		def_compute_target = config.get("default_compute_target", "local")
		nested_key = f"{wildcards.source}_{wildcards.model_id}"

		pdb_config_entries = (
			config.get("custom_simulations", {})
			.get(wildcards.pdb, {})
			.get(nested_key, [])
		)
		pdb_config = next((entry for entry in pdb_config_entries if entry.get("protocol") == wildcards.protocol), {})
		compute_target = pdb_config.get("compute_target", def_compute_target)

		work_dir = os.path.dirname(output.scheduling)
		os.makedirs(work_dir, exist_ok=True)
		os.makedirs(os.path.dirname(params.log_abs), exist_ok=True)

		exec_dir = work_dir
		tpr_path = os.path.abspath(output.tpr_file)

		# Calculate absolute target total time in picoseconds (-until)
		length_match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns|ps)', wildcards.protocol, re.IGNORECASE)
		if length_match:
			val, unit = float(length_match.group(1)), length_match.group(2).lower()
			if unit == "us":
				until_ps = int(val * 1_000_000)
			elif unit == "ns":
				until_ps = int(val * 1_000)
			else:
				until_ps = int(val)
		else:
			until_ps = 1_000_000  # Fallback to 1,000,000 ps (1 us total time)

		# Build .tpr file if it doesn't exist
		if not os.path.exists(tpr_path):
			if hasattr(input, "cpt") and hasattr(input, "tpr"):
				# Always stage parent files into JOB directory under standardized names
				shutil.copy2(input.cpt, os.path.join(work_dir, "parent.cpt"))
				shutil.copy2(input.tpr, os.path.join(work_dir, "parent.tpr"))

				# Option B: Extend simulation UNTIL total time equals until_ps
				shell("""
					gmx convert-tpr -s {work_dir}/parent.tpr -until {until_ps} -o {tpr_path} > {params.log_abs} 2>&1
				""")
			else:
				# Fresh non-extension runs (e.g., simulated annealing)
				shutil.copy(input.gro_cg, os.path.join(exec_dir, os.path.basename(input.gro_cg)))
				shutil.copy(input.top, os.path.join(exec_dir, os.path.basename(input.top)))
				shutil.copy(input.custom_mdp, os.path.join(exec_dir, os.path.basename(input.custom_mdp)))

				shell("""
					cd {exec_dir}
					gmx grompp -f {params.mdp_abs} \
						-o $(basename {output.tpr_file}) \
						-c $(basename {input.gro_cg}) \
						-r $(basename {input.gro_cg}) \
						-p $(basename {input.top}) -maxwarn 1 > {params.log_abs} 2>&1
				""")

				os.remove(os.path.join(exec_dir, os.path.basename(input.gro_cg)))
				os.remove(os.path.join(exec_dir, os.path.basename(input.top)))
				os.remove(os.path.join(exec_dir, os.path.basename(input.custom_mdp)))

		# Save scheduling info including calculated UNTIL_PS parameter for HPC batch script resolution
		if compute_target in ["local", "HPC"]:
			scheduling_info = {
				"COMPUTE": compute_target,
				"UNTIL_PS": until_ps
			}
			with open(output.scheduling, "w") as f:
				yaml.safe_dump(scheduling_info, f, sort_keys=False)
		else:
			raise ValueError(f"Compute target '{compute_target}' not permitted.")



# job preparation
rule create_custom_md_job:
	input:
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/scheduling.yml"
	output:
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"
	wildcard_constraints:
		protocol = "(?!standard_100ns$)[^/]+"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/create_custom_md_job.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/create_custom_md_job.log"
		)
	run:
		# 1. Setup Isolated Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger_name = f"create_custom_md_job_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"
		logger = logging.getLogger(logger_name)
		logger.setLevel(logging.INFO)
		logger.handlers.clear()  # Prevent duplicate handlers on re-runs

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

		fh = logging.FileHandler(log_path, mode="w")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
		logger.info("Initializing custom MD job creation for target: %s", target_id)

		try:
			# 2. Parse Scheduling Data
			if not os.path.exists(input.scheduling):
				raise FileNotFoundError(f"Scheduling file not found at {input.scheduling}")

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

				# Generate Slurm Batch Script dynamically passing wildcards and wildcards.protocol
				slurm_script = resolve_slurm_job_script(wildcards, protocol=wildcards.protocol)

				os.makedirs(os.path.dirname(script_path), exist_ok=True)
				with open(script_path, "w") as fh:
					fh.write(slurm_script)
				logger.info("Generated Slurm batch script for protocol '%s': %s", wildcards.protocol, script_path)

				# Create Tarball Archive
				job_dir = os.path.dirname(output.job_description)
				job_targz = f"{job_dir}.tar.gz"
				logger.info("Archiving directory '%s' into '%s'...", job_dir, job_targz)

				with tarfile.open(job_targz, "w:gz") as tar:
					for fn in os.listdir(job_dir):
						p = os.path.join(job_dir, fn)
						tar.add(p, arcname=os.path.basename(fn))
				logger.info("Archive created successfully (Size: %.2f KB)", os.path.getsize(job_targz) / 1024.0)

				# Remote Sync via SSH/Rsync to remote protocol path
				remote_dir = os.path.join(HPC_REMOTE_BASE, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}")
				ssh_target = "komondor"
				logger.info("Deploying archive to remote target %s:%s", ssh_target, remote_dir)

				shell("""
					ssh -o BatchMode=yes {ssh_target} "mkdir -p '{remote_dir}'"
					rsync -e "ssh -o BatchMode=yes" -avz "{job_targz}" "{ssh_target}:{remote_dir}/JOB.tar.gz"
					ssh -o BatchMode=yes {ssh_target} "test -f '{remote_dir}/JOB.tar.gz' && echo 'Remote payload verified at {remote_dir}/JOB.tar.gz'"
				""")
				logger.info("Remote transfer and payload verification completed.")

				# Append Submission Metadata
				scheduling_info['JOB_STATUS'] = 'Submitted'
				scheduling_info['REMOTE_DIR'] = str(remote_dir)

				with open(input.scheduling, 'w') as f:
					yaml.safe_dump(scheduling_info, f, sort_keys=False)

				logger.info("Updated scheduling file '%s' with submission metadata.", input.scheduling)

			else:
				logger.error("Invalid COMPUTE target '%s' in %s", compute_target, input.scheduling)
				raise ValueError(f"Unknown compute target in scheduling info: {compute_target}")

			logger.info("Rule create_custom_md_job completed successfully for %s.", target_id)

		except Exception as err:
			logger.exception("Execution failed in create_custom_md_job for %s: %s", target_id, str(err))
			raise
#rule create_custom_md_job:
#	input:
#		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/scheduling.yml"
#	output:
#		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"
#	wildcard_constraints:
#		protocol = "(?!standard_100ns$)[^/]+"
#	log:
#		"logs/{pdb}/{source}/{model_id}/{protocol}/create_custom_md_job.log"
#	params:
#		#log_abs = lambda wildcards, log: os.path.abspath(str(log))
#		log_abs = lambda wildcards: os.path.abspath(
#			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/create_custom_md_job.log"
#		)
#	shell:
#		"python scripts/sims/create_md_job.py - - - -"

# submit job
rule run_custom_md_local:
	input:
		tpr_md = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_{protocol}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md_job.job"
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
	resources:
		gpu = 1
	run:
		# 1. Setup Isolated Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger_name = f"run_custom_md_local_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"
		logger = logging.getLogger(logger_name)
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

		fh = logging.FileHandler(log_path, mode="a")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}_md"
		
		job_dir = os.path.abspath(os.path.dirname(input.job_description))
		local_dir = os.path.abspath(os.path.dirname(output.done))
		step_log_file = os.path.join(job_dir, "simulation_steps.log")

		logger.info("Starting Custom MD execution phase for target: %s", target_id)

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
				base_protocol = wildcards.protocol.split("_")[0]

				# Check if this is an extension or an active checkpoint resume
				if os.path.exists(cpt_path):
					logger.info("Active checkpoint detected at %s. Resuming local MD run...", cpt_path)
					log_sim_step("MD_RESUME", step_log_file, f"Resuming from checkpoint {prefix}.cpt")
					shell(f"""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-cpi {prefix}.cpt \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")
				else:
					logger.info("Launching fresh local MD run for protocol %s...", wildcards.protocol)
					log_sim_step("MD_START", step_log_file, f"Starting fresh mdrun for {prefix}")
					shell(f"""
						cd "{job_dir}"
						gmx mdrun -v -ntmpi 1 \
							-deffnm {prefix} \
							-nb gpu -pme gpu >> "{log_path}" 2>&1
					""")

				# Copy outputs from JOB/ directory to the main protocol directory
				protocol_dir = os.path.abspath(os.path.dirname(input.job_description) + "/..")

				# Copy outputs to top-level protocol directory
				shell(f"""
					cd "{job_dir}"
					[ -f "{prefix}.xtc" ] && cp -f "{prefix}.xtc" "{local_dir}/{prefix}.xtc" || true
					[ -f "{prefix}.tpr" ] && cp -f "{prefix}.tpr" "{local_dir}/{prefix}.tpr" || true
					[ -f "{prefix}.gro" ] && cp -f "{prefix}.gro" "{local_dir}/{prefix}.gro" || true
					[ -f "{prefix}.cpt" ] && cp -f "{prefix}.cpt" "{local_dir}/{prefix}.cpt" || true
				""")
				# Write sentinel file
				with open(output.done, "w") as f:
					f.write(f"Local Custom MD completed for {target_id}.\n")

				log_sim_step("MD_COMPLETE", step_log_file, "Local custom run finished successfully")
			else:
				logger.info("Target %s is configured for HPC mode (job description: '%s'). Skipping local execution.", target_id, job_description)

		except Exception as err:
			log_sim_step("ERROR", step_log_file, str(err))
			logger.exception("Execution failed in run_custom_md_local for %s: %s", target_id, str(err))
			raise


rule run_custom_md_HPC:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_{protocol}_md.tpr",
		#tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/JOB/{pdb}_{source}_{model_id}_md.tpr",
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
		gpu = 0,
		mem_mb = 500
	run:
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger_name = f"run_custom_md_HPC_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"
		logger = logging.getLogger(logger_name)
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
		fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
		ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}_md"
		job_dir = os.path.dirname(input.job_description)
		local_dir = os.path.abspath(os.path.dirname(output.done))
		step_log_file = os.path.join(job_dir, "simulation_steps.log")

		clean_base = HPC_REMOTE_BASE.lstrip("~/")
		remote_dir = os.path.join(clean_base, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}")
		#job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"
		job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"

		hpc = HPCJobManager(
			ssh_target="komondor",
			remote_dir=remote_dir,
			local_dir=local_dir,
			job_name=job_name,
			job_script=os.path.basename(input.job_description),
			log_path=log_path,
			logger=logger,
			step_log_file=step_log_file,
			hpc_user=HPC_USER,
		)

		try:
			#def md_progress(ssh_target, r_dir):
			#	return get_remote_gromacs_progress(ssh_target, r_dir, prefix)

			hpc.execute_pipeline(
				completion_check_file=f"{prefix}.gro",
				target_subdir=".",
				poll_interval_sec=900,
				get_progress_fn=get_simulation_progress,
				unpack_job_archive=True
			)

			with open(output.done, "w") as f:
				f.write(f"HPC MD simulation completed for {target_id}.\n")

		except Exception as err:
			logger.exception("Execution failed in run_custom_md_HPC for %s: %s", target_id, str(err))
			raise

rule finalize_custom_md:
	wildcard_constraints:
		#protocol = "^(?!standard_100ns$).+"
		protocol = "(?!standard_100ns$)[^/]+"

	input:
		done_flag = det_compute_scheduling
	output:
		done = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/md_completed.txt",
		xtc_md = protected("results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_{protocol}_md.xtc"),
		tpr_md = protected("results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_{protocol}_md.tpr"),
		cpt_md = protected("results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/{pdb}_{source}_{model_id}_{protocol}_md.cpt")  # <-- REQUIRED FOR MULTI-STAGE EXTENSIONS
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/finalize_custom_md.log"
	params:
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/finalize_custom_md.log"
		)
	run:
		# SHORTCIRCUIT
		if os.path.exists(output.done) and os.path.exists(output.xtc_md) and os.path.exists(output.tpr_md):
			print(f"Target {wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol} already finalized locally. Skipping HPC pull.")
			return


		# Setup Logger
		log_path = params.log_abs
		os.makedirs(os.path.dirname(log_path), exist_ok=True)

		logger_name = f"finalize_custom_md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}"
		logger = logging.getLogger(logger_name)
		logger.setLevel(logging.INFO)
		logger.handlers.clear()

		fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

		fh = logging.FileHandler(log_path, mode="w")
		fh.setFormatter(fmt)
		logger.addHandler(fh)

		ch = logging.StreamHandler()
		ch.setFormatter(fmt)
		logger.addHandler(ch)

		target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
		prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}_md"

		logger.info("Finalizing custom MD output artifacts for target: %s", target_id)

		try:
			protocol_dir = os.path.abspath(f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}")
			job_dir = os.path.join(protocol_dir, "JOB")
			results_dir = os.path.abspath(os.path.dirname(output.done))
			os.makedirs(results_dir, exist_ok=True)

			search_dirs = [protocol_dir, job_dir]

			def find_and_copy(ext):
				target_filename = f"{prefix}{ext}"
				dest_path = os.path.join(results_dir, target_filename)

				for s_dir in search_dirs:
					src_path = os.path.join(s_dir, target_filename)
					if os.path.exists(src_path):
						logger.info("Found %s -> Copying to %s", src_path, dest_path)
						shutil.copy2(src_path, dest_path)
						return dest_path
				return None

			# Copy mandatory files
			copied_xtc = find_and_copy(".xtc")
			copied_tpr = find_and_copy(".tpr")

			# Copy auxiliary files if available (.gro, .cpt, .edr, .log)
			for ext in [".gro", ".cpt", ".edr", ".log"]:
				find_and_copy(ext)

			if not copied_xtc or not os.path.exists(output.xtc_md):
				raise FileNotFoundError(f"Required trajectory file '{prefix}.xtc' was not found in {search_dirs}")

			if not copied_tpr or not os.path.exists(output.tpr_md):
				raise FileNotFoundError(f"Required topology file '{prefix}.tpr' was not found in {search_dirs}")

			# Write sentinel file
			with open(output.done, "w") as f:
				f.write(f"Custom MD protocol '{wildcards.protocol}' successfully finalized for {target_id}.\n")

			logger.info("Finalization complete for %s. All output artifacts in %s", target_id, results_dir)

		except Exception as err:
			logger.exception("Execution failed in finalize_custom_md for %s: %s", target_id, str(err))
			raise

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
