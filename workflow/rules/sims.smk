import os
import shutil
import tarfile
import yaml
from dotenv import load_dotenv
load_dotenv() # Load environment variables from .env file

HPC_HOST = os.getenv("HPC_HOST")
HPC_USER = os.getenv("HPC_USER")
HPC_REMOTE_BASE = os.getenv("HPC_REMOTE_BASE")
HPC_PARTITION = os.getenv("HPC_PARTITION")
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
		pdb_config = config.get("custom_simulations", {}).get(wildcards.pdb, {})
		compute_target = pdb_config.get("compute_target", def_compute_target)
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
		#gro_cg = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		#top = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
		#mdp = "config/gromacs_settings/interruptable_config_ultimate/md.mdp",
		scheduling = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/scheduling.yml"
	output:
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job"
	log:
		"logs/{pdb}/{source}/{model_id}/create_md_job.log",
	params:
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/create_md_job.log"
		)
	run:

		with open(input.scheduling, 'r') as f:
			scheduling_info = yaml.load(f, Loader=yaml.FullLoader)

		if scheduling_info.get("COMPUTE") == "local":
			print(f"🐍 [Snakemake] No scheduling for local jobs - {wildcards.pdb} {wildcards.source} {wildcards.model_id} to run locally...🐍")
			with open(output.job_description, 'w') as f:
				f.write('BLANK')

		elif scheduling_info.get("COMPUTE") == "HPC":
			submit_hpc = SUBMIT_HPC == '1'

			if not submit_hpc:
				print(f"🐍 [Snakemake] HPC submission disabled. Job for {wildcards.pdb} {wildcards.source} {wildcards.model_id} will be created but not submitted.🐍")
				raise ValueError("HPC submission is disabled. Set SUBMIT_HPC=1 in .env to enable submission.")

			print(f"🐍 [Snakemake] Creating job for {HPC_USER} - {wildcards.pdb} {wildcards.source} {wildcards.model_id} to run on HPC...🐍")

			job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
			script_path = os.path.abspath(output.job_description)

			deffnm = f'{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md'
			content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}_%j.out
#SBATCH --partition={HPC_PARTITION}
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time={HPC_TIME}

gmx mdrun -v -deffnm {deffnm} -ntmpi 1 -nb gpu -pme gpu
	"""

			with open(script_path, 'w') as fh:
				fh.write(content)

			# Compress the JOB directory into a tarball for transfer to the HPC
			job_dir = os.path.dirname(output.job_description)
			job_targz = f'{job_dir}.tar.gz'
			with tarfile.open(job_targz, "w:gz") as tar:
				for fn in os.listdir(job_dir):
					p = os.path.join(job_dir, fn)
					tar.add(p, arcname=os.path.basename(fn))

			remote_dir = os.path.join(HPC_REMOTE_BASE, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns")
			ssh_target = f"{HPC_USER}@{HPC_HOST}"

			# Require successful transfer before writing the submitted status.
			shell("""
				ssh {ssh_target} "mkdir -p '{remote_dir}'"
				rsync -avz {job_targz} {ssh_target}:{remote_dir}/JOB.tar.gz
				ssh {ssh_target} "test -f '{remote_dir}/JOB.tar.gz' && echo 'HPC transfer confirmed: JOB.tar.gz exists at {remote_dir}/JOB.tar.gz'"
			""")

			with open(input.scheduling, 'a') as f:
				f.write("JOB_STATUS: Submitted\n")
				f.write(f"REMOTE_DIR: {remote_dir}\n")

		else:
			raise ValueError(f"Unknown compute target in scheduling info: {scheduling_info.get('COMPUTE')}")



		#os.chmod(script_path, 0o755)



# Rule 5: Run 100n Molecular Dynamics
rule run_molecular_dynamics:
	input:
		tpr_file = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md.tpr",
		job_description = "results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/JOB/{pdb}_{source}_{model_id}_md_job.job",
	output:
		# sentinel method (marks completion in the md_100ns folder)
		done = 'results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/md_completed.txt'
	log:
		"logs/{pdb}/{source}/{model_id}/production_mdrun.log"
	params:
		#log_abs = lambda wildcards, input, output, log: os.path.abspath(str(log[0]))
		#log_abs = lambda log: os.path.abspath(str(log))
		log_abs = lambda wildcards: os.path.abspath(
			f"logs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/production_mdrun.log"
		)

		#log_abs = lambda wildcards, input, log: os.path.abspath(str(log))
	resources:
		gpu = 1
	run:
		with open(input.job_description, 'r') as f:
			job_description = f.read().strip()

		if job_description == 'BLANK': # run locally 
			print(f"🐍 [Snakemake] Running MD locally for {wildcards.pdb} {wildcards.source} {wildcards.model_id}...🐍")

			#log_abs = os.path.abspath(str(log))
			#mdp_abs = os.path.abspath(str(input.mdp))
			job_dir = os.path.dirname(input.job_description)
			cpt_path = os.path.join(job_dir, f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md.cpt")
			prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
			
			# tpr compilation is handled earlier in the pipeline (schedule/create job)
			if os.path.exists(cpt_path):
				print(f"TODO: edit naming to be descriptive of current work (pdb + origin) 🐍 [Snakemake] Active checkpoint detected. Resuming {wildcards.pdb}...🐍")

				shell("""
					cd {job_dir}
					gmx mdrun -v -ntmpi 1 \
					-deffnm {prefix} \
					-cpi $(basename {cpt_path}) \
					-nb gpu -pme gpu >> {params.log_abs} 2>&1
				""")
			else:
				print(f"🐍 [Snakemake] Launching fresh execution tree for {wildcards.pdb} {wildcards.source} {wildcards.model_id}...🐍")

				shell("""
					cd {job_dir}
					gmx mdrun -v -ntmpi 1 \
					-deffnm {prefix} \
					-nb gpu -pme gpu >> {params.log_abs} 2>&1
				""")

			# symlinking I think!!! NOTE: be wary
			shell("""
				cd {job_dir}
				[ -f {prefix}.xtc ] && ln -sf {prefix}.xtc ../md.xtc
				[ -f {prefix}.tpr ] && ln -sf {prefix}.tpr ../md.tpr
			""")

			with open(output.done, "w") as f:
				f.write("MD simulation finished to completion.")
		else:
			print(f"🐍 [Snakemake] Submitting MD job to HPC for {wildcards.pdb} {wildcards.source} {wildcards.model_id}...🐍")

			ssh_target = f"{HPC_USER}@{HPC_HOST}"
			remote_dir = os.path.join(HPC_REMOTE_BASE, f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/standard_100ns")
			job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
			job_script = os.path.basename(input.job_description)
			prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
			local_dir = os.path.abspath(os.path.dirname(output.done))			#tar_gz = os.path.join(REMOTE_, f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md.tar.gz")

			shell("""
				SSH_TARGET="{ssh_target}"
				REMOTE_DIR="{remote_dir}"
				JOB_NAME="{job_name}"
				JOB_SCRIPT="{job_script}"
				PREFIX="{prefix}"
				LOCAL_DIR="{local_dir}"

				# 1. Check if the job finished previously on HPC (.gro coordinate output exists)
				REMOTE_DONE=$(ssh "$SSH_TARGET" "test -f '$REMOTE_DIR/$PREFIX.gro' && echo 'YES' || echo 'NO'")

				if [ "$REMOTE_DONE" = "YES" ]; then
					echo "🐍 [Snakemake] MD simulation already completed on remote HPC." >> "{params.log_abs}"
				else
					# 2. Check if job is currently active in Slurm queue
					JOB_ID=$(ssh "$SSH_TARGET" "squeue -u {HPC_USER} -n $JOB_NAME -h -o %i")

					if [ -n "$JOB_ID" ]; then
						echo "🐍 [Snakemake] Job $JOB_NAME (ID: $JOB_ID) is already in squeue. Monitoring..."
						# Poll every 30 seconds until job leaves queue
						ssh "$SSH_TARGET" "while squeue -j $JOB_ID -h | grep -q $JOB_ID; do sleep 30; done" >> "{params.log_abs}" 2>&1
					else
						echo "🐍 [Snakemake] Submitting new Slurm job to HPC..."
						# Extract transfer bundle and submit with --wait flag
						ssh "$SSH_TARGET" "cd '$REMOTE_DIR' && [ -f JOB.tar.gz ] && tar -xzf JOB.tar.gz && sbatch --wait $JOB_SCRIPT" >> "${params.log_abs}" 2>&1
					fi
				fi

				# 3. Tar up output trajectory, energy, and coordinate files on HPC
				echo "🐍 [Snakemake] Archiving MD results on HPC..." >> "{params.log_abs}"
				ssh "$SSH_TARGET" "cd '$REMOTE_DIR' && tar -czf md_results.tar.gz $PREFIX.*" >> "{params.log_abs}" 2>&1

				# 4. Sync tarball back to local project directory
				echo "🐍 [Snakemake] Transferring compressed results back to workstation..." >> "{params.log_abs}"
				rsync -avz "$SSH_TARGET:$REMOTE_DIR/md_results.tar.gz" "$LOCAL_DIR/" >> "{params.log_abs}" 2>&1

				# 5. Extract results locally and map symlinks for step 6 (trjconv)
				tar -xzf "$LOCAL_DIR/md_results.tar.gz" -C "$LOCAL_DIR/"
				
				if [ -f "$LOCAL_DIR/$PREFIX.xtc" ]; then
					ln -sf "$PREFIX.xtc" "$LOCAL_DIR/md.xtc"
				fi
				if [ -f "$LOCAL_DIR/$PREFIX.tpr" ]; then
					ln -sf "$PREFIX.tpr" "$LOCAL_DIR/md.tpr"
				fi
			""")

			with open(output.done, "w") as f:
				f.write(f"HPC MD simulation completed and output fetched successfully from {HPC_HOST}.")
			#job_dir_targz = f'{params.remote_dir}/JOB.tar.gz'



			# --- REMOTE HPC EXECUTION ---
			#shell("""				
			#	# 0. Check if targz already on server
			#
			#	# 1. Ensure clean remote directory and sync tarball via multiplexed SSH
			#	#ssh {ssh_target} "mkdir -p {params.remote_dir}"
			#	#rsync -avz {input.bundle} komondor:{params.remote_dir}/JOB.tar.gz
			#
			#	# 2. Extract, submit job, and poll for completion remotely
			#	ssh {ssh_target} "cd {remote_dir} && tar -xzf {jobdir} && sbatch --wait btestls64.job"
			#
			#	# 2.5: tar up results
			#	# 3. Download results back to workstation
			#	#rsync -avz {ssh_target}:{params.remote_dir}/md_prod.xtc {output.xtc}
			#	#rsync -avz {ssh_target}:{params.remote_dir}/md_prod.edr {output.edr}
			#
			#	# 4. Immediate Remote Wipe (Keep scratch 100% clean)
			#	#ssh {ssh_target} "rm -rf {params.remote_dir}"
			#""")


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
