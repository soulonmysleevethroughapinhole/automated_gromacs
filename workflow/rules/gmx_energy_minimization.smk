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
		tar = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/FRAMES_compressed.tar.gz",
		em_completed = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/em_completed.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/energy_minimization/run_HPC_gmx_em.log"
	shell:
		"python scripts/smk/gmx_em/run_HPC_gmx_em.py {input.job_description} {output.em_completed} > {log} 2>&1"