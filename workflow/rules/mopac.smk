rule create_frame_mopac_hof_calc_job:
	input:
		em_framestar = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz",
		em_completed = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/em_completed.txt"
	output:
		job_description = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/mopac_job_description.txt"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/mopac/create_mopac_job.log"
	shell:
		"python scripts/smk/mopac/create_mopac_hof_job.py {input.em_framestar} {output.job_description} > {log} 2>&1"


rule run_mopac_hof_calc:
	input:
		job_description = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/mopac_job_description.txt"
	output:
		heat_of_formation = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv",
		calc_tar = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/mopac_results.tar.gz"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/mopac/run_mopac_hof_calc.log"
	shell:
		"python scripts/smk/mopac/run_mopac_hof_calc.py {input.job_description} {output.heat_of_formation} > {log} 2>&1"