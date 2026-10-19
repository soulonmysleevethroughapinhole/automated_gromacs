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
		#calc_tar = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/mopac_results.tar.gz"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/mopac/run_mopac_hof_calc.log"
	shell:
		"python scripts/smk/mopac/run_mopac_hof_calc.py {input.job_description} {output.heat_of_formation} {wildcards.protocol} > {log} 2>&1"



rule stage_em_pdb_for_mopac:
	input:
		em_structure = "results/gmx_em/{pdb}/{source}/em_{model_id}.pdb",
	output:
		canonical_mop = "results/mopac/{pdb}/{source}/{model_id}/SINGLE/em_{model_id}.mop"
	log:
		"logs/{pdb}/{source}/{model_id}/mopac/stage_em_canonical_pdb_for_mopac.log"
	shell:
		"python scripts/smk/mopac/stage_em_canonical_pdb_for_mopac.py {input.em_structure} {output.canonical_mop} > {log} 2>&1"


rule calculate_hof_for_pdb:
	input:
		canonical_mop = "results/mopac/{pdb}/{source}/{model_id}/SINGLE/em_{model_id}.mop"
	output:
		heat_of_formation = "results/mopac/{pdb}/{source}/{model_id}/SINGLE/heat_of_formation.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/mopac/calculate_hof_for_mop.log"
	shell:
		"python scripts/smk/mopac/calculate_hof_for_mop.py {input.canonical_mop} {output.heat_of_formation} > {log} 2>&1"

rule hof_relative_and_delta:
	input:
		single_heat_of_formation = "results/mopac/{pdb}/{source}/{model_id}/SINGLE/heat_of_formation.csv",
		heat_of_formation = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv",
	output:
		hofcsv = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation_relative_and_delta.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/mopac/hof_relative_and_delta.log"
	shell:
		"""
		python scripts/smk/mopac/calc_hof_relative_and_delta.py \
			{input.single_heat_of_formation} \
			{input.heat_of_formation} \
			{output.hofcsv} > {log} 2>&1
		"""
