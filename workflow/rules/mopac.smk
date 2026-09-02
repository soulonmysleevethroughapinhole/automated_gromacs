

rule create_frame_hof_calc:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	output:
		heat_of_formation = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/heat_of_formation.txt"

# TODO: array job on hpc 
