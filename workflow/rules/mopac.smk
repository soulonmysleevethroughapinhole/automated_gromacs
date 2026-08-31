

rule calculate_heat_of_formation:
    input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
