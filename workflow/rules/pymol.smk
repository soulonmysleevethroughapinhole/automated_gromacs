

rule calculate_rmsd:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz",
		# TODO: switch to energy minimized inputs!
		canonical_structure = "results/structures/{pdb}/empirical/canonical_structure.pdb"
	output:
		rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/rmsd.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/calculate_rmsd.log"
	shell:
		"python scripts/calculate_rmsd.py {input.tar} {input.canonical_structure} {output.rmsd} > {log} 2>&1"