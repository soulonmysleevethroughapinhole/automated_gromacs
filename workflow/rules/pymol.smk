

rule calculate_rmsd:
	input:
		#tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz",
		tar = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz",
		# TODO: switch to energy minimized inputs!
		#canonical_structure = "results/structures/{pdb}/empirical/canonical_structure.pdb"
		reference_structure = "results/gmx_em/{pdb}/empirical/em_canonical_structure.pdb"
	output:
		rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/calculate_rmsd.log"
	shell:
		"python scripts/smk/pymol/calculate_rmsd_to_ref.py {input.tar} {input.reference_structure} {output.rmsd} > {log} 2>&1"

