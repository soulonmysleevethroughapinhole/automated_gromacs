


rule unpack_frames_tar:
	input:
		em_tar = "results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz",
		raw_tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz",
	output:
		em_dir = directory("data/simulations/{pdb}/{source}/{model_id}/{protocol}/em_500fr_target/EM_FRAMES"),
		raw_dir = directory("data/simulations/{pdb}/{source}/{model_id}/{protocol}/md_500fr_target/RAW_FRAMES"),
		frame_0 = "data/simulations/{pdb}/{source}/{model_id}/{protocol}/em_500fr_target/EM_FRAMES/frame_0.pdb",
		frame_0_raw = "data/simulations/{pdb}/{source}/{model_id}/{protocol}/md_500fr_target/RAW_FRAMES/frame_0.pdb"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/unpack_frames_tar.log"
	shell:
		"""
		mkdir -p $(dirname {output.em_dir})
		mkdir -p $(dirname {output.raw_dir})
        tar -xzf {input.em_tar} -C $(dirname {output.em_dir}) > {log} 2>&1
        tar -xzf {input.raw_tar} -C $(dirname {output.raw_dir}) > {log} 2>&1
		"""

rule calculate_rmsd_to_em_empirical_CAN_STR:
	input:
		em_dir = "data/simulations/{pdb}/{source}/{model_id}/{protocol}/em_500fr_target/EM_FRAMES",
		#canonical_structure = "results/structures/{pdb}/empirical/canonical_structure.pdb"
		em_canonical_structure = "results/gmx_em/{pdb}/empirical/em_canonical_structure.pdb"
	output:
		rmsd_to_em_canonical = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd_to_em_empirical.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/calculate_rmsd_to_em_canonical.log"
	shell:
		"python scripts/smk/pymol/calculate_rmsd_to_ref.py {input.em_dir} {input.em_canonical_structure} {output.rmsd_to_em_canonical} {wildcards.protocol} > {log} 2>&1"

rule calculate_rmsd_to_em_SELF_frame_0:
	input:
		em_dir = "data/simulations/{pdb}/{source}/{model_id}/{protocol}/em_500fr_target/EM_FRAMES",
		frame_0 = "data/simulations/{pdb}/{source}/{model_id}/{protocol}/em_500fr_target/EM_FRAMES/frame_0.pdb"
	output:
		rmsd_to_frame_0 = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd_to_frame_0.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/calculate_rmsd_to_em_frame0.log"
	shell:
		"python scripts/smk/pymol/calculate_rmsd_to_ref.py {input.em_dir} {input.frame_0} {output.rmsd_to_frame_0} {wildcards.protocol} > {log} 2>&1"


# --- STEP 7: AUTOMATED PYMOL MOVIE GENERATION ---
rule generate_pymol_movie:
	input:
		tar = "results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	output:
		movie = "results/movies/{pdb}/{source}/{model_id}/{protocol}_md_trajectory.mov"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/pymol_render.log"
	shell:
		"python scripts/smk/pymol/render_pymol_movie.py {input.tar} {output.movie} > {log} 2>&1"
