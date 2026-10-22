# NOTE:
# analysis.smk is temporarily kept frozen
# this is a layer for extracting the data

rule process_raw_gmx_data:
	input:
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.xvg",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.xvg",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.xvg",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.xvg",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp_num.xvg",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv",
		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_energies_parsed.csv"

	output:
		touch("data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_processed.txt")

rule process_raw_mopac_data:
    input:
        "data/mopac_processed/{pdb}/{source}/{model_id}/{protocol}/mopac_summary.csv",
        "data/mopac_processed/{pdb}/{source}/{model_id}/{protocol}/mopac_energies_parsed.csv"

    output:
        touch("data/mopac_processesd/{pdb}/{source}/{model_id}/{protocol}/mopac_processed.txt")

