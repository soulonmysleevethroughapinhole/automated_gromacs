#TODO: plddt vs rmsf
#TODO: pae - domain contact/sasa
#TODO: dH vs dRMSD
#TODO:
#TODO:
import pandas as pd
import sys

rule plot_run_hof_n_rmsd_vs_time:
	input:
		rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv",
		hof = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv"
	output:
		plot = "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_n_hof_vs_time.png"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/plot_rmsd_n_hof_vs_time.log"
	shell:
		"python scripts/smk/analysis/plot_run_hof_n_rmsd_vs_time.py {input.rmsd} {input.hof} {output.plot} {wildcards.pdb} {wildcards.source} {wildcards.model_id} > {log} 2>&1"


rule plot_rmsd_vs_hof:
	input:
		rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv",
		hof = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv"
	output:
		plot = "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_vs_hof.png"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/plot_rmsd_vs_hof.log"
	shell:
		"python scripts/smk/analysis/plot_rmsd_vs_hof.py {input.rmsd} {input.hof} {output.plot} {wildcards.pdb} {wildcards.source} {wildcards.model_id} > {log} 2>&1"


#rule absolute_analysis_pdb:
#	input:
#		"results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_n_hof_vs_time.png",
#		"results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_vs_hof.png",
#		"results/movies/{pdb}/{source}/{model_id}/{protocol}_md_trajectory.mov"
#		#"results/movies/{pdb}/{source}/{model_id}/{protocol}_md_trajectory.mov"
#	output:
#		touch("results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/analysis_complete.txt")

# Helper function to resolve filename stems across standard vs custom protocols
def get_md_file_prefix(wildcards):
	if wildcards.protocol == "standard_100ns":
		return f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
	return f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}_md"


rule GMX_rmsd:
	input:
		xtc_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.xtc",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		rmsd = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.log"
	shell:
		"""
		echo "Backbone Backbone" | gmx rms \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o {output.rmsd} > {log} 2>&1
		"""


rule GMX_rmsf:
	input:
		xtc_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.xtc",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		rmsf = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.log"
	shell:
		"""
		echo "C-alpha" | gmx rmsf \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o {output.rmsf} \
			-res > {log} 2>&1
		"""


rule GMX_gyrate:
	input:
		xtc_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.xtc",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		gyrate = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.log"
	shell:
		"""
		echo "Protein" | gmx gyrate \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o {output.gyrate} > {log} 2>&1
		"""


rule GMX_sasa:
	input:
		xtc_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.xtc",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		sasa = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.log"
	shell:
		"""
		echo "Protein" | gmx sasa \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o {output.sasa} > {log} 2>&1
		"""

rule GMX_dssp:
	input:
		xtc_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.xtc",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		dssp_dat = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp.dat",
		dssp_num = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp_num.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp.log"
	shell:
		"""
		echo "MainChain" | gmx dssp \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o {output.dssp_dat} \
			-num {output.dssp_num} > {log} 2>&1 || \
		echo "Protein" | gmx do_dssp \
			-s {input.tpr_md} \
			-f {input.xtc_md} \
			-o results/analysis/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/gmx_dssp.xpm > {log} 2>&1
		"""

def parse_xvg(filepath, value_col_name="val"):
	data = []
	with open(filepath, "r") as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith(("#", "@")):
				continue
			parts = [float(x) for x in line.split()]
			if len(parts) >= 2:
				data.append(parts[:2])  # grabs time/idx and first primary metric
	df = pd.DataFrame(data, columns=["time_ps", value_col_name])
	return df


rule GMX_analysis:
	input:
		rmsd = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.xvg",
		rmsf = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.xvg",
		gyrate = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.xvg",
		sasa = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.xvg",
		dssp_num = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp_num.xvg"
	output:
		csv = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.log"
	run:
		import pandas as pd

		def parse_xvg(filepath, value_col_name="val"):
			data = []
			with open(filepath, "r") as f:
				for line in f:
					line = line.strip()
					if not line or line.startswith(("#", "@")):
						continue
					parts = [float(x) for x in line.split()]
					if len(parts) >= 2:
						data.append(parts[:2])
			
			df = pd.DataFrame(data, columns=["time_raw", value_col_name])
			
			# Detect if time is in ps (> 500) or ns and standardize to time_ns rounded to 3 decimals
			if df["time_raw"].max() > 500:
				df["time_ns"] = (df["time_raw"] / 1000.0).round(3)
			else:
				df["time_ns"] = df["time_raw"].round(3)
				
			return df[["time_ns", value_col_name]]

		try:
			df_rmsd = parse_xvg(input.rmsd, "rmsd_nm")
			df_gyrate = parse_xvg(input.gyrate, "gyrate_nm")
			df_sasa = parse_xvg(input.sasa, "sasa_nm2")
			df_dssp = parse_xvg(input.dssp_num, "dssp_structures")

			# Merge on standardized time_ns
			df_summary = df_rmsd.merge(df_gyrate, on="time_ns", how="outer")
			df_summary = df_summary.merge(df_sasa, on="time_ns", how="outer")
			df_summary = df_summary.merge(df_dssp, on="time_ns", how="outer")

			df_summary.sort_values("time_ns", inplace=True)
			df_summary.to_csv(output.csv, index=False)

			with open(log[0], "w") as l:
				l.write(f"Successfully compiled {len(df_summary)} time points into {output.csv}\n")

		except Exception as e:
			with open(log[0], "w") as l:
				l.write(f"Error parsing XVG files: {str(e)}\n")
			raise e

rule analyze_integrated_system:
	input:
		gmx = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv",
		rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv",
		hof = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv"
	output:
		csv = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/integrated_summary.csv",
		plot = "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/integrated_energy_landscape.png"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/integrated_analysis.log"
	shell:
		"""
		python scripts/smk/analysis/analyze_trajectory_energy.py \
			--gmx {input.gmx} \
			--rmsd {input.rmsd} \
			--hof {input.hof} \
			--out-csv {output.csv} \
			--out-plot {output.plot} \
			--title "{wildcards.pdb} ({wildcards.source}/{wildcards.model_id} - {wildcards.protocol})" > {log} 2>&1
		"""
rule absolute_analysis_pdb:
	input:
		"results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_n_hof_vs_time.png",
		"results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_vs_hof.png",
		"results/movies/{pdb}/{source}/{model_id}/{protocol}_md_trajectory.mov",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.xvg",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.xvg",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.xvg",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.xvg",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp_num.xvg",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv",
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/integrated_summary.csv"
	output:
		touch("results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/analysis_complete.txt")