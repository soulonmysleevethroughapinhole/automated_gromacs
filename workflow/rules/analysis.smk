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

#rule plot_rmsd_vs_ddh:


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
rule GMX_MM_energies:
	input:
		edr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/md_results/{get_md_file_prefix(wildcards)}.edr",
		tpr_md = lambda wildcards: f"results/gromacs/{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}/JOB/{get_md_file_prefix(wildcards)}.tpr"
	output:
		energies = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_energies.xvg"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/gmx_energies.log"
	shell:
		"""
		echo "9 10 11 12 13 0" | gmx energy \
			-s {input.tpr_md} \
			-f {input.edr_md} \
			-o {output.energies} > {log} 2>&1
		"""


def parse_mm_energy_xvg(xvg_file):
	import pandas as pd
	headers = []
	data = []
	
	with open(xvg_file, 'r') as f:
		for line in f:
			line = line.strip()
			# Capture column legend names dynamically from @ s0 legend "..."
			if line.startswith('@ s') and 'legend' in line:
				label = line.split('legend')[-1].replace('"', '').strip()
				headers.append(label)
			elif line.startswith(('@', '#')):
				continue
			elif line:
				parts = [float(x) for x in line.split()]
				if parts:
					data.append(parts)
					
	if not data:
		raise ValueError(f"No numeric data found in {xvg_file}")
		
	num_cols = len(data[0])
	# Assign header names or default if legends were missing
	if headers and len(headers) == num_cols - 1:
		columns = ['Time_ps'] + headers
	else:
		columns = ['Time_ps'] + [f"col_{i}" for i in range(1, num_cols)]
		
	return pd.DataFrame(data, columns=columns)


rule parse_MM_energies:
	input:
		energies = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_energies.xvg"
	output:
		parsed_energies = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_energies_parsed.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/{protocol}/parse_gmx_energies.log"
	run:
		import pandas as pd
		try:
			df_energies = parse_mm_energy_xvg(input.energies)
			df_energies.to_csv(output.parsed_energies, index=False)

			with open(log[0], "w") as l:
				l.write(f"Successfully parsed energies into {output.parsed_energies}\n")

		except Exception as e:
			with open(log[0], "w") as l:
				l.write(f"Error parsing energies: {str(e)}\n")
			raise e


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


rule gmx_em_canonical_str_energies:
	input:
		em_edr = "results/gmx_em/{pdb}/{source}/em_{model_id}.edr"
	output:
		xvg = "results/gmx_em/{pdb}/{source}/em_{model_id}_energies.xvg",
	log:
		"logs/{pdb}/{source}/{model_id}/gmx_em_canonical_str_energies.log"
	shell:
		"""
		echo "Potential 0" | gmx energy \
			-f {input.em_edr} \
			-o {output.xvg} > {log} 2>&1
		"""

def parse_em_energy_xvg(xvg_file):
	import pandas as pd
	headers = []
	data = []
	
	with open(xvg_file, 'r') as f:
		for line in f:
			line = line.strip()
			# Parse dynamic column legend labels (e.g., Potential)
			if line.startswith('@ s') and 'legend' in line:
				label = line.split('legend')[-1].replace('"', '').strip()
				headers.append(label)
			elif line.startswith(('@', '#')):
				continue
			elif line:
				parts = [float(x) for x in line.split()]
				if parts:
					data.append(parts)
					
	if not data:
		raise ValueError(f"No numeric data found in {xvg_file}")
		
	num_cols = len(data[0])
	if headers and len(headers) == num_cols - 1:
		columns = ['Time_ps'] + headers
	else:
		columns = ['Time_ps'] + [f"col_{i}" for i in range(1, num_cols)]
		
	return pd.DataFrame(data, columns=columns)


rule parse_em_canonical_str_energies:
	input:
		xvg = "results/gmx_em/{pdb}/{source}/em_{model_id}_energies.xvg"
	output:
		csv = "results/gmx_em/{pdb}/{source}/em_{model_id}_energies.csv"
	log:
		"logs/{pdb}/{source}/{model_id}/parse_em_canonical_str_energies.log"
	run:
		try:
			df = parse_em_energy_xvg(input.xvg)
			df.to_csv(output.csv, index=False)
			
			with open(log[0], "w") as l:
				l.write(f"Successfully parsed canonical baseline energies to {output.csv}\n")
		except Exception as e:
			with open(log[0], "w") as l:
				l.write(f"Error parsing EM energies: {str(e)}\n")
			raise e


# TODO:
# implement structure selection (so that less frames have to be QMd)
# rmsd to 0
# rmsd to ref separation
# better logic handling files for analysis
#
#
#

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
		"results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_energies_parsed.csv"

	output:
		touch("results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/analysis_complete.txt")