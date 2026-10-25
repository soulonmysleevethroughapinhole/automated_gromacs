import os
import pandas as pd
import numpy as np

# NOTE:
# analysis.smk is temporarily kept frozen
# this is a layer for extracting the data

#rule process_raw_gmx_data:
#	input:
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsd.xvg",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_rmsf.xvg",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_gyrate.xvg",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_sasa.xvg",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_dssp_num.xvg",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv",
#		"data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_energies_parsed.csv"
#
#	output:
#		touch("data/gmx_processed/{pdb}/{source}/{model_id}/{protocol}/gmx_processed.txt")

#rule process_raw_mopac_data:
#    input:
#        "data/mopac_processed/{pdb}/{source}/{model_id}/{protocol}/mopac_summary.csv",
#        "data/mopac_processed/{pdb}/{source}/{model_id}/{protocol}/mopac_energies_parsed.csv"
#
#    output:
#        touch("data/mopac_processesd/{pdb}/{source}/{model_id}/{protocol}/mopac_processed.txt")

rule united_simulation_data:
	input:
		gmx_summary = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_summary.csv",
		gmx_energies_summary = "results/analysis/{pdb}/{source}/{model_id}/{protocol}/gmx_energies_parsed.csv",
		mopac_summary = "data/simulations/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation_relative_and_delta.csv",
		rmsd_to_em_canonical = "data/simulations/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd_to_em_empirical.csv",
		rmsd_to_frame_0 = "data/simulations/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd_to_frame_0.csv"
	output:
		united_data = "data/united/{pdb}/{source}/{model_id}/{protocol}/united_data.csv"
	run:
		os.makedirs(os.path.dirname(output.united_data), exist_ok=True)

		# ------------------------------------------------------------------
		# 1. Base Backbone: GROMACS Summary (time_ns anchor)
		# ------------------------------------------------------------------
		df_base = pd.read_csv(input.gmx_summary)
		df_base.columns = [c.strip().lower() for c in df_base.columns]

		if 'time_ns' not in df_base.columns:
			# Fallback if first col is time without explicit name
			df_base.rename(columns={df_base.columns[0]: 'time_ns'}, inplace=True)

		df_base['time_ns'] = df_base['time_ns'].astype(float).round(4)
		df_base = df_base.sort_values('time_ns').reset_index(drop=True)

		max_time_ns = df_base['time_ns'].max()

		# Helper: standardizes time across frame-based or ps-based dataframes
		def normalize_dataframe_time(df, prefix=""):
			df.columns = [c.strip().lower() for c in df.columns]

			# Detect time vs frame columns
			time_col = next((c for c in df.columns if 'time' in c), None)
			frame_col = next((c for c in df.columns if 'frame' in c or 'index' in c), None)

			if time_col:
				# Convert ps -> ns if values far exceed max simulation ns
				if df[time_col].max() > (max_time_ns * 10):
					df['time_ns'] = (df[time_col] / 1000.0).round(4)
				else:
					df['time_ns'] = df[time_col].astype(float).round(4)
				if time_col != 'time_ns':
					df.drop(columns=[time_col], inplace=True)
			elif frame_col:
				# Calculate time_ns linearly from frame index if time column missing
				num_frames = len(df)
				if num_frames > 1:
					df['time_ns'] = np.linspace(0.0, max_time_ns, num_frames).round(4)
				else:
					df['time_ns'] = 0.0
				df.drop(columns=[frame_col], inplace=True)
			else:
				# Fallback to index-based linear interpolation
				df['time_ns'] = np.linspace(0.0, max_time_ns, len(df)).round(4)

			# Prefix value columns to avoid collisions across files
			renames = {c: f"{prefix}_{c}" for c in df.columns if c != 'time_ns'}
			df.rename(columns=renames, inplace=True)
			return df.sort_values('time_ns').reset_index(drop=True)

		# ------------------------------------------------------------------
		# 2. Process Auxiliary DataSets
		# ------------------------------------------------------------------
		df_eng = normalize_dataframe_time(pd.read_csv(input.gmx_energies_summary), prefix="energy")
		df_mopac = normalize_dataframe_time(pd.read_csv(input.mopac_summary), prefix="mopac")
		df_rmsd_em = normalize_dataframe_time(pd.read_csv(input.rmsd_to_em_canonical), prefix="pymol_em")
		df_rmsd_f0 = normalize_dataframe_time(pd.read_csv(input.rmsd_to_frame_0), prefix="pymol_f0")

		# ------------------------------------------------------------------
		# 3. Robust Alignment via pd.merge_asof (Nearest Time Matching)
		# ------------------------------------------------------------------
		# merge_asof matches nearest timestamps, preventing floating point 
		# offset bugs (e.g. 1.0000001 ns vs 1.0 ns) from causing NaN rows.
		merged = pd.merge_asof(df_base, df_eng, on='time_ns', direction='nearest')
		merged = pd.merge_asof(merged, df_mopac, on='time_ns', direction='nearest')
		merged = pd.merge_asof(merged, df_rmsd_em, on='time_ns', direction='nearest')
		merged = pd.merge_asof(merged, df_rmsd_f0, on='time_ns', direction='nearest')

		# ------------------------------------------------------------------
		# 4. Attach Metadata Columns & Export
		# ------------------------------------------------------------------
		merged.insert(0, 'protocol', wildcards.protocol)
		merged.insert(0, 'model_id', wildcards.model_id)
		merged.insert(0, 'source', wildcards.source)
		merged.insert(0, 'pdb', wildcards.pdb)

		merged.to_csv(output.united_data, index=False)