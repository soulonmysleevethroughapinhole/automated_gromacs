#!/usr/bin/env python3
import sys
import os
import pandas as pd

def compute_relative_and_delta_hof(single_baseline_csv, frames_csv, output_csv):
	"""
	Calculates Delta H and Relative HoF normalized to canonical EM baseline.
	
	Formulas:
	  - dh (Delta H) = HOF@frame - HOF@em_canonical
	  - relative_hof = dh / HOF@em_canonical
	"""
	if not os.path.exists(single_baseline_csv):
		raise FileNotFoundError(f"Baseline HoF CSV not found: {single_baseline_csv}")
	if not os.path.exists(frames_csv):
		raise FileNotFoundError(f"Frames HoF CSV not found: {frames_csv}")

	# 1. Load Data
	df_base = pd.read_csv(single_baseline_csv)
	df_frames = pd.read_csv(frames_csv)

	if df_base.empty:
		raise ValueError(f"Baseline CSV is empty: {single_baseline_csv}")
	if df_frames.empty:
		raise ValueError(f"Frames CSV is empty: {frames_csv}")

	# 2. Extract Canonical Baseline Values (HOF@em_canonical)
	if "heat_of_formation_kcal_mol" in df_base.columns:
		hof_canonical_kcal = df_base["heat_of_formation_kcal_mol"].iloc[0]
		hof_canonical_kj = df_base["heat_of_formation_kj_mol"].iloc[0]
	else:
		raise KeyError("Could not find 'heat_of_formation_kcal_mol' in baseline CSV.")

	if hof_canonical_kcal == 0:
		raise ZeroDivisionError("HOF@em_canonical is 0, cannot divide by zero for relative_hof calculation.")

	# 3. Calculate dh = HOF@frame - HOF@em_canonical
	df_frames["dh_kcal_mol"] = df_frames["heat_of_formation_kcal_mol"] - hof_canonical_kcal
	df_frames["dh_kj_mol"] = df_frames["heat_of_formation_kj_mol"] - hof_canonical_kj

	# 4. Calculate relative_hof = dh / HOF@em_canonical
	df_frames["relative_hof"] = df_frames["dh_kcal_mol"] / hof_canonical_kcal
	df_frames["relative_hof_percent"] = df_frames["relative_hof"] * 100.0

	# Include metadata reference column
	df_frames["hof_em_canonical_kcal_mol"] = hof_canonical_kcal
	df_frames["hof_em_canonical_kj_mol"] = hof_canonical_kj

	# 5. Export Output
	os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
	df_frames.to_csv(output_csv, index=False)

	print(f"✅ Canonical Baseline HOF@em_canonical: {hof_canonical_kcal:.4f} kcal/mol")
	print(f"✅ Successfully computed dh and relative_hof datasets -> {output_csv}")

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python calc_hof_relative_and_delta.py <single_baseline_csv> <frames_csv> <output_csv>")
		sys.exit(1)

	single_csv = sys.argv[1]
	frames_csv = sys.argv[2]
	out_csv = sys.argv[3]

	try:
		compute_relative_and_delta_hof(single_csv, frames_csv, out_csv)
	except Exception as e:
		print(f"❌ Error computing relative/delta HoF: {str(e)}", file=sys.stderr)
		sys.exit(1)