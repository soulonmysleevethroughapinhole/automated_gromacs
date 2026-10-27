#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd

def compute_relative_and_delta_hof(single_csv: str, frames_csv: str, output_csv: str):
	df_em_empirical = pd.read_csv(single_csv)
	df_frames = pd.read_csv(frames_csv)

	# 1. Extract Baseline HOF as a scalar float (kJ/mol)
	em_empirical_hof = float(df_em_empirical["heat_of_formation_kj_mol"].iloc[0])
	frame_0_hof = float(df_frames["hof"].iloc[0])
	#if em_empirical_hof == 0:
	#	raise ZeroDivisionError("Baseline HOF is 0 kJ/mol, cannot calculate relative HoF.")

	# 2. Attach Baseline Metadata
	df_frames["em_canonical_hof_metadata"] = em_empirical_hof
	df_frames["frame_0_hof_kj_mol"] = frame_0_hof

	# 3. Calculate Delta H & Relative HOF
	# Canonical Relative Delta Hf: Hf(t) - Hf(em_empirical)
	df_frames["canonical_relative_dH_f"] = df_frames["hof"] - em_empirical_hof

	# Trajectory Relative Delta Hf: Hf(t) - Hf(t=0)
	df_frames["traj_relative_dH_f"] = df_frames['hof'] - frame_0_hof

	#df_frames["relative_hof"] = df_frames["dh_kj_mol"] / base_kj
	#df_frames["relative_hof_percent"] = df_frames["relative_hof"] * 100.0

	# 5. Export
	out_path = Path(output_csv)
	out_path.parent.mkdir(parents=True, exist_ok=True)
	df_frames.to_csv(out_path, index=False)
	
	print(f"✅ empirical HOF: {em_empirical_hof:.2f} kJ/mol | Output -> {out_path}")

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python calc_hof_relative_and_delta.py <baseline_csv> <frames_csv> <output_csv>")
		sys.exit(1)

	compute_relative_and_delta_hof(sys.argv[1], sys.argv[2], sys.argv[3])