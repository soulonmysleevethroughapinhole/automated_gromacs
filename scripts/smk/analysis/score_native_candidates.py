#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd


def min_max_normalize(series: pd.Series) -> pd.Series:
	"""Min-Max normalizes a series to [0, 1]. Returns 0 if constant."""
	min_val = series.min()
	max_val = series.max()
	if max_val - min_val == 0:
		return pd.Series(0.0, index=series.index)
	return (series - min_val) / (max_val - min_val)


def score_trajectory(input_csv: str, output_scored_csv: str, output_top_csv: str, top_n: int = 5):
	df = pd.read_csv(input_csv)

	if df.empty:
		raise ValueError(f"Input dataset is empty: {input_csv}")

	# 1. Compute Intramolecular Energy (E_intra) [kJ/mol]
	bonded_cols = ["energy_angle", "energy_proper dih.", "energy_per. imp. dih."]
	prot_nonbonded_cols = [
		"energy_coul-14:protein-protein",
		"energy_lj-14:protein-protein",
		"energy_coul-sr:protein-protein",
		"energy_lj-sr:protein-protein",
	]

	df["e_intra_kj_mol"] = df[bonded_cols].sum(axis=1) + df[prot_nonbonded_cols].sum(axis=1)

	# 2. Compute Solvation Interaction Energy (E_solvation) [kJ/mol]
	solvation_cols = ["energy_coul-sr:protein-sol", "energy_lj-sr:protein-sol"]
	df["e_solvation_kj_mol"] = df[solvation_cols].sum(axis=1)

	# 3. Per-Trajectory Min-Max Normalization [0.0 - 1.0]
	df["norm_e_intra"] = min_max_normalize(df["e_intra_kj_mol"])
	df["norm_rg"] = min_max_normalize(df["gyrate_nm"])
	df["norm_rmsd_ca"] = min_max_normalize(df["pymol_em_rmsd_ca_a"])

	# 4. Native-Likeness Index (NLI) Calculation
	# Lower NLI = more native-like (lower energy, compact, low C-alpha RMSD drift)
	df["nli_score"] = (
		(0.50 * df["norm_e_intra"]) +
		(0.25 * df["norm_rg"]) +
		(0.25 * df["norm_rmsd_ca"])
	)

	# 5. Export Full Scored Dataset (ordered by time_ns)
	df_scored = df.sort_values("time_ns").reset_index(drop=True)
	out_scored_path = Path(output_scored_csv)
	out_scored_path.parent.mkdir(parents=True, exist_ok=True)
	df_scored.to_csv(out_scored_path, index=False)

	# 6. Extract Top Candidates (ordered by nli_score ascending)
	df_top = df.sort_values("nli_score").head(top_n).reset_index(drop=True)
	df_top["nli_rank"] = df_top.index + 1

	out_top_path = Path(output_top_csv)
	out_top_path.parent.mkdir(parents=True, exist_ok=True)
	df_top.to_csv(out_top_path, index=False)

	print(f"✅ Scored trajectory saved to -> {out_scored_path}")
	print(f"✅ Top {top_n} candidates saved to -> {out_top_path}")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python score_native_candidates.py <input_csv> <output_scored_csv> <output_top_csv> [top_n]")
		sys.exit(1)

	input_csv = sys.argv[1]
	out_scored = sys.argv[2]
	out_top = sys.argv[3]
	n_candidates = int(sys.argv[4]) if len(sys.argv) > 4 else 5

	score_trajectory(input_csv, out_scored, out_top, top_n=n_candidates)