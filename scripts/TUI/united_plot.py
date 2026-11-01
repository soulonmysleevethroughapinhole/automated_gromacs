from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Dynamically locate and concatenate all individual united_data.csv files
united_dir = Path("data/united")
csv_files = list(united_dir.glob("*/*/*/*/united_data.csv"))

if not csv_files:
	raise FileNotFoundError(
		f"No united_data.csv files found under {united_dir.resolve()}"
	)

print(f"Found {len(csv_files)} united_data.csv files. Concatenating...")
df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)

# 2. Filter for standard 100ns protocol
df_100ns = df[df["protocol"] == "standard_100ns"].copy()

# Sort chronologically per trajectory
df_100ns = df_100ns.sort_values(by=["pdb", "source", "model_id", "time_ns"])

# 3. Extract last 10 entries per trajectory run
df_last10 = (
	df_100ns.groupby(["pdb", "source", "model_id"])
	.tail(10)
	.reset_index(drop=True)
)

# 4. Plotting setup
hue_order = ["empirical", "alphafold", "colabfold"]
palette = {"empirical": "#2b5c8f", "alphafold": "#d95f02", "colabfold": "#7570b3"}

# Only keep sources present in the loaded dataset
existing_sources = [s for s in hue_order if s in df_last10["source"].unique()]

plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

# Main Boxplot
ax = sns.boxplot(
	data=df_last10,
	x="pdb",
	#y="rmsd_nm",
	y="pymol_f0_rmsd_ca_a",
	hue="source",
	hue_order=existing_sources,
	palette=palette,
	showmeans=True,
	meanprops={
		"marker": "o",
		"markerfacecolor": "white",
		"markeredgecolor": "black",
		"markersize": "6",
	},
)

# Overlay data points
sns.stripplot(
	data=df_last10,
	x="pdb",
	y="pymol_f0_rmsd_ca_a",
	hue="source",
	hue_order=existing_sources,
	dodge=True,
	color="black",
	alpha=0.35,
	jitter=0.2,
	size=4,
	legend=False,
)

plt.title(
	"Trajectory Stability Analysis (Terminal 10 Frames of 100ns Baseline)",
	fontsize=14,
	pad=15,
)
plt.xlabel("PDB Entry", fontsize=12, labelpad=10)
plt.ylabel("RMSD from $t=0$ (nm)", fontsize=12, labelpad=10)
plt.legend(title="Source", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

# Save & render
output_plot = "md_stability_rmsd_boxplots.png"
plt.savefig(output_plot, dpi=300)
print(f"Plot successfully saved to {output_plot}")
plt.show()