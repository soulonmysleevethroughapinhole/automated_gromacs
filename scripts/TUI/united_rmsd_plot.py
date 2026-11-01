from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Dynamically locate all individual RMSD files
base_dir = Path("data/simulations/rmsd")
csv_files = list(base_dir.glob("*/*/*/*/pymol_results/frames/rmsd_to_frame_0.csv"))

if not csv_files:
	raise FileNotFoundError(f"No rmsd_to_frame_0.csv files found under {base_dir.resolve()}")

print(f"Found {len(csv_files)} RMSD files. Reading and aggregating...")

data_frames = []
for file_path in csv_files:
	parts = file_path.parts
	pdb = parts[3]
	source = parts[4]
	model_id = parts[5]
	protocol = parts[6]

	if protocol != "standard_100ns":
		continue

	df_single = pd.read_csv(file_path, usecols=["time_ns", "rmsd_ca_A"])
	df_single["pdb"] = pdb
	df_single["source"] = source
	df_single["model_id"] = model_id
	df_single["protocol"] = protocol
	
	data_frames.append(df_single)

if not data_frames:
	raise ValueError("No standard_100ns trajectories found among the CSV files.")

df_all = pd.concat(data_frames, ignore_index=True)

# 2. Extract last 10 frames per trajectory run
df_sorted = df_all.sort_values(by=["pdb", "source", "model_id", "time_ns"])
df_last10 = (
	df_sorted.groupby(["pdb", "source", "model_id"])
	.tail(10)
	.reset_index(drop=True)
)

# Sort PDBs alphabetically for consistent x-axis ordering
pdbs = sorted(df_last10["pdb"].unique())

# 3. High-contrast, distinct color palette
hue_order = ["empirical", "alphafold", "colabfold"]
palette = {
	"empirical": "#0F4C81",  # Deep Royal Navy
	"alphafold": "#00A896",  # Vibrant Cyber Teal
	"colabfold": "#8A2BE2",  # Electric Violet
}
existing_sources = [s for s in hue_order if s in df_last10["source"].unique()]

# 4. Render Plot
plt.figure(figsize=(14, 7), dpi=300)
sns.set_theme(style="white")

# Main Boxplot
ax = sns.boxplot(
	data=df_last10,
	x="pdb",
	y="rmsd_ca_A",
	order=pdbs,
	hue="source",
	hue_order=existing_sources,
	palette=palette,
	width=0.6,
	boxprops=dict(alpha=0.85, edgecolor="black", linewidth=1.2),
	medianprops=dict(color="black", linewidth=1.5),
	whiskerprops=dict(color="black", linewidth=1.2),
	capprops=dict(color="black", linewidth=1.2),
	showmeans=True,
	meanprops={
		"marker": "o",
		"markerfacecolor": "white",
		"markeredgecolor": "black",
		"markersize": 6,
		"markeredgewidth": 1.2,
	},
)

# Stripplot overlay for individual frame distributions
sns.stripplot(
	data=df_last10,
	x="pdb",
	y="rmsd_ca_A",
	order=pdbs,
	hue="source",
	hue_order=existing_sources,
	dodge=True,
	color="black",
	alpha=0.4,
	jitter=0.18,
	size=3.5,
	linewidth=0.3,
	edgecolor="black",
	legend=False,
	ax=ax,
)

# Add vertical separator lines between PDBs
for i in range(len(pdbs) - 1):
	ax.axvline(i + 0.5, color="#CCCCCC", linestyle="--", linewidth=1.2, zorder=0)

# Add horizontal gridlines for readability
ax.yaxis.grid(True, linestyle=":", alpha=0.6, color="#888888")
ax.set_axisbelow(True)

# Formatting
plt.title("Baseline Trajectory Stability (Terminal 10 Frames of 100ns)", fontsize=15, fontweight="bold", pad=15)
plt.xlabel("PDB Entry", fontsize=12, fontweight="bold", labelpad=10)
plt.ylabel(r"C-$\alpha$ RMSD from $t=0$ ($\AA$)", fontsize=12, fontweight="bold", labelpad=10)

# Legend styling
plt.legend(
	title="Source",
	title_fontsize="11",
	loc="upper left",
	bbox_to_anchor=(1.01, 1),
	frameon=True,
	facecolor="white",
	edgecolor="#CCCCCC",
)

plt.tight_layout()

# Save output image
out_file = "md_rmsd_shortcut_finalized_boxplots.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
print(f"Finalized boxplot saved to {out_file}")
plt.show()