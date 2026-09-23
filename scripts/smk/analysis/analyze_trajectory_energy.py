#!/usr/bin/env python3
"""
Integrated Trajectory Analysis & Energy Landscape Plotter
"""

import sys
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def setup_plotting_style():
    sns.set_theme(style="ticks", palette="muted")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "figure.dpi": 300
    })

def find_column_by_keywords(df, keywords, default_name):
    """Find a column in df matching any keyword, otherwise return None."""
    for col in df.columns:
        col_lower = str(col).lower()
        if any(kw in col_lower for kw in keywords):
            return col
    return None

def load_and_merge_data(gmx_summary_path, pymol_rmsd_path, mopac_hof_path):
    print(f"[*] Loading GMX summary: {gmx_summary_path}")
    df_gmx = pd.read_csv(gmx_summary_path)

    print(f"[*] Loading PyMOL frame RMSD: {pymol_rmsd_path}")
    df_pymol = pd.read_csv(pymol_rmsd_path)

    print(f"[*] Loading MOPAC Heat of Formation: {mopac_hof_path}")
    df_hof = pd.read_csv(mopac_hof_path)

    # Standardize column headers to lowercase
    df_pymol.columns = [str(c).lower().strip() for c in df_pymol.columns]
    df_hof.columns = [str(c).lower().strip() for c in df_hof.columns]

    # Handle frame column
    if "frame" not in df_pymol.columns:
        df_pymol["frame"] = df_pymol.index
    if "frame" not in df_hof.columns:
        df_hof["frame"] = df_hof.index

    # Locate RMSD and HoF columns dynamically
    rmsd_col = find_column_by_keywords(df_pymol, ["rmsd"], df_pymol.columns[1] if len(df_pymol.columns) > 1 else df_pymol.columns[0])
    hof_col = find_column_by_keywords(df_hof, ["heat", "hof", "energy", "formation", "kcal"], df_hof.columns[1] if len(df_hof.columns) > 1 else df_hof.columns[0])

    if not rmsd_col or not hof_col:
        raise ValueError(f"Could not identify RMSD or HoF columns. PyMOL cols: {list(df_pymol.columns)}, HoF cols: {list(df_hof.columns)}")

    df_pymol.rename(columns={rmsd_col: "pymol_rmsd"}, inplace=True)
    df_hof.rename(columns={hof_col: "heat_of_formation"}, inplace=True)

    # Merge frame data
    df_frames = pd.merge(df_pymol[["frame", "pymol_rmsd"]], df_hof[["frame", "heat_of_formation"]], on="frame", how="inner")

    # Map frame indices linearly onto time_ns from GMX summary
    total_frames = len(df_frames)
    max_time_ns = df_gmx["time_ns"].max() if "time_ns" in df_gmx.columns else df_gmx.iloc[:, 0].max()
    min_time_ns = df_gmx["time_ns"].min() if "time_ns" in df_gmx.columns else df_gmx.iloc[:, 0].min()

    df_frames["time_ns"] = np.linspace(min_time_ns, max_time_ns, total_frames).round(3)

    # Merge with GMX summary using nearest time point matching
    time_col_gmx = "time_ns" if "time_ns" in df_gmx.columns else df_gmx.columns[0]
    df_gmx_sorted = df_gmx.sort_values(time_col_gmx).copy()
    if time_col_gmx != "time_ns":
        df_gmx_sorted.rename(columns={time_col_gmx: "time_ns"}, inplace=True)

    df_merged = pd.merge_asof(
        df_frames.sort_values("time_ns"),
        df_gmx_sorted,
        on="time_ns",
        direction="nearest"
    )

    return df_merged

def calculate_correlations(df):
    metrics = [c for c in ["pymol_rmsd", "rmsd_nm", "heat_of_formation", "gyrate_nm", "sasa_nm2", "dssp_structures"] if c in df.columns]
    pearson_corr = df[metrics].corr(method="pearson")
    spearman_corr = df[metrics].corr(method="spearman")
    return pearson_corr, spearman_corr

def generate_analysis_plots(df, output_plot_path, system_title=""):
    setup_plotting_style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Integrated Trajectory Analysis: {system_title}", fontweight="bold", y=0.98)

    rmsd_val = "rmsd_nm" if "rmsd_nm" in df.columns else "pymol_rmsd"

    # Panel A: Time Evolution
    ax1 = axes[0, 0]
    ax1_twin = ax1.twinx()

    l1 = ax1.plot(df["time_ns"], df["heat_of_formation"], color="#d95f02", alpha=0.85, label="HoF (kcal/mol)")
    l2 = ax1_twin.plot(df["time_ns"], df[rmsd_val], color="#7570b3", alpha=0.85, label="RMSD (nm)")

    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("Heat of Formation (kcal/mol)", color="#d95f02")
    ax1_twin.set_ylabel("RMSD (nm)", color="#7570b3")
    ax1.set_title("A. Thermodynamic & Structural Drift")

    # Panel B: Energy Landscape
    ax2 = axes[0, 1]
    sns.scatterplot(data=df, x=rmsd_val, y="heat_of_formation", hue="time_ns", palette="viridis", ax=ax2, alpha=0.8)

    if "heat_of_formation" in df.columns and not df["heat_of_formation"].isna().all():
        min_row = df.loc[df["heat_of_formation"].idxmin()]
        ax2.annotate(
            f"Global Min\nFrame {int(min_row['frame'])}\n({min_row['heat_of_formation']:.1f} kcal/mol)",
            xy=(min_row[rmsd_val], min_row["heat_of_formation"]),
            xytext=(min_row[rmsd_val] + 0.05, min_row["heat_of_formation"] + 5),
            arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=5),
            fontsize=8, fontweight="bold", bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.6)
        )

    ax2.set_xlabel("RMSD (nm)")
    ax2.set_ylabel("Heat of Formation (kcal/mol)")
    ax2.set_title(r"B. Energy Landscape ($\Delta H_f$ vs. RMSD)")

    # Panel C: Structural Compactness vs SASA
    ax3 = axes[1, 0]
    if "gyrate_nm" in df.columns and "sasa_nm2" in df.columns:
        sns.regplot(data=df, x="gyrate_nm", y="sasa_nm2", ax=ax3, color="#1b9e77")
        ax3.set_xlabel("Radius of Gyration (nm)")
        ax3.set_ylabel(r"SASA ($nm^2$)")
        ax3.set_title("C. Compactness vs. Solvent Exposure")
    else:
        ax3.text(0.5, 0.5, "Gyrate/SASA Data Missing", ha="center", va="center")

    # Panel D: Correlation Heatmap
    ax4 = axes[1, 1]
    p_corr, s_corr = calculate_correlations(df)
    sns.heatmap(s_corr, annot=True, fmt=".2f", cmap="coolwarm", vmin=-1, vmax=1, ax=ax4)
    ax4.set_title("D. Spearman Rank Correlation Matrix")

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_plot_path)), exist_ok=True)
    plt.savefig(output_plot_path, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze GMX summary, PyMOL RMSD, and MOPAC HoF.")
    parser.add_argument("--gmx", required=True, help="Path to gmx_summary.csv")
    parser.add_argument("--rmsd", required=True, help="Path to pymol rmsd.csv")
    parser.add_argument("--hof", required=True, help="Path to MOPAC heat_of_formation.csv")
    parser.add_argument("--out-csv", dest="out_csv", required=True, help="Path to output merged CSV")
    parser.add_argument("--out-plot", dest="out_plot", required=True, help="Path to output PNG plot")
    parser.add_argument("--title", default="", help="System title for plot headers")

    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_plot)), exist_ok=True)

    df = load_and_merge_data(args.gmx, args.rmsd, args.hof)
    df.to_csv(args.out_csv, index=False)
    print(f"[+] Saved integrated dataset to: {args.out_csv}")

    p_corr, s_corr = calculate_correlations(df)
    print("\n=== Spearman Rank Correlations ===")
    print(s_corr.round(3))

    if "heat_of_formation" in df.columns and not df["heat_of_formation"].isna().all():
        min_idx = df["heat_of_formation"].idxmin()
        max_idx = df["heat_of_formation"].idxmax()
        rmsd_col = "rmsd_nm" if "rmsd_nm" in df.columns else "pymol_rmsd"
        print("\n=== Key Thermodynamic States ===")
        print(f"Global Minimum Energy Frame : Frame {int(df.loc[min_idx, 'frame'])} | HoF: {df.loc[min_idx, 'heat_of_formation']:.2f} kcal/mol | RMSD: {df.loc[min_idx, rmsd_col]:.3f} nm")
        print(f"Maximum Energy State Frame : Frame {int(df.loc[max_idx, 'frame'])} | HoF: {df.loc[max_idx, 'heat_of_formation']:.2f} kcal/mol | RMSD: {df.loc[max_idx, rmsd_col]:.3f} nm")

    generate_analysis_plots(df, args.out_plot, system_title=args.title)

if __name__ == "__main__":
    main()