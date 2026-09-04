#!/usr/bin/env python3
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def load_and_merge_data(rmsd_csv: str, hof_csv: str) -> pd.DataFrame:
    """Loads RMSD and Heat of Formation CSVs and merges them on time_ns."""
    if not os.path.exists(rmsd_csv):
        raise FileNotFoundError(f"RMSD CSV not found: {rmsd_csv}")
    if not os.path.exists(hof_csv):
        raise FileNotFoundError(f"HOF CSV not found: {hof_csv}")

    df_rmsd = pd.read_csv(rmsd_csv)
    df_hof = pd.read_csv(hof_csv)

    # Standardize RMSD column name (supports 'rmsd_ca_A' or 'RMSD')
    if 'rmsd_ca_A' in df_rmsd.columns:
        df_rmsd = df_rmsd.rename(columns={'rmsd_ca_A': 'RMSD'})
    elif 'RMSD' not in df_rmsd.columns:
        raise KeyError(f"Neither 'rmsd_ca_A' nor 'RMSD' column found in {rmsd_csv}")

    if 'time_ns' not in df_rmsd.columns or 'time_ns' not in df_hof.columns:
        raise KeyError("Both CSV files must contain a 'time_ns' column for merging.")

    merged_df = pd.merge(df_hof[['time_ns', 'hof']], df_rmsd[['time_ns', 'RMSD']], on='time_ns')
    return merged_df.sort_values('time_ns').dropna()


def plot_correlation(df: pd.DataFrame, pdb: str, source: str, model_id: str, out_png: str):
    """Generates a scatterplot of RMSD vs Heat of Formation with linear regression & R-squared."""
    x = df['hof']
    y = df['RMSD']

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_squared = r_value**2

    plt.figure(figsize=(8, 6))
    sns.regplot(
        x=x, y=y, data=df,
        scatter_kws={'alpha': 0.65, 'color': '#1f77b4', 's': 30},
        line_kws={'color': '#d62728', 'linewidth': 2, 'label': 'Linear Fit'}
    )

    # R-squared Annotation Box
    plt.text(
        0.05, 0.93, f'$R^2 = {r_squared:.3f}$',
        transform=plt.gca().transAxes,
        fontsize=13, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.85)
    )

    plt.ylabel(r'C$\alpha$ Root Mean Square Deviation ($\AA$)', fontsize=12)
    plt.xlabel(r'Heat of Formation ($\text{kJ/mol}$)', fontsize=12)
    plt.title(f'RMSD vs. Heat of Formation: {pdb} ({source} - {model_id})', fontsize=13, fontweight='bold')

    plt.grid(True, linestyle='--', alpha=0.3)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"✅ Saved correlation plot ($R^2 = {r_squared:.4f}$): {out_png}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python plot_rmsd_vs_hof.py <rmsd_csv> <hof_csv> <output_png> [pdb] [source] [model_id]")
        sys.exit(1)

    rmsd_csv = sys.argv[1]
    hof_csv = sys.argv[2]
    out_png = sys.argv[3]

    pdb = sys.argv[4] if len(sys.argv) > 4 else "PDB"
    source = sys.argv[5] if len(sys.argv) > 5 else "Source"
    model_id = sys.argv[6] if len(sys.argv) > 6 else "Model"

    df = load_and_merge_data(rmsd_csv, hof_csv)
    plot_correlation(df, pdb, source, model_id, out_png)


if __name__ == "__main__":
    main()