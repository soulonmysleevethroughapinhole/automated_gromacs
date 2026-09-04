#!/usr/bin/env python3
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt


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


def plot_trajectory_dynamics(df: pd.DataFrame, pdb: str, source: str, model_id: str, out_png: str):
    """Plots C-alpha RMSD and Heat of Formation over simulation time on dual y-axes."""
    fig, ax1 = plt.subplots(figsize=(10, 5))

    color_rmsd = '#1f77b4'  # Steel Blue
    ax1.set_xlabel('Time [ns]', fontsize=12)
    ax1.set_ylabel(r'C$\alpha$ RMSD ($\AA$)', color=color_rmsd, fontsize=12)
    ax1.plot(df['time_ns'], df['RMSD'], color=color_rmsd, label=r'C$\alpha$ RMSD', linewidth=1.8)
    ax1.tick_params(axis='y', labelcolor=color_rmsd)
    ax1.grid(True, linestyle='--', alpha=0.3)

    ax2 = ax1.twinx()
    color_hof = '#d62728'  # Crimson Red
    ax2.set_ylabel(r'$\Delta_f H$ [kJ/mol]', color=color_hof, fontsize=12)
    ax2.plot(df['time_ns'], df['hof'], color=color_hof, label=r'$H_f$', linewidth=1.8, linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color_hof)

    plt.title(f'Trajectory Dynamics: {pdb} ({source} - {model_id})', fontsize=13, fontweight='bold')

    # Combined Legend
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', frameon=True)

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"✅ Saved dynamics plot: {out_png}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python plot_run_hof_n_rmsd_vs_time.py <rmsd_csv> <hof_csv> <output_png> [pdb] [source] [model_id]")
        sys.exit(1)

    rmsd_csv = sys.argv[1]
    hof_csv = sys.argv[2]
    out_png = sys.argv[3]

    pdb = sys.argv[4] if len(sys.argv) > 4 else "PDB"
    source = sys.argv[5] if len(sys.argv) > 5 else "Source"
    model_id = sys.argv[6] if len(sys.argv) > 6 else "Model"

    df = load_and_merge_data(rmsd_csv, hof_csv)
    plot_trajectory_dynamics(df, pdb, source, model_id, out_png)


if __name__ == "__main__":
    main()