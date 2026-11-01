import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

# ---------------------------------------------------------
# 1. DATA COLLECTION
# ---------------------------------------------------------
base_dir = Path("data/simulations/rmsd")
print(f"🔍 Searching for CSV files in {base_dir.resolve()}...")

# Find all matching CSV files
# Pattern: {pdb}/{source}/{model}/{protocol}/pymol_results/frames/rmsd_to_frame_0.csv
csv_files = list(base_dir.glob("*/*/*/standard_100ns/pymol_results/frames/rmsd_to_frame_0.csv"))
print(f"✅ Found {len(csv_files)} simulation files.")

# List to hold our extracted data
data_records = []

for file in csv_files:
    # Extract PDB and Source from the path structure
    # relative_to removes the base_dir, leaving: ('6WEY', 'colabfold', 'model_1', ...)
    parts = file.relative_to(base_dir).parts
    pdb_id = parts[0]
    source = parts[1] # e.g., 'colabfold' or 'empirical'
    
    try:
        # Read the CSV
        df = pd.read_csv(file)
        
        # Ensure we sort by time just in case frames are out of order
        df = df.sort_values('time_ns')
        
        # Get the last 10 frames
        last_10_rmsd = df['rmsd_ca_A'].tail(10)
        
        # Add each of the 10 values as a separate row for seaborn
        for rmsd_val in last_10_rmsd:
            data_records.append({
                'PDB': pdb_id,
                'Source': source.capitalize(), # Make it look nice (e.g., 'Colabfold', 'Empirical')
                'RMSD': rmsd_val
            })
    except Exception as e:
        print(f"⚠️ Could not read {file}: {e}")

# Create the master DataFrame
plot_df = pd.DataFrame(data_records)

if plot_df.empty:
    print("❌ No data was extracted. Check the directory structure and CSV contents.")
else:
    print(f"✅ Successfully extracted {len(plot_df)} data points for plotting.")

    # ---------------------------------------------------------
    # 2. PLOTTING
    # ---------------------------------------------------------
    print("🎨 Generating box plot...")
    
    # Set a clean, academic style
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define a professional color palette (e.g., Blue for Empirical, Orange/Red for AI)
    # Seaborn will assign colors alphabetically by Source, so we can map them explicitly if needed
    palette = {"Empirical": "#1f77b4", "Colabfold": "#ff7f0e", "Alphafold": "#d62728"}
    # Fallback to default if sources have different names
    available_sources = plot_df['Source'].unique()
    custom_palette = {src: palette.get(src, sns.color_palette()[i]) for i, src in enumerate(available_sources)}

    # 1. Draw the Box Plot
    # linewidth=1.5 makes the box edges thick
    # fliersize=3 keeps outliers small so they don't dominate
    sns.boxplot(
        data=plot_df, 
        x='PDB', 
        y='RMSD', 
        hue='Source',
        palette=custom_palette,
        width=0.6,
        linewidth=1.5,
        fliersize=3,
        medianprops={"linewidth": 3, "color": "black"}, # THICK median line for visibility!
        ax=ax,
        zorder=2 # Ensure boxes are drawn above grid
    )

    # 2. Overlay the Raw Data Points (The "Thin Box" Savior)
    # This plots the actual 10 points per simulation, proving the data is tightly clustered
    sns.stripplot(
        data=plot_df,
        x='PDB',
        y='RMSD',
        hue='Source',
        palette=custom_palette,
        dodge=True, # Aligns points with their respective boxes
        jitter=0.15, # Slight horizontal spread so points don't perfectly overlap
        size=5,
        alpha=0.6,
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
        zorder=3 # Draw points ON TOP of the boxes
    )

    # ---------------------------------------------------------
    # 3. FORMATTING & STYLING
    # ---------------------------------------------------------
    ax.set_xlabel('Protein Target (PDB ID)', fontsize=13, fontweight='bold', labelpad=10)
    ax.set_ylabel('RMSD to Frame 0 (Å)', fontsize=13, fontweight='bold', labelpad=10)
    ax.set_title('Conformational Stability: Last 10 Frames of MD Simulations', 
                 fontsize=15, fontweight='bold', pad=15)

    # Tilt the X-axis annotations (PDB IDs) for readability
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)

    # Clean up the legend (remove duplicate hue labels if stripplot creates them)
    handles, labels = ax.get_legend_handles_labels()
    # Keep only the first N handles/labels where N is the number of unique sources
    n_sources = len(plot_df['Source'].unique())
    ax.legend(handles[:n_sources], labels[:n_sources], 
              title='Simulation Source', 
              title_fontsize=12, 
              fontsize=11, 
              loc='upper right', 
              frameon=True, 
              edgecolor='black')

    # Remove top and right spines for a cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    
    # Save and show
    output_file = 'rmsd_last_10_frames_boxplot.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ Plot successfully saved as '{output_file}'")
    
    plt.show()