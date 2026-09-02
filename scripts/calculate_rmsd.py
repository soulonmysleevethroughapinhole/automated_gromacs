#!/usr/bin/env python3
import sys
import os
import re
import glob
import tarfile
import tempfile
import pandas as pd
import pymol
from pymol import cmd

def extract_frame_num(filepath):
    """
    Extracts the integer frame index from filenames like 'frame_0.pdb' or 'frame_100.pdb'.
    Falls back to regex integer extraction if naming varies.
    """
    filename = os.path.basename(filepath)
    match = re.search(r'\d+', filename)
    if match:
        return int(match.group(0))
    return -1

def main():
    if len(sys.argv) < 4:
        print("Usage: python calculate_rmsd.py <FRAMES_compressed.tar.gz> <canonical_structure.pdb> <output_rmsd.csv>")
        sys.exit(1)

    tar_path = os.path.abspath(sys.argv[1])
    canonical_pdb = os.path.abspath(sys.argv[2])
    output_csv = os.path.abspath(sys.argv[3])

    if not os.path.exists(tar_path):
        print(f"Error: Trajectory archive not found at {tar_path}")
        sys.exit(1)

    if not os.path.exists(canonical_pdb):
        print(f"Error: Canonical reference PDB not found at {canonical_pdb}")
        sys.exit(1)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # 1. Unpack compressed PDB frames into a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Extracting {tar_path} into temporary workspace...")
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=temp_dir)

        # Locate extracted PDB files
        pdb_files = glob.glob(os.path.join(temp_dir, "**", "*.pdb"), recursive=True)
        # Exclude reference structure if accidentally packed inside
        pdb_files = [f for f in pdb_files if os.path.abspath(f) != canonical_pdb]

        if not pdb_files:
            print(f"Error: No frame PDB files found inside {tar_path}")
            sys.exit(1)

        # Pair each file with its numerical index and sort ascending
        framed_pdbs = [(f, extract_frame_num(f)) for f in pdb_files]
        framed_pdbs.sort(key=lambda x: x[1])

        print(f"Found {len(framed_pdbs)} frames. First frame: {framed_pdbs[0][1]}, Last frame: {framed_pdbs[-1][1]}")

        # 2. Initialize Headless PyMOL Engine
        pymol.finish_launching(['pymol', '-qc'])  # -q: quiet, -c: command-line
        cmd.reinitialize()

        # Load reference canonical structure
        ref_obj = "canonical_ref"
        cmd.load(canonical_pdb, ref_obj)

        data_rms = []

        # 3. Calculate C-alpha RMSD frame by frame
        for pdb_path, frame_num in framed_pdbs:
            frame_obj = f"frame_{frame_num}"
            cmd.load(pdb_path, frame_obj)

            # cmd.align returns a tuple: (RMSD_after_refinement, n_atoms, n_cycles, RMSD_before_refinement, ...)
            # Aligning C-alpha atoms between reference and trajectory frame
            align_result = cmd.align(f"{frame_obj} and name CA", f"{ref_obj} and name CA")
            rmsd_val = align_result[0]

            # Store tuple: (frame_index, rmsd)
            data_rms.append((frame_num, rmsd_val))

            cmd.delete(frame_obj)

        cmd.delete(ref_obj)
        cmd.quit()

    # 4. Format into DataFrame and export CSV
    df = pd.DataFrame(data_rms, columns=['frame', 'RMSD_CA'])
    
    # Strictly enforce sorting by numerical frame index (0, 1, ..., 100)
    df = df.sort_values('frame').reset_index(drop=True)

    # Optional: If 100 frames span 100 ns, compute time_ns explicitly
    # Assumes uniform spacing if total frames match protocol steps
    df['time_ns'] = df['frame'] * (100.0 / max(1, len(df) - 1)) if len(df) > 1 else 0.0

    df.to_csv(output_csv, index=False)
    print(f"Successfully saved RMSD series ({len(df)} entries) to {output_csv}")

if __name__ == "__main__":
    main()