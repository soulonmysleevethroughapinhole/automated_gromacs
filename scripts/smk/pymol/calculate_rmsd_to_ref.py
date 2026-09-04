#!/usr/bin/env python3
import sys
import os
import re
import glob
import tarfile
import tempfile
import pandas as pd
from pathlib import Path

# Initialize PyMOL in headless/quiet mode before importing cmd
import pymol
pymol.finish_launching(['pymol', '-cqk'])
from pymol import cmd


def extract_frame_num(filename: str) -> int:
    """Extracts integer frame/time index from strings like 'frame_10.pdb' or 'frame10.pdb'."""
    match = re.search(r'\d+', os.path.basename(filename))
    return int(match.group(0)) if match else -1


def calculate_rmsd_pymol(tar_path: str, ref_pdb_path: str, output_csv: str):
    if not os.path.exists(tar_path):
        raise FileNotFoundError(f"Input frame tar archive not found: {tar_path}")
    if not os.path.exists(ref_pdb_path):
        raise FileNotFoundError(f"Reference canonical PDB structure not found: {ref_pdb_path}")

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)

    # Re-initialize PyMOL session cleanly
    cmd.reinitialize()

    # Load reference structure
    ref_obj = "reference_structure"
    cmd.load(ref_pdb_path, object=ref_obj)

    data_rms = []

    with tempfile.TemporaryDirectory() as temp_dir:
        # Extract frame archive safely
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=temp_dir, filter='data')

        # Locate frame PDBs (supports nested EM_FRAMES folder structure)
        frame_pdbs = sorted(
            glob.glob(os.path.join(temp_dir, "**", "*.pdb"), recursive=True),
            key=extract_frame_num
        )

        if not frame_pdbs:
            raise FileNotFoundError(f"No PDB frames found inside {tar_path}")

        print(f"Loaded {len(frame_pdbs)} frames from {tar_path} for RMSD calculation against {ref_pdb_path}")

        for pdb_path in frame_pdbs:
            frame_num = extract_frame_num(pdb_path)
            if frame_num < 0:
                continue

            frame_obj = f"frame_{frame_num}"
            cmd.load(pdb_path, object=frame_obj)

            # 1. C-alpha RMSD (standard structural drift metric)
            ca_rms = cmd.align(f"{frame_obj} and name CA", f"{ref_obj} and name CA")[0]

            # 2. Backbone RMSD (C, CA, N, O)
            bb_rms = cmd.align(f"{frame_obj} and name C+CA+N+O", f"{ref_obj} and name C+CA+N+O")[0]

            # 3. All Heavy Atoms RMSD (non-hydrogen)
            heavy_rms = cmd.align(f"{frame_obj} and not elem H", f"{ref_obj} and not elem H")[0]

            data_rms.append({
                "frame": frame_num,
                "time_ns": float(frame_num),  # 1 frame per ns
                "rmsd_ca_A": round(ca_rms, 4),
                "rmsd_bb_A": round(bb_rms, 4),
                "rmsd_heavy_A": round(heavy_rms, 4)
            })

            cmd.delete(frame_obj)

    cmd.delete(ref_obj)

    # Build DataFrame, sort chronologically, and save CSV
    df = pd.DataFrame(data_rms)
    df = df.sort_values(by="time_ns").reset_index(drop=True)
    df.to_csv(output_csv, index=False)

    print(f"✅ RMSD calculations successfully saved to {output_csv} ({len(df)} entries)")


def main():
    if len(sys.argv) < 4:
        print("Usage: python calculate_rmsd.py <em_frames_tar> <canonical_pdb> <output_csv>")
        sys.exit(1)

    tar_path = sys.argv[1]
    ref_pdb_path = sys.argv[2]
    output_csv = sys.argv[3]

    calculate_rmsd_pymol(tar_path, ref_pdb_path, output_csv)


if __name__ == "__main__":
    main()