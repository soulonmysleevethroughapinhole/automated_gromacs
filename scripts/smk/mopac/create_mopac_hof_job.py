#!/usr/bin/env python3
import sys
import os
import textwrap
import tarfile
import glob
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

HPC_ACCOUNT = os.getenv("HPC_ACCOUNT")
HPC_PARTITION = os.getenv("HPC_PARTITION", "cpu")
HPC_TIME = os.getenv("HPC_TIME", "02:00:00")
HPC_MOPAC_HOME = os.getenv("HPC_MOPAC_HOME")

def count_frames_in_tar(tar_path):
    """Inspects tarball or directory locally to determine total frame count."""
    if os.path.isfile(tar_path) and tarfile.is_tarfile(tar_path):
        with tarfile.open(tar_path, "r:*") as tar:
            members = [m.name for m in tar.getmembers() if m.name.endswith(".pdb")]
            return len(members)
    elif os.path.isdir(tar_path):
        pdbs = glob.glob(os.path.join(tar_path, "**", "*.pdb"), recursive=True)
        return len(pdbs)
    raise ValueError(f"Cannot count number of frames in tar or path: {tar_path}")

def generate_mopac_slurm_job(job_output_path, prefix, total_frames):
    max_array_index = max(0, total_frames - 1)

    slurm_script = textwrap.dedent(f"""\
#!/bin/bash
#SBATCH --job-name=mopac_{prefix}
#SBATCH --output=mopac_%a_%j.out
#SBATCH --error=mopac_%a_%j.err
#SBATCH --account={HPC_ACCOUNT}
#SBATCH --partition={HPC_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-{max_array_index}%32
#SBATCH --time="00:15:00"
#SBATCH --no-requeue

# Load MOPAC module
module --force purge
module load gcc/12.2.0 intel/tbb intel/compiler-rt intel/mkl

export mopachome={HPC_MOPAC_HOME}
export PATH="${{mopachome}}/bin:${{PATH}}"
export LD_LIBRARY_PATH="${{mopachome}}/lib64:${{mopachome}}/lib:${{LD_LIBRARY_PATH}}"

# Ensure EM_FRAMES directory exists on HPC before file globbing
if [ ! -d "EM_FRAMES" ]; then
    mkdir -p EM_FRAMES
    TAR_FILE=$(ls *.tar.gz 2>/dev/null | head -n 1)
    if [ -n "$TAR_FILE" ]; then
        tar -xzf "$TAR_FILE" -C EM_FRAMES --strip-components=1 2>/dev/null || tar -xzf "$TAR_FILE"
    fi
fi

# Dynamically map SLURM_ARRAY_TASK_ID to sorted PDB frame files
shopt -s extglob nullglob
FRAME_FILES=($(ls -v EM_FRAMES/*.pdb 2>/dev/null))

if [ ${{#FRAME_FILES[@]}} -eq 0 ]; then
    echo "❌ ERROR: No PDB frame files found in EM_FRAMES/"
    exit 1
fi

INPUT_PDB="${{FRAME_FILES[${{SLURM_ARRAY_TASK_ID}}]}}"

if [ -z "$INPUT_PDB" ] || [ ! -f "$INPUT_PDB" ]; then
    echo "Task ID ${{SLURM_ARRAY_TASK_ID}} out of range or missing file, skipping."
    exit 0
fi

FRAME_NAME=$(basename "$INPUT_PDB" .pdb)
mkdir -p MOPAC_STAGED
MOP_FILE="MOPAC_STAGED/${{FRAME_NAME}}.mop"

# 1. Convert PDB structure into MOPAC input file (.mop)
cat << 'EOF' > "$MOP_FILE"
PM7 1SCF EPS=78.3 MOZYME

EOF

# Append atomic coordinates starting from line 5
grep -E "^ATOM|^HETATM" "$INPUT_PDB" >> "$MOP_FILE"

# 2. Execute MOPAC calculation
cd MOPAC_STAGED
mopac "${{FRAME_NAME}}.mop" > /dev/null 2>&1
cd ..

echo "MOPAC calculation completed for frame ${{FRAME_NAME}}"
    """)

    os.makedirs(os.path.dirname(job_output_path), exist_ok=True)
    with open(job_output_path, "w") as f:
        f.write(slurm_script)

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_mopac_hof_job.py <input_tar> <output_job_description>")
        sys.exit(1)

    input_tar = sys.argv[1]
    output_job = sys.argv[2]

    # Safely extract job prefix from output_job path: results/mopac/{pdb}/{source}/{model_id}/...
    try:
        p = Path(output_job)
        parts = p.parts
        if "mopac" in parts:
            idx = parts.index("mopac")
            # Joins remaining directory elements up to 'frames'
            pdb = parts[idx + 1]
            source = parts[idx + 2]
            model_id = parts[idx + 3]
            protocol = parts[idx + 4]
            prefix = f"{pdb}_{source}_{model_id}_{protocol}"
        else:
            prefix = p.parent.parent.name
    except Exception:
        prefix = "default_mopac"

    total_frames = count_frames_in_tar(input_tar)
    generate_mopac_slurm_job(output_job, prefix, total_frames)
    print(f"Successfully generated MOPAC Slurm Job script for {prefix} ({total_frames} frames): {output_job}")

if __name__ == "__main__":
    main()