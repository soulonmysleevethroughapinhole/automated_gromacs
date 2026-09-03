#!/usr/bin/env python3
import sys
import os
import textwrap
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

HPC_ACCOUNT = os.getenv("HPC_ACCOUNT", "c_drugint")
HPC_PARTITION = os.getenv("HPC_PARTITION", "cpu")
HPC_TIME = os.getenv("HPC_TIME", "02:00:00")
HPC_MOPAC_HOME = os.getenv("HPC_MOPAC_HOME")

def generate_mopac_slurm_job(job_output_path, prefix):
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
#SBATCH --array=0-100%32
#SBATCH --time="00:15:00"
#SBATCH --no-requeue

# Load MOPAC module
module --force purge
#module load mopac/22.1.1-gcc-12.2.0 || module load mopac 2>/dev/null || true
module load gcc/12.2.0 intel/tbb intel/compiler-rt intel/mkl

export mopachome={HPC_MOPAC_HOME}
export PATH="${{mopachome}}/bin:${{PATH}}"
export LD_LIBRARY_PATH="${{mopachome}}/lib64:${{mopachome}}/lib:${{LD_LIBRARY_PATH}}"

FRAME_ID=${{SLURM_ARRAY_TASK_ID}}
INPUT_PDB="EM_FRAMES/frame_${{FRAME_ID}}.pdb"

if [ ! -f "$INPUT_PDB" ]; then
    INPUT_PDB="EM_FRAMES/frame${{FRAME_ID}}.pdb"
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "PDB Frame $INPUT_PDB not found, skipping."
    exit 0
fi

mkdir -p MOPAC_STAGED
MOP_FILE="MOPAC_STAGED/frame_${{FRAME_ID}}.mop"

# 1. Convert PDB structure into MOPAC input file (.mop)
cat << 'EOF' > "$MOP_FILE"
PM7 1SCF EPS=78.3 MOZYME

EOF

# Append atomic coordinates starting from line 5
grep -E "^ATOM|^HETATM" "$INPUT_PDB" >> "$MOP_FILE"

# 2. Execute MOPAC calculation
cd MOPAC_STAGED
mopac "frame_${{FRAME_ID}}.mop" > /dev/null 2>&1
cd ..

echo "MOPAC calculation completed for frame ${{FRAME_ID}}"
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
            prefix = f"{pdb}_{source}_{model_id}"
        else:
            prefix = p.parent.parent.name
    except Exception:
        prefix = "default_mopac"

    generate_mopac_slurm_job(output_job, prefix)
    print(f"Successfully generated MOPAC Slurm Job script: {output_job}")

if __name__ == "__main__":
    main()