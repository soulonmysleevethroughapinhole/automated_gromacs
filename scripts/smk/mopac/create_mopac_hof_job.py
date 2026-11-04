#!/usr/bin/env python3
import sys
import os
import textwrap
import tarfile
from pathlib import Path
from dotenv import load_dotenv
import re

load_dotenv()

HPC_ACCOUNT = os.getenv("HPC_ACCOUNT")
HPC_PARTITION = os.getenv("HPC_PARTITION", "cpu")
HPC_TIME = os.getenv("HPC_TIME", "02:15:00")
HPC_MOPAC_HOME = os.getenv("HPC_MOPAC_HOME")

def parse_length_ns(protocol_str: str) -> float:
	"""Extracts simulation duration in nanoseconds from protocol string."""
	match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns)', protocol_str, re.IGNORECASE)
	if match:
		val, unit = float(match.group(1)), match.group(2).lower()
		return val * 1000.0 if unit == "us" else val
	raise ValueError(f"Could not parse simulation time length from protocol string '{protocol_str}'")


def get_protocol_frame_limit(protocol_str: str) -> int:
	"""
	Returns protocol-dependent max frame target:
	- standard_100ns (<= 100ns) -> 100 frames
	- long runs (> 100ns, e.g. 1us) -> 500 frames
	"""
	length_ns = parse_length_ns(protocol_str)
	return 100 if length_ns <= 100.0 else 500


def validate_input_frames(tar_path: str, protocol_str: str, total_frames: int) -> int:
	"""
	Validates frame count against protocol expectations and returns the target task count.
	"""
	expected_limit = get_protocol_frame_limit(protocol_str)

	if total_frames < expected_limit:
		raise ValueError(
			f"❌ INCOMPLETE TRAJECTORY: '{tar_path}' contains {total_frames} frames, "
			f"but protocol '{protocol_str}' requires {expected_limit} frames. "
			f"Simulation/EM likely crashed or ended prematurely."
		)

	return expected_limit

def count_frames_in_tar(tar_path):
	"""Inspects tarball to determine total frame count safely without memory spikes."""
	try:
		with tarfile.open(tar_path, "r:*") as tar:
			count = sum(1 for member in tar if member.name.endswith(".pdb"))
			if count == 0:
				raise ValueError(f"No PDB frames found inside archive '{tar_path}'.")
			return count
	except (tarfile.TarError, EOFError, OSError) as e:
		sys.stderr.write(f"[ERROR] Corrupted or incomplete tarball '{tar_path}': {e}\n")
		raise

def generate_mopac_slurm_job(job_output_path, prefix, total_frames, target_frames):
	max_array_index = max(0, target_frames - 1)

	raw_job_name = f"mopac_{prefix}"
	job_name = raw_job_name[:60] if len(raw_job_name) > 60 else raw_job_name

	slurm_script = textwrap.dedent(f"""\
#!/bin/bash
#SBATCH --job-name={job_name}
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

# Load MOPAC environment
module --force purge
module load gcc/12.2.0 intel/tbb intel/compiler-rt intel/mkl 2>/dev/null || true

export mopachome={HPC_MOPAC_HOME}
export PATH="${{mopachome}}/bin:${{PATH}}"
export LD_LIBRARY_PATH="${{mopachome}}/lib64:${{mopachome}}/lib:${{LD_LIBRARY_PATH}}"

# Atomic Directory Unpack Lock
if [ ! -d "EM_FRAMES" ]; then
	if mkdir EM_FRAMES.lock 2>/dev/null; then
		if [ ! -d "EM_FRAMES" ]; then
			mkdir -p EM_FRAMES
			TAR_FILE=$(ls *.tar.gz 2>/dev/null | head -n 1)
			if [ -n "$TAR_FILE" ]; then
				tar -xzf "$TAR_FILE" -C EM_FRAMES --strip-components=1 2>/dev/null || tar -xzf "$TAR_FILE" -C EM_FRAMES
			fi
		fi
		rm -rf EM_FRAMES.lock
	else
		while [ ! -d "EM_FRAMES" ] || [ -d EM_FRAMES.lock ]; do
			sleep 1
		done
	fi
fi

shopt -s extglob nullglob
FRAME_FILES=($(ls -v EM_FRAMES/*.pdb 2>/dev/null))
TOTAL_FRAMES=${{#FRAME_FILES[@]}}

if [ $TOTAL_FRAMES -eq 0 ]; then
	echo "❌ ERROR: No PDB frame files found in EM_FRAMES/"
	exit 1
fi

TARGET_TASKS={target_frames}

if [ "$TARGET_TASKS" -eq "$TOTAL_FRAMES" ]; then
	FRAME_IDX=$SLURM_ARRAY_TASK_ID
else
	FRAME_IDX=$(( SLURM_ARRAY_TASK_ID * (TOTAL_FRAMES - 1) / (TARGET_TASKS - 1) ))
fi

INPUT_PDB="${{FRAME_FILES[${{FRAME_IDX}}]}}"

if [ -z "$INPUT_PDB" ] || [ ! -f "$INPUT_PDB" ]; then
	echo "Task ID ${{SLURM_ARRAY_TASK_ID}} mapped to index ${{FRAME_IDX}} out of range or missing file, skipping."
	exit 0
fi

FRAME_NAME=$(basename "$INPUT_PDB" .pdb)
mkdir -p MOPAC_STAGED
MOP_FILE="MOPAC_STAGED/${{FRAME_NAME}}.mop"

cat << EOF > "$MOP_FILE"
PM7 1SCF EPS=78.3 MOZYME GEO-OK
MOPAC Single-Point Heat of Formation for ${{FRAME_NAME}}

EOF

grep -E "^ATOM|^HETATM" "$INPUT_PDB" >> "$MOP_FILE"

cd MOPAC_STAGED
mopac "${{FRAME_NAME}}.mop" > "${{FRAME_NAME}}_mopac_exec.log" 2>&1
cd ..

if [ -f "MOPAC_STAGED/${{FRAME_NAME}}.out" ]; then
	echo "✅ MOPAC calculation completed for frame ${{FRAME_NAME}}"
else
	echo "❌ ERROR: MOPAC failed to generate output for frame ${{FRAME_NAME}}"
fi

COMPLETED_COUNT=$(ls -1 MOPAC_STAGED/*.out 2>/dev/null | wc -l)
if [ "$COMPLETED_COUNT" -ge "{target_frames}" ]; then
	touch MOPAC_COMPLETE.txt
	echo "🎉 All {target_frames} MOPAC calculations complete. Touched MOPAC_COMPLETE.txt"
fi
	""")

	os.makedirs(os.path.dirname(job_output_path), exist_ok=True)
	with open(job_output_path, "w") as f:
		f.write(slurm_script)

def main():
	if len(sys.argv) < 3:
		print("Usage: python create_mopac_hof_job.py <input_tar> <output_job_description>")
		sys.exit(1)

	input_tar = os.path.abspath(sys.argv[1])
	output_job = os.path.abspath(sys.argv[2])

	p = Path(output_job)
	parts = p.parts

	try:
		if "mopac" in parts:
			idx = parts.index("mopac")
		elif "gromacs" in parts:
			idx = parts.index("gromacs")
		else:
			idx = 0

		protocol_idx = parts.index("JOB") - 1 if "JOB" in parts else -2

		pdb = parts[idx + 1]
		source = parts[idx + 2]
		model_id = "_".join(parts[idx + 3 : protocol_idx])
		protocol = parts[protocol_idx]
		prefix = f"{pdb}_{source}_{model_id}_{protocol}"
	except (ValueError, IndexError):
		prefix = Path(output_job).stem.replace("_mopac_job", "").replace("_job", "")

	total_frames = count_frames_in_tar(input_tar)
	generate_mopac_slurm_job(output_job, prefix, total_frames, total_frames)
	print(f"Generated Slurm MOPAC Job script for {prefix} ({total_frames} frames): {output_job}")


if __name__ == "__main__":
	main()