#!/usr/bin/env python3
import sys
import os
import textwrap
from pathlib import Path
from dotenv import load_dotenv
import tarfile
import re

load_dotenv()

HPC_ACCOUNT = os.getenv("HPC_ACCOUNT")
HPC_PARTITION = os.getenv("HPC_PARTITION")
HPC_GMX_HOME = os.getenv("HPC_GMX_HOME")

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
	"""Inspects tarball to dynamically determine frame count without memory overhead."""
	try:
		with tarfile.open(tar_path, "r:*") as tar:
			count = sum(1 for member in tar if member.name.endswith(".pdb"))
			if count == 0:
				raise ValueError(f"No .pdb frames found inside archive '{tar_path}'.")
			return count
	except (tarfile.TarError, EOFError, OSError) as e:
		sys.stderr.write(f"[ERROR] Corrupted or incomplete tarball '{tar_path}': {e}\n")
		raise

def generate_em_slurm_job(job_output_path, prefix, total_frames):
	MAX_ALLOWED_TASKS = 500
	target_frames = min(total_frames, MAX_ALLOWED_TASKS)
	max_array_index = max(0, target_frames - 1)

	raw_job_name = f"em_{prefix}"
	job_name = raw_job_name[:60] if len(raw_job_name) > 60 else raw_job_name

	slurm_script = textwrap.dedent(f"""\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=em_%a_%j.out
#SBATCH --error=em_%a_%j.err
#SBATCH --account={HPC_ACCOUNT}
#SBATCH --partition={HPC_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-{max_array_index}%32
#SBATCH --time="00:15:00"
#SBATCH --no-requeue

# Load GROMACS modules
module --force purge
module load PrgEnv-gnu/8.6.0
module load cray-mpich
module load cray-fftw
module load gcc/12.2.0

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export gmxhome={HPC_GMX_HOME}
export PATH="${{gmxhome}}/bin:${{PATH}}"
export LD_LIBRARY_PATH="${{gmxhome}}/lib64:${{gmxhome}}/lib:${{LD_LIBRARY_PATH}}"

# Atomic Directory Unpack Lock
if [ -f "FRAMES_compressed.tar.gz" ] && [ ! -d "FRAMES" ]; then
	if mkdir FRAMES.lock 2>/dev/null; then
		if [ ! -d "FRAMES" ]; then
			mkdir -p FRAMES
			tar -xzf FRAMES_compressed.tar.gz -C .
		fi
		rm -rf FRAMES.lock
	else
		while [ ! -d "FRAMES" ] || [ -d FRAMES.lock ]; do
			sleep 1
		done
	fi
fi

shopt -s extglob nullglob
FRAME_FILES=($(ls -v FRAMES/*.pdb 2>/dev/null))
TOTAL_FRAMES=${{#FRAME_FILES[@]}}

if [ $TOTAL_FRAMES -eq 0 ]; then
	echo "❌ ERROR: No PDB frame files found in FRAMES/"
	exit 1
fi

# Compute strided frame index if total frames > 500
TARGET_TASKS={target_frames}
FRAME_IDX=$(( SLURM_ARRAY_TASK_ID * TOTAL_FRAMES / TARGET_TASKS ))
INPUT_PDB="${{FRAME_FILES[${{FRAME_IDX}}]}}"

if [ -z "$INPUT_PDB" ] || [ ! -f "$INPUT_PDB" ]; then
	echo "Task ID ${{SLURM_ARRAY_TASK_ID}} mapped to index ${{FRAME_IDX}} out of range or file missing, skipping."
	exit 0
fi

# Extract base filename for job isolation
FRAME_NAME=$(basename "$INPUT_PDB" .pdb)
WORK_DIR="work_${{FRAME_NAME}}"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"
set -e

cat << 'EOF' > emw_steep_hydr.mdp
cpp                 =  /lib/cpp
define              =  -DFLEXIBLE -DPOSRES
constraints         =  none
integrator          =  steep
dt                  =  0.002
nsteps              =  250000
nstlist             =  10
ns_type             =  grid
rlist               =  0.9
coulombtype         =  PME
rcoulomb            =  0.9
rvdw                =  0.9
fourierspacing      =  0.12
pme_order           =  4
ewald_rtol          =  1e-5
optimize_fft        =  yes
emtol               =  100
emstep              =  0.5
EOF

cat << 'EOF' > emw_cg_hydr.mdp
cpp                 =  /lib/cpp
define              =  -DFLEXIBLE -DPOSRES
constraints         =  none
integrator          =  cg
dt                  =  0.001
nsteps              =  50000
nstlist             =  10
ns_type             =  grid
rlist               =  0.9
coulombtype         =  PME
rcoulomb            =  0.9
rvdw                =  0.9
fourierspacing      =  0.12
pme_order           =  4
ewald_rtol          =  1e-5
optimize_fft        =  yes
emtol               =  10
emstep              =  0.05
EOF

gmx_mpi pdb2gmx -f "../$INPUT_PDB" -o "${{FRAME_NAME}}_processed.gro" -p "topol.top" -ff amber99sb-ildn -water tip3p -ignh >> em_frame.log 2>&1
gmx_mpi editconf -f "${{FRAME_NAME}}_processed.gro" -o "${{FRAME_NAME}}_newbox.gro" -d 1.0 -bt cubic >> em_frame.log 2>&1
gmx_mpi solvate -cp "${{FRAME_NAME}}_newbox.gro" -cs spc216.gro -o "${{FRAME_NAME}}_solv.gro" -p topol.top >> em_frame.log 2>&1

gmx_mpi grompp -v -f emw_steep_hydr.mdp -c "${{FRAME_NAME}}_solv.gro" -r "${{FRAME_NAME}}_solv.gro" -o ion_prep.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
echo "SOL" | gmx_mpi genion -s ion_prep.tpr -o ion_b4em.gro -p topol.top -pname NA -nname CL -neutral >> em_frame.log 2>&1

gmx_mpi grompp -v -f emw_steep_hydr.mdp -c ion_b4em.gro -r ion_b4em.gro -o st.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
gmx_mpi mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log -ntomp 4 >> em_frame.log 2>&1

gmx_mpi grompp -v -f emw_cg_hydr.mdp -c after_st.gro -r after_st.gro -o cg.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
gmx_mpi mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log -ntomp 4 >> em_frame.log 2>&1

mkdir -p ../EM_FRAMES
echo "Protein" | gmx_mpi trjconv -s cg.tpr -f after_cg.gro -o "../EM_FRAMES/${{FRAME_NAME}}.pdb" -pbc mol -ur compact >> em_frame.log 2>&1

cd ..
rm -rf "$WORK_DIR"

COMPLETED_COUNT=$(ls -1 EM_FRAMES/*.pdb 2>/dev/null | wc -l)
if [ "$COMPLETED_COUNT" -ge "{target_frames}" ]; then
	touch EM_COMPLETE.txt
fi
	""")

	os.makedirs(os.path.dirname(job_output_path), exist_ok=True)
	with open(job_output_path, "w") as f:
		f.write(slurm_script)


def main():
	if len(sys.argv) < 3:
		print("Usage: python create_gmx_em_job.py <input_tar> <output_job_script>")
		sys.exit(1)

	input_tar = os.path.abspath(sys.argv[1])
	output_job = os.path.abspath(sys.argv[2])

	p = Path(input_tar)
	parts = p.parts

	try:
		gmx_idx = parts.index("gromacs")
		if "md_results" in parts:
			protocol_idx = parts.index("md_results") - 1
		elif "frames" in parts:
			protocol_idx = parts.index("frames") - 2
		else:
			protocol_idx = -3

		protocol = parts[protocol_idx]
		pdb = parts[gmx_idx + 1]
		source = parts[gmx_idx + 2]
		model_id = "_".join(parts[gmx_idx + 3 : protocol_idx])
		prefix = f"{pdb}_{source}_{model_id}_{protocol}"
	except (ValueError, IndexError):
		job_filename = Path(output_job).name
		prefix = job_filename.replace("_em_job.job", "").replace("_md_job.job", "")

	total_frames = count_frames_in_tar(input_tar)
	generate_em_slurm_job(output_job, prefix, total_frames)
	print(f"Generated Slurm EM Job script for {prefix} ({total_frames} frames): {output_job}")


if __name__ == "__main__":
	main()