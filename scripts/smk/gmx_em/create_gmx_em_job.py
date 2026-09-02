#!/usr/bin/env python3
import sys
import os
import textwrap
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()  # Load environment variables from .env file

HPC_HOST = os.getenv("HPC_HOST")
HPC_USER = os.getenv("HPC_USER")
HPC_ACCOUNT = os.getenv("HPC_ACCOUNT")
HPC_REMOTE_BASE = os.getenv("HPC_REMOTE_BASE")
HPC_PARTITION = os.getenv("HPC_PARTITION")
HPC_GMX_HOME = os.getenv("HPC_GMX_HOME")
HPC_TIME = os.getenv("HPC_TIME", "00:10:00")
SUBMIT_HPC = os.getenv("SUBMIT_HPC", "0")
SSH_DAEMON_PORT = os.getenv("SSH_DAEMON_PORT", "22")

def generate_em_slurm_job(job_output_path, prefix):
    slurm_script = textwrap.dedent(f"""\
#!/bin/bash
#SBATCH --job-name=em_{prefix}
#SBATCH --output=em_%a_%j.out
#SBATCH --error=em_%a_%j.err
#SBATCH --account={HPC_ACCOUNT}
#SBATCH --partition={HPC_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-100%32
#SBATCH --time={HPC_TIME}
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

FRAME_ID=${{SLURM_ARRAY_TASK_ID}}
INPUT_PDB="FRAMES/frame_${{FRAME_ID}}.pdb"
WORK_DIR="work_frame_${{FRAME_ID}}"

if [ ! -f "$INPUT_PDB" ]; then
    # Fallback check for ununderscored frame naming convention (frame0.pdb)
    INPUT_PDB="FRAMES/frame${{FRAME_ID}}.pdb"
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "Frame $INPUT_PDB not found, skipping."
    exit 0
fi

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ----------------------------------------------------------------------
# Write MDP Files
# ----------------------------------------------------------------------

# 1. Steepest Descent MDP Parameters
cat << 'EOF' > emw_steep_hydr.mdp
;
;	EM input altalanos
;
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
fourier_nx          =  0
fourier_ny          =  0
fourier_nz          =  0
pme_order           =  4
ewald_rtol          =  1e-5
optimize_fft        =  yes
;
;	Energy minimizing stuff
;
emtol               =  100
emstep              =  0.5
EOF

# 2. Conjugate Gradient MDP Parameters
cat << 'EOF' > emw_cg_hydr.mdp
;
;	EM input altalanos
;
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
fourier_nx          =  0
fourier_ny          =  0
fourier_nz          =  0
pme_order           =  4
ewald_rtol          =  1e-5
optimize_fft        =  yes
;
;	Energy minimizing stuff
;
emtol               =  10
emstep              =  0.05
EOF

# ----------------------------------------------------------------------
# Pipeline Execution
# ----------------------------------------------------------------------

# 1. Generate Topology & Hydrogen Placement
gmx_mpi pdb2gmx -f "../$INPUT_PDB" -o "frame_${{FRAME_ID}}_processed.gro" -p "topol.top" -ff amber99sb-ildn -water tip3p -ignh >> em_frame.log 2>&1

# 2. Define Box and Solvate
gmx_mpi editconf -f "frame_${{FRAME_ID}}_processed.gro" -o "frame_${{FRAME_ID}}_newbox.gro" -d 1.0 -bt cubic >> em_frame.log 2>&1
gmx_mpi solvate -cp "frame_${{FRAME_ID}}_newbox.gro" -cs spc216.gro -o "frame_${{FRAME_ID}}_solv.gro" -p topol.top >> em_frame.log 2>&1

# 3. Add Ions (Neutralize)
gmx_mpi grompp -v -f emw_steep_hydr.mdp -c "frame_${{FRAME_ID}}_solv.gro" -r "frame_${{FRAME_ID}}_solv.gro" -o ion_prep.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
echo "SOL" | gmx_mpi genion -s ion_prep.tpr -o ion_b4em.gro -p topol.top -pname NA -nname CL -neutral >> em_frame.log 2>&1

# 4. Stage 1: Steepest Descent Minimization
gmx_mpi grompp -v -f emw_steep_hydr.mdp -c ion_b4em.gro -r ion_b4em.gro -o st.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
gmx_mpi mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log -ntomp 4 >> em_frame.log 2>&1

# 5. Stage 2: Conjugate Gradient Minimization
gmx_mpi grompp -v -f emw_cg_hydr.mdp -c after_st.gro -r after_st.gro -o cg.tpr -p topol.top -maxwarn 2 >> em_frame.log 2>&1
gmx_mpi mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log -ntomp 4 >> em_frame.log 2>&1

# 6. Extract Final Cleaned Protein Structure
mkdir -p ../EM_FRAMES
echo "Protein" | gmx_mpi trjconv -s cg.tpr -f after_cg.gro -o "../EM_FRAMES/frame_${{FRAME_ID}}.pdb" -pbc mol -ur compact >> em_frame.log 2>&1

cd ..
rm -rf "$WORK_DIR"
    """)

    os.makedirs(os.path.dirname(job_output_path), exist_ok=True)
    with open(job_output_path, "w") as f:
        f.write(slurm_script)

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_gmx_em_job.py <input_tar> <output_job_script>")
        sys.exit(1)

    input_tar = sys.argv[1]
    output_job = sys.argv[2]
    
    # Robust wildcard parsing using pathlib
    p = Path(input_tar)
    parts = p.parts  # Converts path to normalized tuple of components
    
    # Search for 'gromacs' in path to locate relative indices safely
    try:
        gmx_idx = parts.index("gromacs")
        pdb = parts[gmx_idx + 1]
        source = parts[gmx_idx + 2]
        model_id = parts[gmx_idx + 3]
        prefix = f"{pdb}_{source}_{model_id}"
    except (ValueError, IndexError):
        # Fallback to output script filename parsing
        job_filename = Path(output_job).name
        prefix = job_filename.replace("_em_job.job", "").replace("_md_job.job", "")

    generate_em_slurm_job(output_job, prefix)
    print(f"Generated Slurm EM Job script for {prefix}: {output_job}")

if __name__ == "__main__":
    main()