#!/bin/bash
set -e          # Exit on any command failure
set -o pipefail # Catch errors in piped commands (like the echo | gmx)

# --- Configuration ---
# Check if an argument was provided
if [ -z "$1" ]; then
    echo "ERROR: No PDB target specified. Usage: ./script.sh <target_name>"
    exit 1
fi

pdb="$1"
filename="${pdb}.pdb"

cleanup() {
    echo -e "\n\n[!] Interrupted by user. Exiting safely..."
    exit 130
}
trap cleanup SIGINT

# --- Helper Function for MDRun ---
# This checks if a checkpoint exists for the specific step to resume it
run_md() {
    local name=$1
    if [ -f "${name}.cpt" ]; then
        echo "--> Resuming ${name} from checkpoint..."
        gmx mdrun -v -deffnm "$name" -cpi -nb gpu -pme gpu
    else
        echo "--> Starting ${name}..."
        gmx mdrun -v -deffnm "$name" -nb gpu -pme gpu
    fi
}

if [ ! -f "repaired.pdb" ]; then
    pdbfixer $filename --output=repaired.pdb --add-atoms=heavy --add-residues
fi
filename="repaired.pdb"

# --- 1. System Preparation ---
if [ ! -f "${pdb}_solv.gro" ]; then
    echo "--> Processing PDB with pdb2gmx..."
    if ! gmx pdb2gmx -water tip3p -ff amber99sb-ildn -ignh -f $filename -o "${pdb}_processed.gro" -p topol.top; then
        echo "ERROR: pdb2gmx failed. This usually means your PDB is missing heavy atoms (like CG in GLU)."
        echo "Please repair the PDB using PDBFixer or similar tools and try again."
        exit 1
    fi
    gmx editconf -f "${pdb}_processed.gro" -o "${pdb}_newbox.gro" -d 1.0 -bt cubic
    #gmx solvate -cp "${pdb}_newbox.gro" -cs -o "${pdb}_solv.gro" -p topol.top
    gmx solvate -cp "${pdb}_newbox.gro" -cs spc216.gro -o "${pdb}_solv.gro" -p topol.top
fi

# --- 2. Ions ---
if [ ! -f "ion_b4em.gro" ]; then
    gmx grompp -v -f emw_steep_br.mdp -c "${pdb}_solv.gro" -o em_setup.tpr -p topol.top -maxwarn 1
    echo SOL | gmx genion -s em_setup.tpr -o ion_b4em.gro -p topol.top -pname NA -nname CL -neutral
fi

# --- 3. Minimization (Steepest Descent) ---
if [ ! -f "after_st.gro" ]; then
    gmx grompp -v -f emw_steep_br.mdp -c ion_b4em.gro -o st.tpr -p topol.top -maxwarn 1
    gmx mdrun -v -s st.tpr -o st.trr -c after_st.gro -g st.log
fi

# --- 4. Minimization (Conjugate Gradient) ---
if [ ! -f "after_cg.gro" ]; then
    gmx grompp -v -f emw_cg_br.mdp -c after_st.gro -o cg.tpr -p topol.top -maxwarn 1
    gmx mdrun -v -s cg.tpr -o cg.trr -c after_cg.gro -g cg.log
fi

# --- 5. Production MD ---
if [ ! -f "md.gro" ]; then
    # Only grompp if the TPR doesn't exist yet
    if [ ! -f "md.tpr" ]; then
        gmx grompp -f md_br.mdp -o md.tpr -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1
    fi
    # Use our helper function to handle resumption
    run_md "md"
fi

# --- 6. Post-Processing ---
if [ ! -d "FRAMES" ]; then
    echo "--> Correcting Periodic Boundary Conditions (PBC)..."
    
    # 1. Make molecules whole
    #echo "Protein" | gmx trjconv -f md.trr -s md.tpr -o md_whole.xtc -pbc mol -ur compact
    echo "Protein"  | gmx trjconv -f md.xtc -s md.tpr -o md_whole.xtc -pbc mol -ur compact
    # 2. Center and fit the protein to remove rotational/translational drift
    echo "Protein Protein Protein" | gmx trjconv -f md_whole.xtc -s md.tpr -o md_clean.xtc -center -fit rot+trans

    echo "--> Extracting cleaned frames..."
    mkdir -p FRAMES
    # 3. Extract individual frames from the cleaned trajectory
    echo "Protein" | gmx trjconv -f md_clean.xtc -o FRAMES/frame.pdb -s md.tpr -sep
    
    tar czf FRAMES_compressed.tar.gz FRAMES
    
    # Optional clean-up of massive intermediate trajectory files
    rm md_whole.xtc md_clean.xtc
fi

echo "Workflow complete!"
