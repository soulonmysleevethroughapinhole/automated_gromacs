# https://www.sciencedirect.com/science/article/pii/S0223523426004836
import sys
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pathlib import Path

def fix_missing_atoms(input_pdb, output_pdb, ph=7.4):
    """Uses PDBFixer to fill missing atoms and add hydrogens to the PDB file."""
    if not os.path.exists(input_pdb):
        raise FileNotFoundError(f"Input file not found: {input_pdb}")
    
    
    fixer = PDBFixer(filename=input_pdb)

    # mutations to original
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # missing 
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)
    try:
        Path(output_pdb).parent.mkdir(exist_ok=True, parents=True)

        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"Fixed missing atoms and saved cleaned file to: {output_pdb}")
    except Exception as e:
        raise IOError(f"Failed to write output file: {e}")
    
def main():
    if len(sys.argv) != 3:
        print("Usage: python prep_structure.py <input.pdb> <output.pdb>")
        sys.exit(1)
    try:
        fix_missing_atoms(sys.argv[1], sys.argv[2])
    except Exception as e:
        print(f"❌ Pipeline Error: {e}")
        # Return a non-zero exit code so Snakemake marks the rule as FAILED
        sys.exit(1)

if __name__ == '__main__':
    main()