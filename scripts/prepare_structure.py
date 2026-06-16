import os
import glob

#another method of doing it?  TODO check for edge cases to see which works better
import pras

# Load, fix missing main-chain loops, fix side-chains, and save
structure = pras.StructuralBiology("'data/proteins/empirical/sars_cov_2_wuhan/7UU9/7UU9.pdb")
structure.repair_missing_residues()
structure.write_pdb("wuhan_fixed_pras.pdb")

#snakemake.input[0]

# Collect all raw CIF files for SARS-CoV-2 Wuhan strain proteins, 
# later may extend to new protein groups listed in config.yaml
#RAW_CIFS = glob.glob("data/proteins/empirical/sars_cov_2_wuhan/*/*.cif")
#TARGETS = [os.path.basename(os.path.dirname(cif)) for cif in RAW_CIFS]

#p#rint(RAW_CIFS)