Snakemake based workflow for 

Requires conda, GROMACS, MOPAC, and a pymol installation on your machine.

Can be run locally, or be used to orchestrate processes on HPC over ssh (use at yo own risk)



usage:

conda activate snakemake

if the code is modified, run snakemake --touch

snakemake --cores 8
