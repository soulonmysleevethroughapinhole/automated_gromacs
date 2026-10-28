Snakemake based workflow for 

Requires conda, GROMACS, MOPAC, and a pymol installation on your machine.

Can be run locally, or be used to orchestrate processes on HPC over ssh (use at yo own risk)

Connect to HPC over ssh, and leave the connection open while running the snakemake workflow. I think it'd be possible/preferable to mirror the SLURM queue, so that may be an avenue for improvement.



usage:

conda activate snakemake

if the code is modified, run snakemake --touch
note: if simulation is already running on HPC then touching won't complete

if sims are runnig use : snakemake --rerun-triggers mtime input --cores 8 

if all is fine:
snakemake --cores 8
