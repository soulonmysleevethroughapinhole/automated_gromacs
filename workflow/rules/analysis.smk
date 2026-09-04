rule plot_run_hof_n_rmsd_vs_time:
    input:
        rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv",
        hof = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv"
    output:
        plot = "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_n_hof_vs_time.png"
    log:
        "logs/{pdb}/{source}/{model_id}/{protocol}/plot_rmsd_n_hof_vs_time.log"
    shell:
        "python scripts/smk/analysis/plot_run_hof_n_rmsd_vs_time.py {input.rmsd} {input.hof} {output.plot} {wildcards.pdb} {wildcards.source} {wildcards.model_id} > {log} 2>&1"


rule plot_rmsd_vs_hof:
    input:
        rmsd = "results/rmsd/{pdb}/{source}/{model_id}/{protocol}/pymol_results/frames/rmsd.csv",
        hof = "results/mopac/{pdb}/{source}/{model_id}/{protocol}/frames/heat_of_formation.csv"
    output:
        plot = "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_vs_hof.png"
    log:
        "logs/{pdb}/{source}/{model_id}/{protocol}/plot_rmsd_vs_hof.log"
    shell:
        "python scripts/smk/analysis/plot_rmsd_vs_hof.py {input.rmsd} {input.hof} {output.plot} {wildcards.pdb} {wildcards.source} {wildcards.model_id} > {log} 2>&1"


rule absolute_analysis_pdb:
    input:
        "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_n_hof_vs_time.png",
        "results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/rmsd_vs_hof.png"
    output:
        touch("results/analysis/plots/{pdb}/{source}/{model_id}/{protocol}/analysis_complete.txt")