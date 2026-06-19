# this workflow will be used for acquiring the AlphaFold and ColabFold generated protein structures
dimport json

checkpoint screen_structure:
    input:
        empirical: 'data/{pdb}/empirical/clean_structure.pdb',
        predicted: 'data/predictoins/{pdb}/{proctype}/{pdb}_{proctype}_{model_number}.pdb'
    output:
        report = 'results/'

def aggregate_md_inputs(wildcards):
    # 1. Force Snakemake to wait for the checkpoint to complete and evaluate
    checkpoint_output = checkpoints.screen_structure.get(**wildcards).output.report
    
    # 2. Read the evaluation data
    with open(checkpoint_output) as f:
        data = json.load(f)
    
    # 3. Decide the DAG path based on structural compliance
    if data["rmsd"] < config["max_allowed_rmsd"]:
        # Target passed! Route it to the production pipeline
        target_length = config["targets"][wildcards.pdb]["length_ns"]
        return f"results/2_simulations/{wildcards.pdb}/{target_length}ns/md.gro"
    else:
        # Target failed! Route it to an alternative terminal endpoint
        return f"results/1_predicted_models/rejected_{wildcards.pdb}.txt"

rule all:
    input:
        # Dynamically resolved target paths
        expand(aggregate_md_inputs, pdb=TARGETS)