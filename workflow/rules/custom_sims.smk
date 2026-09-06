

import re

def parse_length_ns(protocol_str, cfg_dict):
    """
    Determines length in nanoseconds.
    Checks explicit `length_ns` in config first.
    If missing, parses strings like 'extended_1us' -> 1000.0 or 'extended_500ns' -> 500.0.
    """
    if "length_ns" in cfg_dict:
        return float(cfg_dict["length_ns"])
    
    # Parse length from protocol string (e.g., 'extended_1us' or 'extended_500ns')
    match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns)', protocol_str, re.IGNORECASE)
    if match:
        val, unit = float(match.group(1)), match.group(2).lower()
        return val * 1000.0 if unit == "us" else val
    
    raise ValueError(f"Could not parse simulation length from protocol: {protocol_str}")


def get_md_input_dependencies(wildcards):
    """
    Resolves input structure/checkpoint based on whether the run extends a previous run
    or starts freshly from energy minimization.
    """
    pdb = wildcards.pdb
    source = wildcards.source
    model_id = wildcards.model_id
    protocol = wildcards.protocol

    # Retrieve target entry from config
    target_entry = None
    for entry in cocustom_mdpnfig["custom_simulations"][pdb][source][model_id]:
        if entry["protocol"] == protocol:
            target_entry = entry
            break

    if target_entry is None:
        raise ValueError(f"Protocol '{protocol}' not found for {pdb}/{source}/{model_id}")

    # If extending a previous run, grab its checkpoint, topology, and TPR
    if "md_to_extend" in target_entry:
        parent_protocol = target_entry["md_to_extend"]
        parent_dir = f"results/gromacs/{pdb}/{source}/{model_id}/{parent_protocol}"
        return {
            "cpt": f"{parent_dir}/md.cpt",
            "tpr": f"{parent_dir}/md.tpr",
            "top": f"{parent_dir}/topol.top"
        }
    
    # Standard or Simulated Annealing run: starts freshly from GROMACS EM outputs
    #em_dir = f"results/gromacs/{pdb}/{source}/{model_id}/em_results"
    return {
		gro_cg = f"results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/after_cg.gro",
		top = f"results/gromacs/{pdb}/{source}/{model_id}/standard_100ns/topol.top",
        #"gro": f"{em_dir}/em_canonical_structure.gro",
        #"top": f"{em_dir}/topol.top"
    }


# for extended runs... ahh
# stage resolution map: N -> N-1
#STAGE_PRECEDENCE = {
#	'100ns': None,
#	'1us': '100ns',
#	'2us': '1us',
#	'5us': '2us',
#	'10us': '5us',
#}

