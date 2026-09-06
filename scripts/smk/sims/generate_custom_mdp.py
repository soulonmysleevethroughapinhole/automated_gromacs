#!/usr/bin/env python3
import sys
import yaml
import re

BASE_MDP_TEMPLATE = """\
; Base Production MDP
cpp                 =  /usr/bin/cpp
constraints         =  all-bonds
continuation            = {continuation} ; yes for extensions, no for new runs
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {nsteps}    
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 500000      ; save energies every 1000 ps
nstlog                  = 500000      ; update log file every 1000 ps
nstxout-compressed      = 500000      ; save compressed coordinates every 1000 ps
compressed-x-grps       = System    ; save the whole system

nstlist             =  10
ns_type             =  grid
coulombtype         =  PME
cutoff-scheme       =  Verlet
fourierspacing      =  0.12
rlist               =  1.1
rcoulomb            =  1.1
rvdw                =  1.1

Tcoupl              =  v-rescale
tc-grps		        =  Protein	non-Protein
tau_t               =  0.1	0.1
ref_t               =  300	300

Pcoupl              =  Parrinello-Rahman
tau_p                   = 2.0       ; 2.0 ps recommended for Parrinello-Rahman stability
;tau_p               =  0.5
compressibility     =  4.5e-5
ref_p               =  1.0
refcoord_scaling    =  all

gen_vel             =  {generate_velocities}
gen_seed            =  28480426

;gen_temp            =  298

"""

def parse_length_ns(full_protocol, cfg_dict):
    """
    Parses simulation length in nanoseconds.
    1. Checks explicit `length_ns` key in config dictionary.
    2. Parses the duration suffix from full_protocol (e.g., 'extended_1000ns' -> 1000.0, 'extended_1us' -> 1000.0).
    """
    if "length_ns" in cfg_dict:
        return float(cfg_dict["length_ns"])

    # Match numbers followed by 'us' or 'ns' anywhere in the protocol name
    match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns)', full_protocol, re.IGNORECASE)
    if match:
        val, unit = float(match.group(1)), match.group(2).lower()
        return val * 1000.0 if unit == "us" else val

    raise ValueError(f"Could not parse length from protocol string '{full_protocol}' or config: {cfg_dict}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python generate_mdp.py <full_protocol> <cfg_yaml_path> <out_mdp_path>")
        sys.exit(1)

    full_protocol = sys.argv[1]
    cfg_yaml_path = sys.argv[2]
    out_mdp_path = sys.argv[3]

    base_protocol = full_protocol.split("_")[0]  # e.g., "extended" or "sa"

    with open(cfg_yaml_path, 'r') as f:
        cfg = yaml.safe_load(f)

    length_ns = parse_length_ns(full_protocol, cfg)
    nsteps = int((length_ns * 1000.0) / 0.002) # 2 fs timestep

    # Determine if this run extends an earlier stage
    is_extension = (base_protocol == "extended") or ("md_to_extend" in cfg)

    mdp_text = BASE_MDP_TEMPLATE.format(
        nsteps=nsteps,
        continuation="yes" if is_extension else "no",
        generate_velocities="no" if is_extension else "yes"
    )

    # Append simulated annealing schedule if present
    if "temperature_schedule" in cfg:
        sched = cfg["temperature_schedule"]
        npts = len(sched)
        times = [str(int(pt["time_ns"] * 1000)) for pt in sched]
        temps = [str(pt["temp_k"]) for pt in sched]

        times_str = " ".join(times) + " " + " ".join(times)
        temps_str = " ".join(temps) + " " + " ".join(temps)

        mdp_text += f"""
; Simulated Annealing
annealing               = single single
annealing_npoints       = {npts} {npts}
annealing_time          = {times_str}
annealing_temp          = {temps_str}
"""

    with open(out_mdp_path, 'w') as f:
        f.write(mdp_text)

if __name__ == "__main__":
    main()