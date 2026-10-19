#!/usr/bin/env python3
import sys
import os
import re
import subprocess
import pandas as pd

def run_mopac(mop_path):
    """Executes MOPAC on the given .mop file."""
    # Run MOPAC via CLI (expects 'mopac' binary in PATH)
    cmd = f"mopac '{mop_path}'"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if res.returncode != 0:
        raise RuntimeError(f"MOPAC execution failed with return code {res.returncode}:\n{res.stderr}")

def parse_heat_of_formation(out_path):
    """
    Parses Heat of Formation (kcal/mol or kJ/mol) from MOPAC .out / .arc file.
    Looking for lines like: 'FINAL HEAT OF FORMATION =  -1234.56789 KCAL'
    """
    hof_kcal = None
    hof_kj = None

    if not os.path.exists(out_path):
        raise FileNotFoundError(f"Expected MOPAC output file not found: {out_path}")

    with open(out_path, "r", errors="replace") as f:
        for line in f:
            if "HEAT OF FORMATION" in line and "=" in line:
                # Match numerical value following the '=' sign
                match = re.search(r"=\s*([-\d\.]+)", line)
                if match:
                    val = float(match.group(1))
                    if "KCAL" in line.upper():
                        hof_kcal = val
                        hof_kj = val * 4.184
                    elif "KJ" in line.upper():
                        hof_kj = val
                        hof_kcal = val / 4.184

    if hof_kcal is None:
        raise ValueError(f"Could not extract 'HEAT OF FORMATION' from MOPAC output: {out_path}")

    return hof_kcal, hof_kj

def main(mop_input, csv_output):
    mop_abs = os.path.abspath(mop_input)
    csv_abs = os.path.abspath(csv_output)
    
    if not os.path.exists(mop_abs):
        raise FileNotFoundError(f"Input MOP file not found: {mop_abs}")

    os.makedirs(os.path.dirname(csv_abs), exist_ok=True)
    mop_dir = os.path.dirname(mop_abs)
    base_name = os.path.splitext(os.path.basename(mop_abs))[0]
    out_file = os.path.join(mop_dir, f"{base_name}.out")
    arc_file = os.path.join(mop_dir, f"{base_name}.arc")

    print(f"🚀 Running MOPAC calculation for: {mop_abs}")
    run_mopac(mop_abs)

    # Check .out first, fallback to .arc if .out is missing/truncated
    target_parse_file = out_file if os.path.exists(out_file) else arc_file
    
    print(f"🔍 Parsing Heat of Formation from: {target_parse_file}")
    hof_kcal, hof_kj = parse_heat_of_formation(target_parse_file)

    # Construct DataFrame
    df = pd.DataFrame([{
        "file": os.path.basename(mop_input),
        "heat_of_formation_kcal_mol": hof_kcal,
        "heat_of_formation_kj_mol": hof_kj
    }])

    df.to_csv(csv_abs, index=False)
    print(f"✅ Successfully wrote baseline Heat of Formation ({hof_kcal:.4f} kcal/mol) to {csv_abs}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_hof_for_mop.py <input_mop> <output_csv>")
        sys.exit(1)

    mop_file = sys.argv[1]
    csv_file = sys.argv[2]

    try:
        main(mop_file, csv_file)
    except Exception as e:
        print(f"❌ Error during MOPAC calculation: {str(e)}", file=sys.stderr)
        sys.exit(1)