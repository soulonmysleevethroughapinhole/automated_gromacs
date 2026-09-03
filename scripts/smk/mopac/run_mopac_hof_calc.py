#!/usr/bin/env python3
import sys
import os
import re
import glob
import logging
import tarfile
import tempfile
import subprocess
import pandas as pd
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

sys.path.insert(0, os.path.abspath("scripts/smk/sims"))
from hpc_engine import HPCJobManager

def parse_heat_of_formation(arc_filepath):
    """
    Extracts the Heat of Formation in KJ/MOL from a MOPAC .arc file using regex.
    """
    pattern = re.compile(r'=\s*([-\d.]+)\s*KJ/MOL', re.IGNORECASE)
    try:
        with open(arc_filepath, 'r', errors='ignore') as f:
            for line in f:
                if "HEAT OF FORMATION" in line and "KJ/MOL" in line:
                    match = pattern.search(line)
                    if match:
                        return float(match.group(1))
    except Exception as e:
        print(f"Error reading {arc_filepath}: {e}")
    return None

def extract_frame_num(filename):
    """Extracts integer frame index from strings like 'frame_100.arc' or 'frame100.mop'."""
    match = re.search(r'\d+', os.path.basename(filename))
    return int(match.group(0)) if match else -1

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_mopac_hof_calc.py <job_description> <output_csv>")
        sys.exit(1)

    job_description = os.path.abspath(sys.argv[1])
    output_csv = os.path.abspath(sys.argv[2])
    
    local_dir = os.path.dirname(output_csv)
    log_path = os.path.join(local_dir, "run_mopac_hof_calc.log")
    output_tar = os.path.join(local_dir, "mopac_results.tar.gz")

    # --- Robust Wildcard Extraction ---
    p = Path(job_description)
    parts = p.parts
    try:
        idx = parts.index("mopac")
        pdb = parts[idx + 1]
        source = parts[idx + 2]
        model_id = parts[idx + 3]
        protocol = parts[idx + 4]
    except (ValueError, IndexError):
        print(f"Error: Could not parse wildcards from job description path: {job_description}")
        sys.exit(1)

    target_id = f"{pdb}/{source}/{model_id}"
    prefix = f"{pdb}_{source}_{model_id}"

    # Locate input frame tar archive from gmx_em step
    input_em_tar = os.path.abspath(f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/FRAMES_compressed.tar.gz")
    if not os.path.exists(input_em_tar):
        # Fallback location
        input_em_tar = os.path.abspath(f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz")

    # Setup Logger
    logger = logging.getLogger(f"mopac_{prefix}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
    ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

    ssh_target = "komondor"
    hpc_user = os.environ.get("HPC_USER")
    hpc_base = os.environ.get("HPC_REMOTE_BASE")
    clean_base = hpc_base.lstrip("~/")
    remote_dir = f"{clean_base}/{pdb}/{source}/{model_id}/{protocol}/mopac"

    logger.info("Initializing HPC MOPAC HoF pipeline for %s...", target_id)

    # 1. Sync input EM archive and job script to remote directory on Komondor
    subprocess.run(["ssh", ssh_target, f"mkdir -p '{remote_dir}'"], check=True)
    rsync_cmd = f"rsync -e 'ssh' -avz '{input_em_tar}' '{job_description}' '{ssh_target}:{remote_dir}/'"
    subprocess.run(rsync_cmd, shell=True, check=True)

    # 2. Unpack EM frames on Komondor into EM_FRAMES
    unpack_cmd = f"ssh {ssh_target} \"cd '{remote_dir}' && mkdir -p EM_FRAMES && tar -xzf {os.path.basename(input_em_tar)} -C EM_FRAMES --strip-components=1 2>/dev/null || tar -xzf {os.path.basename(input_em_tar)}\""
    subprocess.run(unpack_cmd, shell=True, check=True)

    # 3. Instantiate Manager & Run Slurm Job
    job_script = os.path.basename(job_description)
    job_name = f"mopac_{prefix}"

    hpc = HPCJobManager(
        ssh_target=ssh_target,
        remote_dir=remote_dir,
        local_dir=local_dir,
        job_name=job_name,
        job_script=job_script,
        log_path=log_path,
        logger=logger,
        hpc_user=hpc_user
    )

    def prepare_and_check_completion():
        pack_cmd = f"ssh {ssh_target} \"cd '{remote_dir}' && [ -d MOPAC_STAGED ] && tar -czf mopac_results.tar.gz MOPAC_STAGED && echo YES || echo NO\""
        res = subprocess.run(pack_cmd, shell=True, capture_output=True, text=True, check=False)
        return "YES" in res.stdout

    if hpc.is_remote_finished("mopac_results.tar.gz") or prepare_and_check_completion():
        logger.info("🎉 MOPAC calculations already finished on HPC.")
    else:
        job_id = hpc.get_active_job_id()
        if not job_id:
            job_id = hpc.submit_job(unpack_job_archive=False)

        hpc.monitor_job(job_id, poll_interval_sec=120)

        if not prepare_and_check_completion():
            raise RuntimeError("MOPAC job completed on HPC but mopac_results.tar.gz was not created.")

    # 4. Pull results back locally
    logger.info("📦 Pulling MOPAC result archive locally...")
    hpc.pull_results(target_subdir="")

    # 5. Extract .arc files locally and compile heat_of_formation.csv
    if not os.path.exists(output_tar):
        raise FileNotFoundError(f"Expected retrieved archive missing: {output_tar}")

    logger.info("Parsing Heat of Formation data from MOPAC .arc outputs...")
    data_hof = []

    with tempfile.TemporaryDirectory() as temp_dir:
        with tarfile.open(output_tar, "r:gz") as tar:
            tar.extractall(path=temp_dir)

        arc_files = glob.glob(os.path.join(temp_dir, "**", "*.arc"), recursive=True)
        if not arc_files:
            logger.error("No .arc files found inside %s", output_tar)
            raise FileNotFoundError(f"No .arc files found in {output_tar}")

        for arc_path in arc_files:
            frame_num = extract_frame_num(arc_path)
            hof_val = parse_heat_of_formation(arc_path)
            
            # Convert frame index to time in ns (Frame 0 = 0.0 ns, Frame 100 = 100.0 ns)
            time_ns = float(frame_num)  # Assuming 1 frame per ns; scales dynamically if frame count varies
            data_hof.append({"frame": frame_num, "time_ns": time_ns, "hof": hof_val})

    # 6. Sort by frame index ascending and save CSV
    df = pd.DataFrame(data_hof)
    df = df.sort_values(by="frame").reset_index(drop=True)
    df.to_csv(output_csv, index=False)

    logger.info("✅ Heat of Formation CSV successfully written to %s (%d entries)", output_csv, len(df))

if __name__ == "__main__":
    main()