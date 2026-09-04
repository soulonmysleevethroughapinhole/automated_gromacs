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
    Extracts Heat of Formation in KJ/MOL from MOPAC .arc files.
    Normalizes non-breaking spaces and raises None if missing/unparsed.
    """
    try:
        with open(arc_filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read().replace('\xa0', ' ')

        for line in content.splitlines():
            if "HEAT OF FORMATION" in line:
                # Search for the float immediately preceding KJ/MOL
                match_kj = re.search(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*KJ/MOL', line, re.IGNORECASE)
                if match_kj:
                    return float(match_kj.group(1))

                # Fallback search for KCAL/MOL (* 4.184)
                match_kcal = re.search(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*KCAL/MOL', line, re.IGNORECASE)
                if match_kcal:
                    return float(match_kcal.group(1)) * 4.184
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
    os.makedirs(local_dir, exist_ok=True)  # Ensure local directory exists before logging

    log_path = os.path.join(local_dir, "run_mopac_hof_calc.log")
    output_tar = os.path.join(local_dir, "mopac_results.tar.gz")

    # --- Dynamic Wildcard Path Extraction ---
    p = Path(job_description)
    parts = p.parts
    try:
        idx_mopac = parts.index("mopac")
        idx_std = parts.index("standard_100ns")
        pdb = parts[idx_mopac + 1]
        source = parts[idx_mopac + 2]
        model_id = "/".join(parts[idx_mopac + 3 : idx_std])
        protocol = parts[idx_std]
    except (ValueError, IndexError) as e:
        print(f"Error parsing path wildcards from {job_description}: {e}")
        sys.exit(1)

    target_id = f"{pdb}/{source}/{model_id}"
    prefix = f"{pdb}_{source}_{model_id.replace('/', '_')}"

    # Setup Logger
    logger = logging.getLogger(f"mopac_{prefix}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
    ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

    # Locate input frame tar archive from gmx_em step
    input_em_tar = os.path.abspath(f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz")
    if not os.path.exists(input_em_tar):
        input_em_tar = os.path.abspath(f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/FRAMES_compressed.tar.gz")

    if not os.path.exists(input_em_tar):
        logger.error("❌ Input EM frames archive not found: %s", input_em_tar)
        raise FileNotFoundError(f"Input EM frames archive not found: {input_em_tar}")

    ssh_target = "komondor"
    hpc_user = os.environ.get("HPC_USER", "c_drukv")
    hpc_base = os.environ.get("HPC_REMOTE_BASE", "c_drugint/valentin/remote_compute/automated_gromacs")
    clean_base = hpc_base.lstrip("~/")
    remote_dir = f"{clean_base}/{pdb}/{source}/{model_id}/{protocol}/mopac"

    logger.info("Initializing HPC MOPAC HoF pipeline for %s...", target_id)

    # 1. Prepare remote folder on Komondor and sync input EM frame archive + job script
    subprocess.run(["ssh", ssh_target, f"mkdir -p '{remote_dir}'"], check=True)
    rsync_cmd = f"rsync -e 'ssh' -avz '{input_em_tar}' '{job_description}' '{ssh_target}:{remote_dir}/'"
    subprocess.run(rsync_cmd, shell=True, check=True)

    # 2. Unpack EM frames on Komondor into EM_FRAMES directory
    unpack_cmd = f"ssh {ssh_target} \"cd '{remote_dir}' && mkdir -p EM_FRAMES && tar -xzf {os.path.basename(input_em_tar)} -C EM_FRAMES --strip-components=1 2>/dev/null || tar -xzf {os.path.basename(input_em_tar)}\""
    subprocess.run(unpack_cmd, shell=True, check=True)

    # 3. Instantiate HPC Manager & Execute Slurm Array Job
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
        # Requires MOPAC_STAGED to contain at least 1 .arc file before tarring
        pack_cmd = (
            f"ssh {ssh_target} \""
            f"cd '{remote_dir}' && "
            f"if [ -d MOPAC_STAGED ] && ls MOPAC_STAGED/*.arc >/dev/null 2>&1; then "
            f"  tar -czf mopac_results.tar.gz MOPAC_STAGED && echo YES; "
            f"else "
            f"  echo NO; "
            f"fi\""
        )
        res = subprocess.run(pack_cmd, shell=True, capture_output=True, text=True, check=False)
        return "YES" in res.stdout

    if prepare_and_check_completion():
        logger.info("🎉 MOPAC calculations already finished on HPC.")
    else:
        job_id = hpc.get_active_job_id()
        if not job_id:
            job_id = hpc.submit_job(unpack_job_archive=False)

        hpc.monitor_job(job_id, poll_interval_sec=120)

        if not prepare_and_check_completion():
            raise RuntimeError("MOPAC job completed on HPC but mopac_results.tar.gz was not generated or contains no .arc files.")

    # 4. Pull results archive back locally
    logger.info("📦 Pulling MOPAC result archive locally...")
    hpc.pull_results(target_subdir="")

    if not os.path.exists(output_tar):
        # Fallback move if pulled directly into parent directory
        alt_tar = os.path.join(os.path.dirname(local_dir), "mopac_results.tar.gz")
        if os.path.exists(alt_tar):
            os.rename(alt_tar, output_tar)
        else:
            raise FileNotFoundError(f"Expected retrieved archive missing: {output_tar}")

    # 5. Extract .arc files locally and compile heat_of_formation.csv
    logger.info("Parsing Heat of Formation data from MOPAC .arc outputs...")
    data_hof = []
    missing_hof_files = []

    with tempfile.TemporaryDirectory() as temp_dir:
        with tarfile.open(output_tar, "r:gz") as tar:
            tar.extractall(path=temp_dir, filter='data')

        arc_files = glob.glob(os.path.join(temp_dir, "**", "*.arc"), recursive=True)
        if not arc_files:
            logger.error("❌ No .arc files found inside %s", output_tar)
            raise FileNotFoundError(f"No .arc files found in {output_tar}")

        for arc_path in sorted(arc_files, key=extract_frame_num):
            fn = os.path.basename(arc_path)
            frame_num = extract_frame_num(arc_path)

            if frame_num < 0:
                logger.error("❌ Could not extract frame integer index from filename: %s", fn)
                missing_hof_files.append(fn)
                continue

            hof_val = parse_heat_of_formation(arc_path)

            if hof_val is None:
                logger.error("❌ Heat of Formation NOT FOUND or UNPARSABLE in %s", fn)
                missing_hof_files.append(fn)
            else:
                time_ns = float(frame_num)
                data_hof.append({"frame": frame_num, "time_ns": time_ns, "hof": hof_val})

    # STRICT CHECK: Fail immediately if ANY frame failed to calculate
    if missing_hof_files:
        err_msg = f"MOPAC Heat of Formation extraction failed for {len(missing_hof_files)} file(s): {', '.join(missing_hof_files)}"
        logger.error("💥 CRITICAL FAILURE: %s", err_msg)
        raise RuntimeError(err_msg)

    # 6. Sort by frame index ascending and save CSV
    df = pd.DataFrame(data_hof)
    df = df.sort_values(by="frame").reset_index(drop=True)
    df.to_csv(output_csv, index=False)

    logger.info("✅ Heat of Formation CSV successfully written to %s (%d frames validated)", output_csv, len(df))


if __name__ == "__main__":
    main()