#!/usr/bin/env python3
import sys
import os
import logging
import subprocess
from pathlib import Path
from dotenv import load_dotenv
load_dotenv() # Load environment variables from .env file

sys.path.insert(0, os.path.abspath("scripts/smk/sims"))
from hpc_engine import HPCJobManager

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_HPC_gmx_em.py <job_description> <em_completed_sentinel>")
        sys.exit(1)

    job_description = os.path.abspath(sys.argv[1])
    em_completed_sentinel = os.path.abspath(sys.argv[2])
    
    local_dir = os.path.dirname(em_completed_sentinel)
    log_path = os.path.join(local_dir, "run_HPC_gmx_em.log")
    
    # --- Robust Wildcard Extraction ---
    p = Path(job_description)
    parts = p.parts
    
    try:
        # Locates 'gmx_em' or 'results' in path to isolate {pdb}/{source}/{model_id}/{protocol}
        idx = parts.index("gmx_em")
        pdb = parts[idx + 1]
        source = parts[idx + 2]
        model_id = parts[idx + 3]
        protocol = parts[idx + 4]
    except (ValueError, IndexError):
        # Fallback if path structure varies
        logger.error("Failed to parse wildcards from job description path: %s", job_description)
        sys.exit(1)

    target_id = f"{pdb}/{source}/{model_id}"
    prefix = f"{pdb}_{source}_{model_id}"
    
    # Locate input frame tar archive
    input_tar = os.path.abspath(f"results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz")
    
    # Setup Logger
    logger = logging.getLogger(f"em_{prefix}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    
    fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
    ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

    ssh_target = "komondor"
    hpc_user = os.getenv("HPC_USER")
    hpc_base = os.getenv("HPC_REMOTE_BASE")
    clean_base = hpc_base.lstrip("~/")
    
    # Clean remote directory path (No absolute local paths injected)
    remote_dir = f"{clean_base}/{pdb}/{source}/{model_id}/{protocol}/gmx_em"

    logger.info("Initializing HPC EM pipeline for %s...", target_id)
    logger.info("Remote Directory: %s", remote_dir)

    # 1. Sync input archive and job script to remote directory
    subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, f"mkdir -p '{remote_dir}'"], 
                   check=True,
                   env=os.environ)
    rsync_cmd = f"rsync -e 'ssh -o BatchMode=yes' -avz '{input_tar}' '{job_description}' '{ssh_target}:{remote_dir}/'"
    subprocess.run(rsync_cmd, shell=True, check=True, env=os.environ)

    # 2. Unpack frame archive on Komondor
    unpack_cmd = f"ssh -o BatchMode=yes {ssh_target} \"cd '{remote_dir}' && tar -xzf FRAMES_compressed.tar.gz\""
    subprocess.run(unpack_cmd, shell=True, check=True, env=os.environ)

    # 3. Instantiate HPC Manager
    job_script = os.path.basename(job_description)
    job_name = f"em_{prefix}"

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

    # 4. Custom post-execution check & tar compression on HPC
    def prepare_and_check_completion():
        pack_cmd = f"ssh -o BatchMode=yes {ssh_target} \"cd '{remote_dir}' && [ -d EM_FRAMES ] && tar -czf EM_FRAMES_compressed.tar.gz EM_FRAMES && echo YES || echo NO\""
        res = subprocess.run(pack_cmd, shell=True, capture_output=True, text=True, check=False, env=os.environ)
        return "YES" in res.stdout

    # 5. Execute pipeline
    if hpc.is_remote_finished("EM_FRAMES_compressed.tar.gz") or prepare_and_check_completion():
        logger.info("🎉 Energy minimization already finished on HPC.")
    else:
        job_id = hpc.get_active_job_id()
        if not job_id:
            job_id = hpc.submit_job(unpack_job_archive=False)
        
        hpc.monitor_job(job_id, poll_interval_sec=120)

        if not prepare_and_check_completion():
            raise RuntimeError("HPC Energy Minimization completed but EM_FRAMES_compressed.tar.gz was not produced.")

    # 6. Retrieve minimized archive locally
    hpc.pull_results(target_subdir="frames")

    # Move pulled archive into expected path
    local_pulled_tar = os.path.join(local_dir, "frames", "EM_FRAMES_compressed.tar.gz")
    expected_output_tar = os.path.join(local_dir, "frames", "FRAMES_compressed.tar.gz")
    if os.path.exists(local_pulled_tar):
        os.rename(local_pulled_tar, expected_output_tar)

    # 7. Write completion sentinel
    with open(em_completed_sentinel, "w") as f:
        f.write(f"Energy minimization completed for {target_id}.\n")

    logger.info("✅ Energy Minimization workflow complete for %s", target_id)

if __name__ == "__main__":
    main()