# workflow/utils/md_engine.py

import os
import time
import logging
import subprocess

def execute_hpc_md(
    wildcards,
    job_description_path,
    output_done_path,
    log_path,
    hpc_user,
    hpc_host,
    hpc_remote_base,
    get_progress_fn=None
):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logger = logging.getLogger(f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_{wildcards.protocol}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    target_id = f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
    prefix = f"{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}_md"
    job_dir = os.path.dirname(job_description_path)
    local_dir = os.path.abspath(os.path.dirname(output_done_path))

    ssh_target = f"{hpc_user}@{hpc_host}"
    clean_base = hpc_remote_base.lstrip("~/")
    remote_dir = os.path.join(
        clean_base,
        f"{wildcards.pdb}/{wildcards.source}/{wildcards.model_id}/{wildcards.protocol}"
    )
    job_name = f"md_{wildcards.pdb}_{wildcards.source}_{wildcards.model_id}"
    job_script = os.path.basename(job_description_path)

    try:
        # Guarantee remote directory exists
        subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, f"mkdir -p '{remote_dir}'"], check=False)

        # 1. Completion Check
        check_cmd = ["ssh", "-o", "BatchMode=yes", ssh_target, f"if [ -f '{remote_dir}/{prefix}.gro' ]; then echo YES; else echo NO; fi"]
        try:
            remote_done = subprocess.check_output(check_cmd, text=True).strip()
        except Exception:
            remote_done = "NO"

        if remote_done != "YES":
            # 2. Check Queue & Slurm History
            squeue_cmd = ["ssh", "-o", "BatchMode=yes", ssh_target, f"squeue -u {hpc_user} -n {job_name} -h -o '%i %t'"]
            q_res = subprocess.run(squeue_cmd, capture_output=True, text=True, check=False)
            q_out = q_res.stdout.strip()

            job_id = None
            if q_out:
                job_id = q_out.split()[0]
            else:
                sacct_cmd = ["ssh", "-o", "BatchMode=yes", ssh_target, f"sacct -u {hpc_user} -n {job_name} --format=State -n | head -n 1"]
                sacct_res = subprocess.run(sacct_cmd, capture_output=True, text=True, check=False)
                if "COMPLETED" not in sacct_res.stdout:
                    # Submit Slurm Job
                    submit_script = f"cd '{remote_dir}' && if [ -f JOB.tar.gz ]; then tar -xzf JOB.tar.gz; fi && sbatch {job_script}"
                    sub_res = subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, submit_script], capture_output=True, text=True, check=False)
                    if sub_res.returncode == 0 and "Submitted batch job" in sub_res.stdout:
                        job_id = sub_res.stdout.strip().split()[-1]
                    else:
                        raise RuntimeError(f"Slurm sbatch submission failed: {sub_res.stderr}")

            # 3. Monitor Loop
            if job_id:
                while True:
                    q_check = subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, f"squeue -j {job_id} -h -o '%t'"], capture_output=True, text=True, check=False)
                    job_state = q_check.stdout.strip().split()[0] if q_check.stdout.strip() else ""
                    if not job_state:
                        break
                    time.sleep(900)

                # Final Verification
                final_done = subprocess.check_output(check_cmd, text=True).strip()
                if final_done != "YES":
                    raise RuntimeError(f"HPC Job {job_id} ended without producing {prefix}.gro")

        # 4. Sync & Retrieve Results
        retrieve_cmd = f"""
            SSH_TARGET="{ssh_target}"
            REMOTE_DIR="{remote_dir}"
            LOCAL_DIR="{local_dir}"
            LOG_FILE="{log_path}"

            ssh -o BatchMode=yes "$SSH_TARGET" "cd '$REMOTE_DIR' && tar -czf md_results.tar.gz '{prefix}'.*" >> "$LOG_FILE" 2>&1
            rsync -avz "$SSH_TARGET:$REMOTE_DIR/md_results.tar.gz" "$LOCAL_DIR/" >> "$LOG_FILE" 2>&1
            mkdir -p "$LOCAL_DIR/md_results"
            tar -xzf "$LOCAL_DIR/md_results.tar.gz" -C "$LOCAL_DIR/md_results"
            rm -f "$LOCAL_DIR/md_results.tar.gz"
            ssh -o BatchMode=yes "$SSH_TARGET" "rm -f '$REMOTE_DIR/md_results.tar.gz'" >> "$LOG_FILE" 2>&1
        """
        subprocess.run(retrieve_cmd, shell=True, check=True)

        with open(output_done_path, "w") as f:
            f.write(f"MD simulation completed for {target_id}.\n")

    finally:
        logger.removeHandler(fh)
        fh.close()