import os
import time
import logging
import subprocess
import yaml
import shutil

def get_simulation_progress(ssh_target: str, remote_dir: str) -> str | None:
    """Helper to parse the latest frame/time progress from remote md.log."""
    cmd = (
        f"ssh -o BatchMode=yes {ssh_target} "
        f"\"[ -f '{remote_dir}/*.log' ] && tail -n 20 '{remote_dir}'/*.log | grep 'Step' | tail -n 1 || echo ''\""
    )
    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
        return res.stdout.strip() if res.stdout else None
    except Exception:
        return None

def execute_hpc_md(
    wildcards,
    job_description_path: str,
    output_done_path: str,
    log_path: str,
    hpc_user: str = None,
    hpc_host: str = "komondor",
    hpc_remote_base: str = None
):
    """
    Orchestrates HPC execution and monitoring for custom simulation protocols.
    """
    pdb = wildcards.pdb
    source = wildcards.source
    model_id = wildcards.model_id
    protocol = wildcards.protocol

    # 1. Setup Isolated Logger
    log_abs = os.path.abspath(log_path)
    os.makedirs(os.path.dirname(log_abs), exist_ok=True)

    logger_name = f"run_custom_md_HPC_{pdb}_{source}_{model_id}_{protocol}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(log_abs, mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    target_id = f"{pdb}/{source}/{model_id}/{protocol}"
    prefix = f"{pdb}_{source}_{model_id}_{protocol}_md"
    job_script_filename = f"{prefix}_job.job"

    job_dir = os.path.abspath(os.path.dirname(job_description_path))
    local_protocol_dir = os.path.abspath(os.path.dirname(output_done_path))
    step_log_file = os.path.join(job_dir, "simulation_steps.log")

    logger.info("Initializing HPC custom MD execution phase for target: %s", target_id)

    # 2. Check execution mode in job description
    if not os.path.exists(job_description_path):
        raise FileNotFoundError(f"Job description file not found: {job_description_path}")

    with open(job_description_path, "r") as f:
        job_description_content = f.read().strip()

    if job_description_content == "BLANK":
        logger.info("Target %s is configured for LOCAL mode. Skipping HPC rule.", target_id)
        return

    # 3. Construct Remote Directory
    remote_base = hpc_remote_base or os.environ.get("HPC_REMOTE_BASE", "/home/johnnys/projects/automated_gromacs")
    remote_dir = os.path.join(remote_base, pdb, source, model_id, protocol)
    ssh_target = hpc_host or "komondor"

    # 4. Instantiate HPCJobManager
    manager = HPCJobManager(
        ssh_target=ssh_target,
        remote_dir=remote_dir,
        local_dir=local_protocol_dir,
        job_name=f"md_{pdb}_{source}_{model_id}_{protocol}",
        job_script=job_script_filename,
        log_path=log_abs,
        logger=logger,
        step_log_file=step_log_file,
        hpc_user=hpc_user or os.environ.get("HPC_USER", "johnnys")
    )

    try:
        # 5. Execute Remote Slurm Pipeline
        # Expects {prefix}.gro upon completion
        manager.execute_pipeline(
            completion_check_file=f"{prefix}.gro",
            #target_subdir=".",  # Pull directly into protocol folder
            target_subdir="md_results_HPC",  # Pull directly into protocol folder
            poll_interval_sec=300, # Poll Slurm status every 5 minutes
            get_progress_fn=get_simulation_progress,
            unpack_job_archive=True
        )

        # 6. Write Local Sentinel File
        with open(output_done_path, "w") as f:
            f.write(f"HPC Custom MD completed for {target_id} at {remote_dir}.\n")

        # 7. Mirror the canonical MD artifacts into the protocol md_results directory so downstream rules
        # (trajectory cleanup, analysis, movies) work without a second staging hop.
        results_dir = os.path.join(local_protocol_dir, "md_results")
        os.makedirs(results_dir, exist_ok=True)
        for suffix in [".xtc", ".tpr", ".gro", ".cpt", ".edr", ".log"]:
            src = os.path.join(local_protocol_dir, f"{prefix}{suffix}")
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(results_dir, os.path.basename(src)))
        with open(os.path.join(results_dir, "md_completed.txt"), "w") as f:
            f.write(f"Custom MD protocol '{protocol}' completed for {target_id}.\n")

        # 8. Also mirror the canonical output contract expected by downstream rules
        # in the same protocol directory used by the DAG.
        canonical_done = os.path.join(local_protocol_dir, "md_results", "md_completed.txt")
        if os.path.exists(canonical_done):
            with open(canonical_done, "w") as f:
                f.write(f"Custom MD protocol '{protocol}' completed for {target_id}.\n")

        logger.info("🎉 HPC Custom MD finished successfully for %s", target_id)

    except Exception as err:
        logger.exception("Execution failed in execute_hpc_md for %s: %s", target_id, str(err))
        raise

class HPCJobManager:
    """
    Manages remote execution on HPC clusters (e.g., Komondor via Slurm).
    Handles:
      - Directory creation & job synchronization
      - Checking existing status (Finished / Active / Submitted)
      - Non-blocking polling loops
      - Synchronizing results back locally via rsync
    """
    def __init__(
        self,
        ssh_target: str,
        remote_dir: str,
        local_dir: str,
        job_name: str,
        job_script: str,
        log_path: str,
        logger: logging.Logger = None,
        step_log_file: str = None,
        hpc_user: str = None,
    ):
        self.ssh_target = ssh_target
        self.remote_dir = remote_dir
        self.local_dir = os.path.abspath(local_dir)
        self.job_name = job_name
        self.job_script = job_script
        self.log_path = os.path.abspath(log_path)
        self.logger = logger or logging.getLogger("HPCJobManager")
        self.step_log_file = step_log_file
        self.hpc_user = hpc_user or os.environ.get("USER", "")

    def _log_step(self, stage: str, message: str):
        """Helper to append structured workflow progress steps."""
        if self.step_log_file and callable(globals().get("log_sim_step")):
            log_sim_step(stage, self.step_log_file, message)

    def _run_ssh(self, command: str, capture_output=True, check=False) -> subprocess.CompletedProcess:
        """Executes SSH command with BatchMode enabled to prevent hanging."""
        cmd = ["ssh", "-o", "BatchMode=yes", self.ssh_target, command]
        return subprocess.run(cmd, capture_output=capture_output, text=True, check=check)

    def is_remote_finished(self, completion_check_file: str) -> bool:
        """Checks if a key output file exists remotely (e.g., .gro file or output archive)."""
        check_cmd = f"[ -f '{self.remote_dir}/{completion_check_file}' ] && echo YES || echo NO"
        res = self._run_ssh(check_cmd)
        return res.returncode == 0 and "YES" in res.stdout

    def get_active_job_id(self) -> str | None:
        """Returns the Slurm Job ID if the job is already QUEUED or RUNNING."""
        squeue_cmd = f"squeue -u {self.hpc_user} -n {self.job_name} -h -o '%i %t'"
        res = self._run_ssh(squeue_cmd)
        q_out = res.stdout.strip()
        if q_out:
            return q_out.split()[0]
        return None

    def submit_job(self, unpack_job_archive: bool = True) -> str:
        """Unpacks job archive (if required) and submits the Slurm script."""
        prep_cmd = f"cd '{self.remote_dir}' && [ -f JOB.tar.gz ] && tar -xzf JOB.tar.gz; sbatch {self.job_script}" if unpack_job_archive else f"cd '{self.remote_dir}' && sbatch {self.job_script}"
        res = self._run_ssh(prep_cmd)
        
        if res.returncode == 0 and "Submitted batch job" in res.stdout:
            job_id = res.stdout.strip().split()[-1]
            self.logger.info("✅ Slurm job submitted! Assigned Job ID: %s", job_id)
            self._log_step("SLURM_SUBMIT", f"Submitted job ID: {job_id}")
            return job_id
        else:
            self.logger.error("❌ Failed to submit Slurm job: %s", res.stderr)
            self._log_step("SLURM_FAILED", f"Submission failed: {res.stderr.strip()}")
            raise RuntimeError(f"Slurm sbatch submission failed: {res.stderr}")

    def monitor_job(self, job_id: str, poll_interval_sec: int = 900, get_progress_fn=None):
        """Polls Slurm queue until job leaves queue, logging heartbeats."""
        self.logger.info("Monitoring Slurm Job %s...", job_id)
        last_progress_msg = ""

        while True:
            if str(job_id).isdigit():
                q_check = self._run_ssh(f"squeue -j {job_id} -h -o '%t'")
                stdout_clean = q_check.stdout.strip()
                job_state = stdout_clean.split()[0] if stdout_clean else ""

                if not job_state:
                    self.logger.info("Job %s left the queue. Verifying completion...", job_id)
                    break

                if job_state == "R":
                    progress = get_progress_fn(self.ssh_target, self.remote_dir) if callable(get_progress_fn) else None
                    if progress and progress != last_progress_msg:
                        self._log_step("HEARTBEAT", f"Job {job_id} running - {progress}")
                        last_progress_msg = progress
                    else:
                        self._log_step("HEARTBEAT", f"Job {job_id} actively executing on HPC")
                else:
                    self._log_step("HEARTBEAT", f"Job {job_id} queued (State: {job_state})")

            time.sleep(poll_interval_sec)

    def pull_results(self, target_subdir: str = "md_results", excludes: list = None):
        """Pulls files from remote_dir into local_dir/target_subdir via rsync."""
        if excludes is None:
            excludes = ["JOB", "JOB.tar.gz"]

        local_target_dir = os.path.join(self.local_dir, target_subdir)
        os.makedirs(local_target_dir, exist_ok=True)

        exclude_flags = " ".join([f'--exclude="{ex}"' for ex in excludes])
        
        shell_cmd = f"""
            SSH_TARGET="{self.ssh_target}"
            REMOTE_DIR="{self.remote_dir}"
            LOCAL_TARGET_DIR="{local_target_dir}"
            LOG_FILE="{self.log_path}"

            rsync -e "ssh -o BatchMode=yes" -avz \
                {exclude_flags} \
                "$SSH_TARGET:$REMOTE_DIR/" \
                "$LOCAL_TARGET_DIR/" >> "$LOG_FILE" 2>&1
        """
        subprocess.run(shell_cmd, shell=True, check=True, executable="/bin/bash")
        self._log_step("RETRIEVE_COMPLETE", f"Downloaded HPC outputs into {target_subdir}")

    def execute_pipeline(
        self,
        completion_check_file: str,
        target_subdir: str = "md_results",
        poll_interval_sec: int = 900,
        get_progress_fn=None,
        unpack_job_archive: bool = True
    ):
        """
        Main entrypoint: Orchestrates directory setup, status checks, 
        submission/monitoring, and automatic retrieval.
        """
        # 0. LOCAL OVERRIDE: If local finalized output already exists, exit immediately!
        local_finalized_sentinel = os.path.join(self.local_dir, "md_results", "md_completed.txt")
        if os.path.exists(local_finalized_sentinel):
            self.logger.info("Local finalized results found at %s. Skipping HPC checks/downloads.", local_finalized_sentinel)
            self._log_step("LOCAL_CHECK", "Local finalized results already present. HPC step bypassed.")
            return
        # 1. Ensure remote target directory exists
        self._run_ssh(f"mkdir -p '{self.remote_dir}'")

        # 2. Check if already finished
        if self.is_remote_finished(completion_check_file):
            self.logger.info("🎉 Task already finished on HPC (Found %s)", completion_check_file)
            self._log_step("CHECK_REMOTE", f"Task finished ({completion_check_file} present)")
        else:
            # 3. Check queue for active job ID or submit new job
            job_id = self.get_active_job_id()
            if job_id:
                self.logger.info("⏳ Slurm Job '%s' already active (Job ID: %s). Attaching monitor...", self.job_name, job_id)
                self._log_step("SLURM_ACTIVE", f"Job ID {job_id} running/queued on HPC")
            else:
                self.logger.info("🚀 Submitting new Slurm job %s to %s...", self.job_name, self.ssh_target)
                job_id = self.submit_job(unpack_job_archive=unpack_job_archive)

            # 4. Monitor loop
            if job_id:
                self.monitor_job(job_id, poll_interval_sec=poll_interval_sec, get_progress_fn=get_progress_fn)

            # 5. Verify completion marker file produced
            if not self.is_remote_finished(completion_check_file):
                self._log_step("SLURM_FAILED", f"Job {job_id} ended without producing {completion_check_file}")
                raise RuntimeError(f"HPC Job {job_id} terminated unexpectedly ({completion_check_file} missing).")

            self._log_step("SLURM_FINISHED", f"Job {job_id} completed successfully")

        # 6. Synchronize output artifacts locally
        self.logger.info("📦 Synchronizing results from HPC...")
        self.pull_results(target_subdir=target_subdir)

