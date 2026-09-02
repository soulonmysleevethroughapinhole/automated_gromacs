import os
import time
import logging
import subprocess

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

