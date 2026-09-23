import os
import time
import logging
import subprocess
import yaml
import shutil

def get_simulation_progress(ssh_target: str, remote_dir: str, job_id: str = None) -> str | None:
    """
    Parses live progress (Step, Time, ns/day, performance, estimated finish)
    from the current Slurm .err or GROMACS .log output on HPC.
    """
    if job_id:
        # Match Slurm stderr file specifically using job_id
        cmd = (
            f"ssh -o BatchMode=yes {ssh_target} "
            f"\"[ -f '{remote_dir}'/*_{job_id}.err ] && "
            f"grep -E 'Step|Vol|Time|ns/day|Hours' '{remote_dir}'/*_{job_id}.err | tail -n 10 || echo ''\""
        )
    else:
        # Fallback to latest md log
        cmd = (
            f"ssh -o BatchMode=yes {ssh_target} "
            f"\"[ -f '{remote_dir}'/*.log ] && "
            f"grep -E 'Step|Vol|Time|ns/day|Hours' '{remote_dir}'/*.log | tail -n 10 || echo ''\""
        )

    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=15)
        return res.stdout.strip() if res.stdout else None
    except Exception:
        return None
    


class HPCJobManager:
    """
    Manages remote execution on HPC clusters (e.g., Komondor via Slurm).
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
        self.local_status_file = os.path.join(self.local_dir, "job_execution_status.log")

    def _run_ssh(self, command: str, capture_output=True, check=False) -> subprocess.CompletedProcess:
        cmd = ["ssh", "-o", "BatchMode=yes", self.ssh_target, command]
        return subprocess.run(cmd, capture_output=capture_output, text=True, check=check)

    def is_remote_finished(self, completion_check_file: str) -> bool:
        check_cmd = f"[ -f '{self.remote_dir}/{completion_check_file}' ] && echo YES || echo NO"
        res = self._run_ssh(check_cmd)
        return res.returncode == 0 and "YES" in res.stdout

    def get_active_job_id(self) -> str | None:
        squeue_cmd = f"squeue -u {self.hpc_user} -n {self.job_name} -h -o '%i %t'"
        res = self._run_ssh(squeue_cmd)
        q_out = res.stdout.strip()
        if q_out:
            return q_out.split()[0]
        return None

    def has_checkpoint(self) -> bool:
        """Checks if a valid GROMACS checkpoint exists on HPC to allow resuming."""
        check_cmd = f"ls '{self.remote_dir}'/*.cpt 2>/dev/null | wc -l"
        res = self._run_ssh(check_cmd)
        try:
            return int(res.stdout.strip()) > 0
        except ValueError:
            return False

    def submit_job(self, unpack_job_archive: bool = True) -> str:
        prep_cmd = (
            f"cd '{self.remote_dir}' && [ -f JOB.tar.gz ] && tar -xzf JOB.tar.gz; sbatch {self.job_script}"
            if unpack_job_archive
            else f"cd '{self.remote_dir}' && sbatch {self.job_script}"
        )
        res = self._run_ssh(prep_cmd)
        
        if res.returncode == 0 and "Submitted batch job" in res.stdout:
            job_id = res.stdout.strip().split()[-1]
            self.logger.info("✅ Slurm job submitted! Assigned Job ID: %s", job_id)
            return job_id
        else:
            self.logger.error("❌ Failed to submit Slurm job: %s", res.stderr)
            raise RuntimeError(f"Slurm sbatch submission failed: {res.stderr}")

    def pull_results(self, target_subdir: str = "md_results", excludes: list = None):
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

    def monitor_job(self, job_id: str, completion_check_file: str, poll_interval_sec: int = 300, get_progress_fn=None):
        """
        Polls Slurm queue. If the job exits due to Slurm walltime limit but a 
        checkpoint exists, automatically resubmits the job until completion.
        """
        self.logger.info("Monitoring Slurm Job %s (Polling interval: %ds)...", job_id, poll_interval_sec)

        while True:
            if str(job_id).isdigit():
                q_check = self._run_ssh(f"squeue -j {job_id} -h -o '%t'")
                stdout_clean = q_check.stdout.strip()
                job_state = stdout_clean.split()[0] if stdout_clean else ""

                if not job_state:
                    self.logger.info("Job %s left the Slurm queue. Verifying completion...", job_id)
                    time.sleep(10)  # Brief buffer for filesystem sync

                    if self.is_remote_finished(completion_check_file):
                        self.logger.info("🎉 Simulation finished successfully on HPC!")
                        break

                    if self.has_checkpoint():
                        self.logger.warning("⏰ Job %s exited (likely hit 7-day walltime limit), but checkpoint file (.cpt) found.", job_id)
                        self.logger.info("🔄 Automatically resubmitting Slurm job to resume from checkpoint...")
                        job_id = self.submit_job(unpack_job_archive=False)
                        self.logger.info("Attached monitor to new Slurm Job ID: %s", job_id)
                        continue
                    else:
                        self.logger.error("❌ Job %s left queue without producing %s or a valid .cpt file.", job_id, completion_check_file)
                        break

                timestamp = time.strftime("%Y-%m-%d %H:%M:%S")

                if job_state == "R":
                    progress = (
                        get_progress_fn(self.ssh_target, self.remote_dir, job_id)
                        if callable(get_progress_fn)
                        else None
                    )

                    status_msg = f"[{timestamp}][JOB {job_id} RUNNING]\n"
                    if progress:
                        status_msg += f"{progress}\n"
                    else:
                        status_msg += "Waiting for trajectory metrics output...\n"

                    self.logger.info("Heartbeat Job %s: Running", job_id)
                    with open(self.local_status_file, "a") as f:
                        f.write(status_msg + "-" * 50 + "\n")

                else:
                    status_msg = f"[{timestamp}][JOB {job_id} QUEUED] State: {job_state}\n"
                    self.logger.info("Heartbeat Job %s: Queued (%s)", job_id, job_state)
                    with open(self.local_status_file, "a") as f:
                        f.write(status_msg + "-" * 50 + "\n")

            time.sleep(poll_interval_sec)

    def execute_pipeline(
        self,
        completion_check_file: str,
        target_subdir: str = "md_results",
        poll_interval_sec: int = 300,
        get_progress_fn=None,
        unpack_job_archive: bool = True
    ):
        local_finalized_sentinel = os.path.join(self.local_dir, "md_results", "md_completed.txt")
        if os.path.exists(local_finalized_sentinel):
            self.logger.info("Local finalized results found. Bypassing HPC monitor.")
            return

        self._run_ssh(f"mkdir -p '{self.remote_dir}'")

        if self.is_remote_finished(completion_check_file):
            self.logger.info("🎉 Task already finished on HPC (%s present)", completion_check_file)
        else:
            job_id = self.get_active_job_id()
            if job_id:
                self.logger.info("⏳ Job '%s' already active (ID: %s). Attaching monitor...", self.job_name, job_id)
            else:
                job_id = self.submit_job(unpack_job_archive=unpack_job_archive)

            if job_id:
                self.monitor_job(
                    job_id=job_id,
                    completion_check_file=completion_check_file,
                    poll_interval_sec=poll_interval_sec,
                    get_progress_fn=get_progress_fn
                )

            if not self.is_remote_finished(completion_check_file):
                raise RuntimeError(f"HPC execution ended without producing {completion_check_file}.")

        self.logger.info("📦 Synchronizing results from HPC...")
        self.pull_results(target_subdir=target_subdir)