#!/usr/bin/env python3
import sys
import os
import logging
import subprocess
import tarfile
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

sys.path.insert(0, os.path.abspath("scripts/smk/sims"))
from hpc_engine import HPCJobManager


def verify_local_tar(tar_path: str) -> bool:
	"""Verifies that a local gzip tarball is uncorrupted and readable."""
	if not os.path.exists(tar_path) or os.path.getsize(tar_path) == 0:
		return False
	try:
		with tarfile.open(tar_path, "r:gz") as tar:
			for _ in tar:
				pass
		return True
	except (tarfile.TarError, EOFError, OSError):
		return False


def main():
	if len(sys.argv) < 3:
		print("Usage: python run_HPC_gmx_em.py <job_description> <em_completed_sentinel>")
		sys.exit(1)

	job_description = os.path.abspath(sys.argv[1])
	em_completed_sentinel = os.path.abspath(sys.argv[2])

	local_dir = os.path.dirname(em_completed_sentinel)
	log_path = os.path.join(local_dir, "run_HPC_gmx_em.log")

	p = Path(job_description)
	parts = p.parts

	try:
		idx = parts.index("gmx_em")
		protocol_idx = parts.index("JOB") - 1 if "JOB" in parts else -2
		protocol = parts[protocol_idx]
		pdb = parts[idx + 1]
		source = parts[idx + 2]
		model_id = "/".join(parts[idx + 3 : protocol_idx])
	except (ValueError, IndexError) as e:
		print(f"Failed to parse wildcards from job description path {job_description}: {e}")
		sys.exit(1)

	target_id = f"{pdb}/{source}/{model_id}/{protocol}"
	prefix = f"{pdb}_{source}_{model_id.replace('/', '_')}_{protocol}"

	logger = logging.getLogger(f"em_{prefix}")
	logger.setLevel(logging.INFO)
	logger.handlers.clear()

	fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
	fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
	ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

	input_tar = os.path.abspath(
		f"results/gromacs/{pdb}/{source}/{model_id}/{protocol}/md_results/frames/FRAMES_compressed.tar.gz"
	)

	ssh_target = os.getenv("SSH_TARGET")
	hpc_user = os.getenv("HPC_USER")
	hpc_base = os.getenv("HPC_REMOTE_BASE")
	#clean_base = hpc_base.lstrip("~/")
	#remote_dir = f"{clean_base}/gmx_em/{pdb}/{source}/{model_id}/{protocol}/"
	remote_dir = os.path.join(hpc_base, 'gmx_em', pdb, source, model_id, protocol)

	logger.info("Initializing HPC EM pipeline for %s...", target_id)
	logger.info("Remote Directory: %s", remote_dir)

	job_script = os.path.basename(job_description)
	raw_job_name = f"em_{prefix}"
	job_name = raw_job_name[:60] if len(raw_job_name) > 60 else raw_job_name
	#job_name = f"em_{prefix}"

	hpc = HPCJobManager(
		ssh_target=ssh_target,
		remote_dir=remote_dir,
		local_dir=local_dir,
		job_name=job_name,
		job_script=job_script,
		log_path=log_path,
		logger=logger,
		hpc_user=hpc_user,
	)

	def prepare_remote_archive():
		"""Atomic compression on remote ONLY if completion marker (EM_COMPLETE.txt) exists."""
		pack_cmd = (
			f"ssh -o BatchMode=yes {ssh_target} "
			f"\"cd '{remote_dir}' && "
			f"if [ -f EM_COMPLETE.txt ]; then "
			f"  if [ ! -f EM_FRAMES_compressed.tar.gz ]; then "
			f"    tar -czf EM_FRAMES_compressed.tar.gz.tmp EM_FRAMES && mv EM_FRAMES_compressed.tar.gz.tmp EM_FRAMES_compressed.tar.gz; "
			f"  fi; "
			f"  echo YES; "
			f"else "
			f"  echo NO; "
			f"fi\""
		)
		res = subprocess.run(pack_cmd, shell=True, capture_output=True, text=True, check=False, env=os.environ)
		return "YES" in res.stdout

	# 1. Check if complete on HPC
	if not prepare_remote_archive():
		# Setup and submit HPC job
		subprocess.run(
			["ssh", "-o", "BatchMode=yes", ssh_target, f"mkdir -p '{remote_dir}'"],
			check=True,
			env=os.environ,
		)
		rsync_cmd = f"rsync -e 'ssh -o BatchMode=yes' -avz '{input_tar}' '{job_description}' '{ssh_target}:{remote_dir}/'"
		subprocess.run(rsync_cmd, shell=True, check=True, env=os.environ)

		unpack_cmd = (
			f"ssh -o BatchMode=yes {ssh_target} "
			f"\"cd '{remote_dir}' && [ ! -d FRAMES ] && tar -xzf FRAMES_compressed.tar.gz || true\""
		)
		subprocess.run(unpack_cmd, shell=True, check=True, env=os.environ)

		job_id = hpc.get_active_job_id()
		if not job_id:
			job_id = hpc.submit_job(unpack_job_archive=False)

		hpc.monitor_job(job_id, completion_check_file="EM_COMPLETE.txt", poll_interval_sec=120)

		if not prepare_remote_archive():
			raise RuntimeError("HPC Energy Minimization finished, but EM_COMPLETE.txt or archive generation failed on remote.")

	# 2. Pull the compressed archive locally
	local_frames_dir = os.path.join(local_dir, "frames")
	os.makedirs(local_frames_dir, exist_ok=True)
	local_tar_path = os.path.join(local_frames_dir, "EM_FRAMES_compressed.tar.gz")

	pull_cmd = f"rsync -e 'ssh -o BatchMode=yes' -avz '{ssh_target}:{remote_dir}/EM_FRAMES_compressed.tar.gz' '{local_tar_path}'"
	res = subprocess.run(pull_cmd, shell=True, executable="/bin/bash")
	if res.returncode not in (0, 24):
		raise subprocess.CalledProcessError(res.returncode, pull_cmd)

	# 3. Strict local integrity check before marking rule complete
	logger.info("Verifying local tarball integrity: %s", local_tar_path)
	if not verify_local_tar(local_tar_path):
		os.remove(local_tar_path)  # Delete corrupted download
		raise ValueError(f"Downloaded archive {local_tar_path} failed gzip integrity check! Removed corrupt file.")

	# 4. Write completion sentinel
	with open(em_completed_sentinel, "w") as f:
		f.write(f"Energy minimization completed for {target_id}.\n")

	logger.info("✅ Energy Minimization workflow complete for %s", target_id)


if __name__ == "__main__":
	main()