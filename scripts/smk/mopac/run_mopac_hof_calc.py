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
import shutil

load_dotenv()

sys.path.insert(0, os.path.abspath("scripts/smk/sims"))
from hpc_engine import HPCJobManager


def parse_heat_of_formation(arc_filepath):
	"""
	Extracts Heat of Formation in KJ/MOL from MOPAC .arc files.
	Normalizes non-breaking spaces and returns None if missing/unparsed.
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


def parse_length_ns(full_protocol):
	"""Extracts protocol duration in ns from strings containing '100ns' or '1us'."""
	match = re.search(r'(\d+(?:\.\d+)?)\s*(us|ns)', full_protocol, re.IGNORECASE)
	if match:
		val, unit = float(match.group(1)), match.group(2).lower()
		return val * 1000.0 if unit == "us" else val

	raise ValueError(f"Could not parse simulation length from protocol string '{full_protocol}'")


def verify_local_tar(tar_path: str) -> bool:
	"""Verifies that a downloaded tarball is readable and uncorrupted."""
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
	if len(sys.argv) < 4:
		print("Usage: python run_mopac_hof_calc.py <job_description> <output_csv> <full_protocol>")
		sys.exit(1)

	job_description = os.path.abspath(sys.argv[1])
	output_csv = os.path.abspath(sys.argv[2])
	full_protocol = sys.argv[3]

	LENGTH_NS = parse_length_ns(full_protocol)

	local_dir = os.path.dirname(output_csv)
	os.makedirs(local_dir, exist_ok=True)

	log_path = os.path.join(local_dir, "run_mopac_hof_calc.log")
	output_tar = os.path.join(local_dir, "mopac_results.tar.gz")

	p = Path(job_description)
	parts = p.parts

	try:
		idx_mopac = parts.index("mopac")
		
		# Strip trailing filename and structural subfolders ('frames', 'JOB')
		subparts = list(parts[idx_mopac + 1:])
		if subparts[-1].endswith(".txt") or subparts[-1].endswith(".job"):
			subparts.pop()
		if subparts and subparts[-1] in ("frames", "JOB"):
			subparts.pop()

		pdb = subparts[0]
		source = subparts[1]
		protocol = subparts[-1]
		model_id = "/".join(subparts[2:-1])
	except (ValueError, IndexError) as e:
		print(f"Error parsing path wildcards from {job_description}: {e}")
		sys.exit(1)

	target_id = f"{pdb}/{source}/{model_id}/{protocol}"
	prefix = f"{pdb}_{source}_{model_id.replace('/', '_')}_{protocol}"

	# Setup Logger
	logger = logging.getLogger(f"mopac_{prefix}")
	logger.setLevel(logging.INFO)
	logger.handlers.clear()

	fmt = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
	fh = logging.FileHandler(log_path, mode="a"); fh.setFormatter(fmt); logger.addHandler(fh)
	ch = logging.StreamHandler(); ch.setFormatter(fmt); logger.addHandler(ch)

	# Locate input frame tar archive from gmx_em step
	#input_em_tar = os.path.		(f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz")
	input_em_tar = f"results/gmx_em/{pdb}/{source}/{model_id}/{protocol}/em_results/frames/EM_FRAMES_compressed.tar.gz"

	if not os.path.exists(input_em_tar):
		logger.error("❌ Input EM frames archive not found: %s", input_em_tar)
		raise FileNotFoundError(f"Input EM frames archive not found: {input_em_tar}")

	ssh_target = os.getenv("SSH_TARGET")
	hpc_user = os.environ.get("HPC_USER")
	hpc_base = os.environ.get("HPC_REMOTE_BASE", "")
	#clean_base = hpc_base.lstrip("~/")
	remote_dir = os.path.join(hpc_base, 'mopac', pdb, source, model_id, protocol)

	#remote_dir = f"{clean_base}/mopac/{pdb}/{source}/{model_id}/{protocol}"
	logger.info("Initializing HPC MOPAC HoF pipeline for %s...", target_id)

	job_script = os.path.basename(job_description)
	raw_job_name = f"mopac_{prefix}"
	job_name = raw_job_name[:60] if len(raw_job_name) > 60 else raw_job_name
	#job_name = f"mopac_{prefix}"

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

	def prepare_remote_archive():
		"""Atomic compression on remote ONLY when MOPAC_COMPLETE.txt exists."""
		pack_cmd = (
			f"ssh -o BatchMode=yes {ssh_target} \""
			f"cd '{remote_dir}' && "
			f"if [ -f MOPAC_COMPLETE.txt ]; then "
			f"  if [ ! -f mopac_results.tar.gz ]; then "
			f"    tar -czf mopac_results.tar.gz.tmp MOPAC_STAGED && mv mopac_results.tar.gz.tmp mopac_results.tar.gz; "
			f"  fi; "
			f"  echo YES; "
			f"else "
			f"  echo NO; "
			f"fi\""
		)
		res = subprocess.run(pack_cmd, shell=True, capture_output=True, text=True, check=False)
		return "YES" in res.stdout

	# 1. Check if complete on HPC
	if not prepare_remote_archive():
		logger.info("Staging inputs and submitting Slurm job on HPC...")
		subprocess.run(["ssh", "-o", "BatchMode=yes", ssh_target, f"mkdir -p '{remote_dir}'"], check=True)
		rsync_cmd = f"rsync -e 'ssh -o BatchMode=yes' -avz '{input_em_tar}' '{job_description}' '{ssh_target}:{remote_dir}/'"
		subprocess.run(rsync_cmd, shell=True, check=True)

		unpack_cmd = f"ssh -o BatchMode=yes {ssh_target} \"cd '{remote_dir}' && mkdir -p EM_FRAMES && tar -xzf {os.path.basename(input_em_tar)} -C EM_FRAMES --strip-components=1 2>/dev/null || tar -xzf {os.path.basename(input_em_tar)} -C EM_FRAMES\""
		subprocess.run(unpack_cmd, shell=True, check=True)

		job_id = hpc.get_active_job_id()
		if not job_id:
			job_id = hpc.submit_job(unpack_job_archive=False)

		hpc.monitor_job(job_id, completion_check_file="MOPAC_COMPLETE.txt", poll_interval_sec=120)

		if not prepare_remote_archive():
			raise RuntimeError("MOPAC job finished on HPC, but MOPAC_COMPLETE.txt or archive generation failed.")

	# 2. Pull results archive locally
	logger.info("📦 Pulling MOPAC result archive locally...")
	hpc.pull_results(target_subdir="")

	# Clean up any un-tarred directories accidentally copied by rsync
	for stray_dir in ["MOPAC_STAGED", "EM_FRAMES", "FRAMES"]:
		stray_path = os.path.join(local_dir, stray_dir)
		if os.path.exists(stray_path) and os.path.isdir(stray_path):
			shutil.rmtree(stray_path)

	if not os.path.exists(output_tar):
		alt_tar = os.path.join(os.path.dirname(local_dir), "mopac_results.tar.gz")
		if os.path.exists(alt_tar):
			os.rename(alt_tar, output_tar)
		else:
			raise FileNotFoundError(f"Expected retrieved archive missing: {output_tar}")

	# 3. Verify Local Archive Integrity
	if not verify_local_tar(output_tar):
		os.remove(output_tar)
		raise ValueError(f"Downloaded archive {output_tar} is corrupted or incomplete. Removed corrupt file.")

	# 4. Extract .arc files locally and compile CSV
	logger.info("Parsing Heat of Formation data from MOPAC .arc outputs...")
	missing_hof_files = []
	parsed_records = []

	with tempfile.TemporaryDirectory() as temp_dir:
		with tarfile.open(output_tar, "r:gz") as tar:
			for member in tar.getmembers():
				if member.name.endswith(".arc"):
					member.name = os.path.basename(member.name)
					tar.extract(member, path=temp_dir, filter='data')

		arc_files = glob.glob(os.path.join(temp_dir, "*.arc"))

		if not arc_files:
			logger.error("❌ No .arc files found inside %s", output_tar)
			raise FileNotFoundError(f"No .arc files found in {output_tar}")

		for arc_path in sorted(arc_files, key=extract_frame_num):
			fn = os.path.basename(arc_path)
			frame_num = extract_frame_num(arc_path)

			if frame_num < 0:
				logger.error("❌ Could not extract frame index from filename: %s", fn)
				missing_hof_files.append(fn)
				continue

			hof_val = parse_heat_of_formation(arc_path)

			if hof_val is None:
				logger.error("❌ Heat of Formation NOT FOUND or UNPARSABLE in %s", fn)
				missing_hof_files.append(fn)
			else:
				parsed_records.append((frame_num, hof_val))

	if missing_hof_files:
		err_msg = f"MOPAC Heat of Formation extraction failed for {len(missing_hof_files)} file(s): {', '.join(missing_hof_files)}"
		logger.error("💥 CRITICAL FAILURE: %s", err_msg)
		raise RuntimeError(err_msg)

	# 5. Compute Frame Time Scale
	min_frame = min(rec[0] for rec in parsed_records)
	max_frame = max(rec[0] for rec in parsed_records)

	frame_range = (max_frame - min_frame) if min_frame > 0 else max_frame
	if frame_range <= 0:
		raise ValueError(f"Invalid frame range ({frame_range}). Cannot compute time scale.")

	time_step_ns = LENGTH_NS / frame_range
	logger.info("Frame range detected: [%d - %d] | Time scale: %.4f ns/frame", min_frame, max_frame, time_step_ns)

	data_hof = [
		{
			"frame": frame_num,
			"time_ns": round((frame_num - min_frame) * time_step_ns, 4),
			"hof": hof_val
		}
		for frame_num, hof_val in parsed_records
	]

	# 6. Save CSV sorted by frame
	df = pd.DataFrame(data_hof)
	df = df.sort_values(by="frame").reset_index(drop=True)
	df.to_csv(output_csv, index=False)

	logger.info("✅ Heat of Formation CSV successfully written to %s (%d frames validated)", output_csv, len(df))


if __name__ == "__main__":
	main()