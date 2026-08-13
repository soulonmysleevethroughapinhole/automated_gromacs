#!/usr/bin/env python3
"""Compute CA-based RMSD between an empirical structure and multiple model files.

Usage:
  python scripts/screen_rmsd.py --input-emp <empirical.pdb> --output-dir <dir> --threshold 1.0 model1.pdb model2.pdb ...

The script writes one JSON report per model into <output-dir>.
"""

import argparse
import json
import os
import subprocess
import tempfile
import numpy as np
from pymol import cmd


def compute_rmsd(a, b):
	if a.size == 0 or b.size == 0:
		return None
	n = min(len(a), len(b))
	a_c = a[:n] - np.mean(a[:n], axis=0)
	b_c = b[:n] - np.mean(b[:n], axis=0)
	return float(np.sqrt(np.sum((a_c - b_c) ** 2) / n))


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input-emp', required=True, help='Empirical structure (PDB or CIF)')
	parser.add_argument('--output-dir', required=True, help='Directory to write per-model JSON reports')
	parser.add_argument('--threshold', type=float, default=1.0, help='RMSD threshold for divergence')
	parser.add_argument('models', nargs='+', help='Model files to evaluate (PDB or CIF)')

	args = parser.parse_args()

	os.makedirs(args.output_dir, exist_ok=True)

	#emp_coords = get_ca_coords(args.input_emp)
	cmd.load(args.input_emp, object='empirical')
	
	for i, model in enumerate(args.models):
		#mod_coords = get_ca_coords(model)
		cmd.load(model, object=i)
		#rmsd = compute_rmsd(emp_coords, mod_coords)
		rmsd = cmd.align(f'{i} and name CA', 'empirical and name CA')[0]
		cmd.delete(i) # to save memoriy

		status = 'unknown'
		if rmsd is None:
			status = 'no_coords'
		else:
			status = 'divergent' if rmsd >= args.threshold else 'too_similar'

		model_name = os.path.basename(model)
		if model_name.lower().endswith('.pdb') or model_name.lower().endswith('.cif'):
			model_name = os.path.splitext(model_name)[0]

		report = {
			'input_empirical': args.input_emp,
			'input_model': model,
			'model': model_name,
			'rmsd': rmsd,
			'structural_evaluation': status,
		}

		out_file = os.path.join(args.output_dir, f"{model_name}.json")
		with open(out_file, 'w') as fh:
			json.dump(report, fh, indent=4)

	cmd.delete('empirical')

if __name__ == '__main__':
	main()


#def get_ca_coords(path):
#	# If CIF provided, convert to PDB using existing helper script
#	cleanup_tmp = False
#	pdb_path = path
#	if path.lower().endswith('.cif'):
#		tmpf = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
#		tmpf.close()
#		subprocess.check_call([
#			'python', 'scripts/convert_cif_to_pdb.py',
#			'--input', path,
#			'--output', tmpf.name,
#		])
#		pdb_path = tmpf.name
#		cleanup_tmp = True
#
#	coords = []
#	with open(pdb_path, 'r') as fh:
#		for line in fh:
#			if (line.startswith('ATOM') or line.startswith('HETATM')) and line[12:16].strip() == 'CA':
#				try:
#					x = float(line[30:38])
#					y = float(line[38:46])
#					z = float(line[46:54])
#				except ValueError:
#					continue
#				coords.append([x, y, z])
#
#	if cleanup_tmp:
#		try:
#			os.unlink(pdb_path)
#		except Exception:
#			pass
#
#	return np.array(coords)
#
#
#
#def parse_args():
#	parser = argparse.ArgumentParser()
#	parser.add_argument("--input", "-i", required=True)
#	parser.add_argument("--output", "-o", required=True)
#	return parser.parse_args()
#
#def get_ca_coords(pdb_path):
#	"""
#	Extracts C-alpha coordinates while ignoring alternate locations (altlocs)
#	and multiple NMR models.
#	"""
#	coords = []
#	seen_residues = set()
#	
#	with open(pdb_path, "r") as f:
#		for line in f:
#			if line.startswith("ENDMDL"): # Stop at first model if multi-model file
#				break
#			if line.startswith("ATOM"):
#				atom_name = line[12:16].strip()
#				alt_loc = line[16].strip()
#				
#				# Grab C-alpha atoms; skip alternate locations (A/B altlocs)
#				if atom_name == "CA" and alt_loc in ("", "A"):
#					chain = line[21]
#					res_num = line[22:26].strip()
#					res_id = f"{chain}_{res_num}"
#					
#					if res_id not in seen_residues:
#						coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
#						seen_residues.add(res_id)
#						
#	return np.array(coords)
#
#
#def calculate_kabsch_rmsd(P, Q):
#	"""
#	Calculates optimal RMSD between two coordinate sets (Nx3) 
#	using the Kabsch SVD algorithm (Translation + Rotation).
#	"""
#	# 1. Translate to origin
#	P_centered = P - np.mean(P, axis=0)
#	Q_centered = Q - np.mean(Q, axis=0)
#	
#	# 2. Compute covariance matrix
#	H = P_centered.T @ Q_centered
#	
#	# 3. Singular Value Decomposition (SVD)
#	U, S, Vt = np.linalg.svd(H)
#	
#	# 4. Compute optimal rotation matrix
#	d = np.linalg.det(Vt.T @ U.T)
#	correction = np.identity(3)
#	correction[2, 2] = np.sign(d) # Handle right/left-handed coordinate reflection
#	
#	R = Vt.T @ correction @ U.T
#	
#	# 5. Rotate P_centered onto Q_centered
#	P_rotated = P_centered @ R.T
#	
#	# 6. Calculate Root Mean Square Deviation
#	diff = P_rotated - Q_centered
#	rmsd = np.sqrt(np.sum(diff ** 2) / len(P))
#	return float(rmsd)
#
##
##checkpoint evaluate_predictions:
##	input:
##		clean_empirical = 'results/structures/{pdb}/empirical/canonical_structure.pdb',
##		models = expand("results/structures/{{pdb}}/{{source}}/model_{model_num}.pdb", model_num=MODELS)
##	output:
##		eval_dir = directory("results/structures/{pdb}/{source}/evaluations")
##	run:
#
#