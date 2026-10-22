#!/usr/bin/env python3

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, PPBuilder, Select

try:
    import pyrosetta
except ImportError as exc:
    raise SystemExit(
        "PyRosetta is required. Install it in the same Python environment used to run this script."
    )

AA_1_TO_3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}


class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


def normalize_residue_name(residue_name: Optional[str]) -> str:
    if residue_name is None:
        return ""
    residue_name = str(residue_name).strip().upper()
    if not residue_name:
        return ""
    if len(residue_name) == 1:
        return residue_name
    if len(residue_name) == 3:
        return protein_letters_3to1.get(residue_name, "X")
    if residue_name in AA_1_TO_3:
        return residue_name
    return residue_name[0]


def load_residue_table(path: Optional[Path]):
    if path is None or not path.exists():
        return None
    if path.suffix.lower() != ".json":
        return None
    with path.open() as handle:
        payload = json.load(handle)
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        for key in ("residue_table", "rows", "data"):
            value = payload.get(key)
            if isinstance(value, list):
                return value
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fix a protein structure by filling missing residues using PyRosetta."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the input structure file (PDB or mmCIF).",
    )
    parser.add_argument(
        "--sequence",
        "-s",
        required=True,
        help="Path to the extracted FASTA sequence for the target protein.",
    )
    parser.add_argument(
        "--residue-table",
        "-t",
        required=True,
        help="Path to the reconstructed residue-table JSON produced by reconstruct_sequence.py.",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Path to the output cleaned PDB file.",
    )
    return parser.parse_args()


def load_fasta_sequence(fasta_path: Path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No sequence records found in {fasta_path}")
    if len(records) == 1:
        return {records[0].id: str(records[0].seq)}
    return {rec.id: str(rec.seq) for rec in records}


def convert_structure_to_pdb(input_path: Path, output_path: Path) -> Path:
    if input_path.suffix.lower() != ".cif":
        return input_path

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(input_path.stem, str(input_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), NonHetSelect())
    return output_path


def parse_pdb_label(pdb_label: Optional[str], previous_chain: str, previous_resi: int) -> tuple[str, int]:
    if not pdb_label:
        return previous_chain, previous_resi + 1 if previous_resi >= 0 else 1

    parts = pdb_label.strip().split()
    if not parts:
        return previous_chain, previous_resi + 1 if previous_resi >= 0 else 1

    chain = None
    resi = None
    for part in parts:
        try:
            resi = int(part)
            break
        except ValueError:
            continue

    for part in parts:
        try:
            int(part)
        except ValueError:
            chain = part
            break

    if chain is None:
        chain = previous_chain
    if resi is None:
        resi = previous_resi + 1 if previous_resi >= 0 else 1

    return chain, resi


def build_gap_list(pose) -> list[tuple[str, int, int, int, int]]:
    gaps = []
    pose2pdb = pose.pdb_info().pose2pdb
    previous_resi = -1
    previous_chain = ""

    for residue in pose.residues:
        seqpos = residue.seqpos()
        pdb_label = pose2pdb(seqpos)
        chain, resi = parse_pdb_label(pdb_label, previous_chain, previous_resi)

        if residue.is_ligand() or residue.is_metal():
            previous_resi = -1
            previous_chain = ""
            continue

        if previous_chain and chain != previous_chain:
            previous_resi = -1
            previous_chain = chain
            continue

        if previous_resi >= 0 and resi != previous_resi + 1:
            gaps.append((chain, previous_resi + 1, resi - 1, seqpos, previous_resi))

        previous_resi = resi
        previous_chain = chain

    return gaps


def get_gap_sequence(rows: list[dict[str, Any]], chain: str, start: int, end: int) -> str:
    def row_residue(row: dict[str, Any]) -> str:
        for key in ("observed_residue", "recorded_residue", "recorded", "real_residue", "realAA", "residue"):
            candidate = row.get(key)
            if candidate is None:
                continue
            normalized = normalize_residue_name(candidate)
            if normalized and normalized != "X":
                return normalized
        return ""

    sequence = ""
    for row in rows:
        row_chain = str(row.get("chainID") or row.get("chain_id") or "A")
        if row_chain != chain:
            continue
        seqnum = row.get("seqnum")
        if seqnum is None:
            continue
        try:
            seqnum = int(seqnum)
        except (TypeError, ValueError):
            continue
        if start <= seqnum <= end:
            sequence += row_residue(row)
    return sequence


def insert_missing_residues(pose, residue_rows):
    gaps = build_gap_list(pose)
    if not gaps:
        print("No gaps detected in the input pose.")
        return pose

    cm = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    rts = cm.residue_type_set("fa_standard")
    rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
    rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue

    for chain, start, end, seqpos_after_gap, previous_resi in reversed(gaps):
        sequence = get_gap_sequence(residue_rows, chain, start, end)
        if not sequence:
            print(f"Missing sequence not found for gap {chain}:{start}-{end}")
            continue

        previous_pose = pose.pdb_info().pdb2pose(chain, previous_resi)
        if previous_pose is None or previous_pose <= 0:
            print(f"Unable to map previous residue {chain}:{previous_resi} to pose index")
            continue

        insert_pos = previous_pose
        for one_letter in sequence:
            residue_type = rts.get_representative_type_name1(one_letter)
            residue = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(residue_type)
            pose.append_polymer_residue_after_seqpos(residue, insert_pos, True)
            insert_pos += 1

        for new_pos in range(previous_pose + 1, previous_pose + 1 + len(sequence)):
            rm_lower(pose.conformation(), new_pos)
            rm_upper(pose.conformation(), new_pos)

        npos = previous_pose + len(sequence)
        loop = pyrosetta.rosetta.protocols.loops.Loop(previous_pose, npos + 1, npos)
        loops = pyrosetta.rosetta.protocols.loops.Loops()
        loops.add_loop(loop)
        lm = pyrosetta.rosetta.protocols.loop_modeler.LoopModeler()
        lm.set_loops(loops)
        lm.enable_centroid_stage()
        lm.enable_fullatom_stage()
        lm.enable_build_stage()
        lm.apply(pose)

        for new_pos in range(previous_pose + 1, previous_pose + 1 + len(sequence)):
            rm_lower(pose.conformation(), new_pos)
            rm_upper(pose.conformation(), new_pos)

        npos = previous_pose + len(sequence)
        loop = pyrosetta.rosetta.protocols.loops.Loop(previous_pose, npos + 1, npos)
        loops = pyrosetta.rosetta.protocols.loops.Loops()
        loops.add_loop(loop)
        lm = pyrosetta.rosetta.protocols.loop_modeler.LoopModeler()
        lm.set_loops(loops)
        lm.enable_centroid_stage()
        lm.enable_fullatom_stage()
        lm.enable_build_stage()
        lm.apply(pose)

    return pose


#def fix_structure_with_pyrosetta(input_path: Path, fasta_path: Path, residue_table_path: Path, output_pdb: Path):
def fix_structure_with_pyrosetta(input_path: Path, residue_table_path: Path, output_pdb: Path):
    if not input_path.exists():
        raise FileNotFoundError(f"Input structure file not found: {input_path}")
    #if not fasta_path.exists():
    #    raise FileNotFoundError(f"Sequence file not found: {fasta_path}")
    if not residue_table_path.exists():
        raise FileNotFoundError(f"Residue table file not found: {residue_table_path}")

    #sequence_records = load_fasta_sequence(fasta_path)
    #print(f"Loaded sequence records from {fasta_path}: {list(sequence_records.keys())}")

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_handle:
        temp_pdb_path = Path(tmp_handle.name)

    try:
        structure_path = convert_structure_to_pdb(input_path, temp_pdb_path)
        if not structure_path.exists():
            raise FileNotFoundError(f"Structure file not found after conversion: {structure_path}")

        pyrosetta.init(extra_options="-mute all")
        residue_rows = load_residue_table(residue_table_path)
        if residue_rows is None:
            raise ValueError(f"Unable to load residue table from {residue_table_path}")

        pose = pyrosetta.pose_from_pdb(str(structure_path))
        repaired_pose = insert_missing_residues(pose, residue_rows)

        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        repaired_pose.dump_pdb(str(output_pdb))
        print(f"Saved PyRosetta-repaired structure to: {output_pdb}")
    finally:
        if temp_pdb_path.exists():
            temp_pdb_path.unlink(missing_ok=True)


def main():
    args = parse_args()
    input_path = Path(args.input)
    #seq_path = Path(args.sequence)
    residue_table_path = Path(args.residue_table)
    output_pdb = Path(args.output)

    try:
        fix_structure_with_pyrosetta(input_path, residue_table_path, output_pdb)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
