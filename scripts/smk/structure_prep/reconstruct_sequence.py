#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import List, Dict, Any, Optional

from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

MOD_MAP = {
    "CSO": "CYS",
    "MSE": "MET",
    "SEC": "CYS",
    "SEP": "SER",
    "TPO": "THR",
    "PTR": "TYR",
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True, help="Path to the input mmCIF file")
    parser.add_argument("--metadata", "-m", required=True, help="Path to the metadata JSON")
    parser.add_argument("--output", "-o", required=True, help="Path to the output FASTA")
    parser.add_argument("--table-out", required=False, help="Optional path for residue table JSON")
    return parser.parse_args()


def get_sequence_from_metadata(metadata: Dict[str, Any]) -> str:
    for key in ("sequence_can", "sequence", "struct_sequence"):
        value = metadata.get(key)
        if not value:
            continue
        if isinstance(value, list):
            value = value[0]
        seq = str(value).replace("\n", "").replace(" ", "")
        if seq:
            return seq
    raise ValueError("No usable sequence found in metadata")


def normalize_residue_name(res_name: str) -> str:
    res_name = res_name.strip().upper()
    if not res_name:
        return "X"
    if res_name in MOD_MAP:
        return MOD_MAP[res_name]
    if len(res_name) == 3:
        return protein_letters_3to1.get(res_name, "X")
    return res_name[0] if len(res_name) == 1 else "X"


def build_residue_table(metadata: Dict[str, Any], structure: PDB.Structure.Structure, chain_id: str = "A") -> List[Dict[str, Any]]:
    full_sequence = get_sequence_from_metadata(metadata)
    full_letters = [aa for aa in full_sequence if aa.isalpha()]

    # Extract missing residue numbers from metadata
    missing_positions = set()
    structural_issues = metadata.get("structural_issues", {}) or {}
    for item in structural_issues.get("missing_residues", []) or []:
        seq_num = item.get("residue_number")
        if seq_num is not None:
            try:
                missing_positions.add(int(seq_num))
            except ValueError:
                pass

    # Pick the chain from the structure
    model = structure[0]
    if chain_id in model:
        chain = model[chain_id]
    else:
        chain = next(iter(model.get_chains()))

    observed = {}
    for residue in chain:
        if not is_aa(residue):
            continue
        observed[residue.id[1]] = residue

    all_seqnums = sorted(set(observed.keys()) | missing_positions)

    rows: List[Dict[str, Any]] = []
    seq_idx = 0

    for seqnum in all_seqnums:
        residue = observed.get(seqnum)
        recorded = None if residue is None else residue.get_resname().strip()
        if seq_idx < len(full_letters):
            real = full_letters[seq_idx]
        else:
            real = "X"

        rows.append({
            "chainID": chain_id,
            "seqnum": seqnum,
            "recorded_residue": recorded,
            "real_residue": real,
        })
        seq_idx += 1

    return rows


def main():
    args = parse_args()

    input_path = Path(args.input)
    metadata_path = Path(args.metadata)
    output_path = Path(args.output)
    table_out_path = Path(args.table_out) if args.table_out else None

    if not input_path.exists():
        raise SystemExit(f"Input mmCIF file not found: {input_path}")
    if not metadata_path.exists():
        raise SystemExit(f"Metadata file not found: {metadata_path}")

    metadata = json.loads(metadata_path.read_text())

    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", str(input_path))

    residue_table = build_residue_table(metadata, structure, chain_id="A")

    # Write FASTA from the canonical full sequence
    full_sequence = get_sequence_from_metadata(metadata)
    record = SeqRecord(
        Seq(full_sequence),
        id=metadata.get("pdb", "unknown"),
        description="reconstructed full sequence for modeling"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([record], output_path, "fasta")

    if table_out_path is not None:
        table_out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(table_out_path, "w") as fh:
            json.dump(residue_table, fh, indent=2)

    print(f"Wrote FASTA to: {output_path}")
    if table_out_path is not None:
        print(f"Wrote residue table to: {table_out_path}")


if __name__ == "__main__":
    main()