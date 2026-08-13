#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    return parser.parse_args()

def build_full_sequence(metadata):
    seq = None

    # Prefer canonical sequence
    if metadata.get("sequence_can"):
        seq = metadata["sequence_can"][0]
    elif metadata.get("sequence"):
        seq = metadata["sequence"][0]

    if not seq:
        raise ValueError("No usable polymer sequence found in metadata")

    # Clean up common mmCIF artifacts
    seq = seq.replace("\n", "").replace(" ", "")

    # Replace non-standard residue codes with parent residues
    normalized = []
    for res in seq:
        # handle multi-letter codes if present
        if len(res) > 1:
            normalized.append(res)
        else:
            normalized.append(res)

    # For this use case, convert known non-standard 3-letter-ish tokens
    # by using the metadata values directly if they are already plain one-letter strings.
    # If you want stronger normalization, expand this logic to handle 3-letter tokens.
    return "".join(MOD_MAP.get(res, res) for res in seq)

def main():
    args = parse_args()
    metadata = json.loads(Path(args.input).read_text())

    seq = build_full_sequence(metadata)

    record = SeqRecord(
        Seq(seq),
        id=metadata.get("pdb", "unknown"),
        description="reconstructed full sequence for modeling"
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([record], args.output, "fasta")

    print(f"Wrote reconstructed sequence to {args.output}")

if __name__ == "__main__":
    main()