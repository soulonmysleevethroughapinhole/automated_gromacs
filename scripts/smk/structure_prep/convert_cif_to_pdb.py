import argparse
from pathlib import Path

try:
    from Bio.PDB import MMCIFParser, PDBIO, Select    
except ImportError as exc:
    raise SystemExit(
        "Biopython is required for CIF->PDB conversion. Install it with `pip install biopython` or via conda."
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an mmCIF file to PDB format."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the input CIF file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Path to the output PDB file.",
    )
    return parser.parse_args()

class NonHetSelect(Select):
    def accept_residue(self, residue):
        # residue.id[0] == " " indicates a standard amino acid/nucleotide
        # HETATMs (like ligands or water) have "H_" or "W" in this position
        return 1 if residue.id[0] == " " else 0
    
def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise SystemExit(f"Input file does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser()
    structure_id = input_path.stem
    structure = parser.get_structure(structure_id, str(input_path))

    io = PDBIO()
    io.set_structure(structure)
    #io.save(str(output_path))
    io.save(str(output_path), NonHetSelect())

    

if __name__ == '__main__':
    main()
