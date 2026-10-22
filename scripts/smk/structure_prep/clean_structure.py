import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select

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

STANDARD_PROTEIN_3LETTER = set(AA_1_TO_3.values()).union(
    {'MSE', 'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH'}
)


def normalize_residue_name(residue_name: Optional[str]) -> str:
    if residue_name is None:
        return ""
    name = str(residue_name).strip().upper()
    if not name:
        return ""
    if len(name) == 1:
        return name if name in AA_1_TO_3 else name
    if len(name) == 3:
        if name in protein_letters_3to1:
            return protein_letters_3to1[name]
        for k, v in AA_1_TO_3.items():
            if v == name:
                return k
    return name[0]


def get_row_target_aa(row: dict[str, Any]) -> str:
    """Extract canonical single-letter AA, prioritizing real_residue."""
    keys_to_check = [
        "real_residue",
        "retromutate_to",
        "canonical_residue",
        "recorded_residue",
        "observed_residue",
        "residue",
    ]
    for key in keys_to_check:
        val = row.get(key)
        if val is not None:
            norm = normalize_residue_name(val)
            if norm and norm != "X":
                return norm
    return ""


def load_residue_table(path: Optional[Path]) -> Optional[list[dict[str, Any]]]:
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


def preprocess_structure(input_path: Path, residue_rows: list[dict[str, Any]], output_pdb_path: Path):
    """
    Pre-processes input PDB/CIF with Biopython:
    1. Replaces non-standard/modified residues (e.g., SNM -> SER, MLY -> LYS) based on real_residue.
    2. Converts HETATM records for modified AAs into standard ATOM records.
    3. Strips non-protein ligands and water.
    """
    target_aa_map = {}
    for row in residue_rows:
        chain = str(row.get("chainID") or row.get("chain_id") or "A")
        seqnum = row.get("seqnum")
        if seqnum is None:
            continue
        try:
            seqnum = int(seqnum)
        except (TypeError, ValueError):
            continue

        target_1letter = get_row_target_aa(row)
        if target_1letter:
            target_3letter = AA_1_TO_3.get(target_1letter, "ALA")
            target_aa_map[(chain, seqnum)] = target_3letter

    if input_path.suffix.lower() in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(input_path.stem, str(input_path))

    class CleanStructureSelect(Select):
        def accept_residue(self, residue):
            hetflag, seqnum, icode = residue.id
            chain = residue.get_parent().id

            # Mutate modified residue to standard AA and convert HETATM -> ATOM
            if (chain, seqnum) in target_aa_map:
                residue.resname = target_aa_map[(chain, seqnum)]
                residue.id = (' ', seqnum, icode)
                return 1

            # Keep standard protein residues
            if hetflag == ' ' and residue.resname.strip().upper() in STANDARD_PROTEIN_3LETTER:
                return 1

            return 0

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path), CleanStructureSelect())


def parse_pdb_label(pdb_label: Optional[str], previous_chain: str, previous_resi: Optional[int]) -> tuple[str, int]:
    if not pdb_label:
        return previous_chain, (previous_resi + 1) if previous_resi is not None else 1

    parts = pdb_label.strip().split()
    if not parts:
        return previous_chain, (previous_resi + 1) if previous_resi is not None else 1

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
        resi = (previous_resi + 1) if previous_resi is not None else 1

    return chain, resi


def build_gap_list(pose, residue_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    gaps = []
    pose2pdb = pose.pdb_info().pose2pdb

    table_chains: dict[str, list[int]] = {}
    for row in residue_rows:
        ch = str(row.get("chainID") or row.get("chain_id") or "A")
        sq = row.get("seqnum")
        if sq is not None:
            try:
                table_chains.setdefault(ch, []).append(int(sq))
            except (TypeError, ValueError):
                pass

    chain_pose_residues: dict[str, list[tuple[int, int]]] = {}
    previous_chain = ""
    previous_resi = None

    for residue in pose.residues:
        seqpos = residue.seqpos()
        if residue.is_ligand() or residue.is_metal():
            continue

        pdb_label = pose2pdb(seqpos)
        chain, resi = parse_pdb_label(pdb_label, previous_chain, previous_resi)
        previous_chain = chain
        previous_resi = resi

        chain_pose_residues.setdefault(chain, []).append((seqpos, resi))

    for chain, pose_res_list in chain_pose_residues.items():
        if not pose_res_list:
            continue

        pose_res_list.sort(key=lambda x: x[0])
        first_pose_pos, first_resi = pose_res_list[0]
        last_pose_pos, last_resi = pose_res_list[-1]

        expected_seqnums = table_chains.get(chain, [])
        if expected_seqnums:
            min_table_seq = min(expected_seqnums)
            max_table_seq = max(expected_seqnums)

            # Detect missing N-terminal residues (e.g. seqnums -2, -1)
            if first_resi > min_table_seq:
                gaps.append({
                    'type': 'n_term',
                    'chain': chain,
                    'start': min_table_seq,
                    'end': first_resi - 1,
                    'anchor_pose_idx': first_pose_pos
                })

            # Detect missing C-terminal residues
            if last_resi < max_table_seq:
                gaps.append({
                    'type': 'c_term',
                    'chain': chain,
                    'start': last_resi + 1,
                    'end': max_table_seq,
                    'anchor_pose_idx': last_pose_pos
                })

        # Detect internal gaps
        for i in range(len(pose_res_list) - 1):
            curr_pose_pos, curr_resi = pose_res_list[i]
            next_pose_pos, next_resi = pose_res_list[i + 1]

            if next_resi > curr_resi + 1:
                gaps.append({
                    'type': 'internal',
                    'chain': chain,
                    'start': curr_resi + 1,
                    'end': next_resi - 1,
                    'anchor_pose_idx': curr_pose_pos
                })

    return gaps


def get_gap_sequence(rows: list[dict[str, Any]], chain: str, start: int, end: int) -> str:
    valid_rows = []
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
            valid_rows.append((seqnum, row))

    valid_rows.sort(key=lambda x: x[0])
    sequence = ""
    for _, row in valid_rows:
        sequence += get_row_target_aa(row)
    return sequence


def insert_missing_residues(pose, residue_rows: list[dict[str, Any]]):
    gaps = build_gap_list(pose, residue_rows)
    if not gaps:
        print("No missing residue gaps detected.")
        return pose

    cm = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    rts = cm.residue_type_set("fa_standard")
    rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
    rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue

    gaps.sort(key=lambda g: g['anchor_pose_idx'], reverse=True)

    for gap in gaps:
        gap_type = gap['type']
        chain = gap['chain']
        start = gap['start']
        end = gap['end']
        anchor_pos = gap['anchor_pose_idx']

        sequence = get_gap_sequence(residue_rows, chain, start, end)
        if not sequence:
            print(f"Missing sequence not found for gap {chain}:{start}-{end}")
            continue

        inserted_pose_indices = []

        if gap_type == 'n_term':
            insert_pos = anchor_pos
            for one_letter in sequence:
                residue_type = rts.get_representative_type_name1(one_letter)
                residue = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(residue_type)
                pose.prepend_polymer_residue_before_seqpos(residue, insert_pos, True)
                inserted_pose_indices.append(insert_pos)
                insert_pos += 1
        else:
            insert_pos = anchor_pos
            for one_letter in sequence:
                residue_type = rts.get_representative_type_name1(one_letter)
                residue = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(residue_type)
                pose.append_polymer_residue_after_seqpos(residue, insert_pos, True)
                insert_pos += 1
                inserted_pose_indices.append(insert_pos)

        # Reassign sequence numbers (including negative seqnums like -2, -1) to PDB header
        if pose.pdb_info():
            for idx, seqnum in enumerate(range(start, end + 1)):
                pose_pos = inserted_pose_indices[idx]
                pose.pdb_info().set_resinfo(res=pose_pos, chain_id=chain, pdb_res=seqnum)

        for new_pos in inserted_pose_indices:
            rm_lower(pose.conformation(), new_pos)
            rm_upper(pose.conformation(), new_pos)

        start_loop = max(1, min(inserted_pose_indices) - 1)
        end_loop = min(pose.total_residue(), max(inserted_pose_indices) + 1)
        cut_point = max(inserted_pose_indices)

        try:
            loop = pyrosetta.rosetta.protocols.loops.Loop(start_loop, end_loop, cut_point)
            loops = pyrosetta.rosetta.protocols.loops.Loops()
            loops.add_loop(loop)
            lm = pyrosetta.rosetta.protocols.loop_modeler.LoopModeler()
            lm.set_loops(loops)
            lm.enable_centroid_stage()
            lm.enable_fullatom_stage()
            lm.enable_build_stage()
            lm.apply(pose)
        except Exception as err:
            print(f"Loop modeling warning for {chain}:{start}-{end}: {err}")

    return pose


def fix_structure_with_pyrosetta(input_path: Path, residue_table_path: Path, output_pdb: Path):
    if not input_path.exists():
        raise FileNotFoundError(f"Input structure file not found: {input_path}")
    if not residue_table_path.exists():
        raise FileNotFoundError(f"Residue table file not found: {residue_table_path}")

    residue_rows = load_residue_table(residue_table_path)
    if residue_rows is None:
        raise ValueError(f"Unable to load residue table from {residue_table_path}")

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_handle:
        temp_pdb_path = Path(tmp_handle.name)

    try:
        preprocess_structure(input_path, residue_rows, temp_pdb_path)

        pyrosetta.init(extra_options="-mute all")
        pose = pyrosetta.pose_from_pdb(str(temp_pdb_path))

        repaired_pose = insert_missing_residues(pose, residue_rows)

        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        repaired_pose.dump_pdb(str(output_pdb))
        print(f"Saved cleaned structure to: {output_pdb}")
    finally:
        if temp_pdb_path.exists():
            temp_pdb_path.unlink(missing_ok=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fix a protein structure using PyRosetta and canonical residue mapping."
    )
    parser.add_argument("--input", "-i", required=True, help="Input PDB/mmCIF file path.")
    parser.add_argument("--residue-table", "-t", required=True, help="Path to <pdb>_canonical_residues.json.")
    parser.add_argument("--sequence", "-s", required=False, help="Path to canonical FASTA sequence file.")
    parser.add_argument("--output", "-o", required=True, help="Output cleaned PDB file path.")
    return parser.parse_args()

def main():
    args = parse_args()
    try:
        fix_structure_with_pyrosetta(Path(args.input), Path(args.residue_table), Path(args.output))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()