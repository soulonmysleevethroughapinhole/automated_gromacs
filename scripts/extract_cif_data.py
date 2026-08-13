#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from typing import Dict, List, Any

from Bio.PDB.MMCIFParser import MMCIF2Dict
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract metadata and structural information from an mmCIF protein structure file."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the input CIF file.",
    )
    parser.add_argument(
        "--metadata",
        "-m",
        required=True,
        help="Path to the output JSON file for metadata.",
    )
    return parser.parse_args()


def _to_list(value):
    """Convert scalar values to lists for consistent handling."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def extract_metadata(cif_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Extract metadata from the CIF dictionary."""
    output_dict = {}

    # Basic identification
    output_dict["pdb"] = cif_dict.get("_entry.id", ["UNKNOWN"])[0]
    output_dict["title"] = cif_dict.get("_struct.title", ["No Title"])[0]

    # Crystallization and measurement data
    output_dict["xtal_pH"] = cif_dict.get("_exptl_crystal_grow.pH")
    output_dict["xtal_temp"] = cif_dict.get("_exptl_crystal_grow.temp")
    output_dict["src_tissue"] = cif_dict.get("_entity_src_gen.gene_src_tissue")
    output_dict["src_cellular"] = cif_dict.get(
        "_entity_src_gen.pdbx_gene_src_cellular_location"
    )

    # Journal/Publication data
    output_dict["jrnl_author"] = cif_dict.get("_citation_author.name", None)
    output_dict["jrnl_title"] = cif_dict.get("_citation.title", [None])[0]
    output_dict["jrnl_ref"] = cif_dict.get("_citation.journal_abbrev", [None])[0]
    output_dict["jrnl_doi"] = cif_dict.get("_citation.pdbx_database_id_DOI", [None])[0]

    # Mutation/Engineering
    entity_details = "".join(cif_dict.get("_entity.details", []))
    output_dict["mutation"] = "MUTANT" in entity_details.upper()

    # Resolution
    res_val = cif_dict.get(
        "_reflns.d_resolution_high",
        cif_dict.get("_refine.ls_d_res_high", ["N/A"]),
    )
    output_dict["resolution"] = res_val[0] if isinstance(res_val, list) else res_val

    # Sequence (if available)
    output_dict["struct_sequence"] = [entry.replace('\n', "") for entry in cif_dict.get("_struct_ref.pdbx_seq_one_letter_code")]
    output_dict["sequence_can"] = [entry.replace('\n', "") for entry in cif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can")]
    output_dict["sequence"] = [entry.replace('\n', "") for entry in cif_dict.get("_entity_poly.pdbx_seq_one_letter_code")]

    return output_dict


def extract_missing_residues(cif_dict: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract missing/unobserved residues from mmCIF 'unobs' table."""
    try:
        res_names = _to_list(
            cif_dict.get("_pdbx_unobs_or_zero_occ_residues.label_comp_id", [])
        )
        chain_ids = _to_list(
            cif_dict.get("_pdbx_unobs_or_zero_occ_residues.auth_asym_id", [])
        )
        seq_nums = _to_list(
            cif_dict.get("_pdbx_unobs_or_zero_occ_residues.auth_seq_id", [])
        )

        missing_list = []
        for res, chain, seq in zip(res_names, chain_ids, seq_nums):
            missing_list.append(
                {"residue_name": res, "chain": chain, "residue_number": seq}
            )

        return missing_list if missing_list else []

    except Exception as e:
        print(f"Error extracting missing residues: {e}")
        return []


def extract_modified_residues(cif_dict: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract modified residues from mmCIF 'mod_residue' table."""
    try:
        mod_names = _to_list(
            cif_dict.get("_pdbx_struct_mod_residue.label_comp_id", [])
        )
        mod_chains = _to_list(
            cif_dict.get("_pdbx_struct_mod_residue.auth_asym_id", [])
        )
        mod_seqs = _to_list(
            cif_dict.get("_pdbx_struct_mod_residue.auth_seq_id", [])
        )
        parent_res = _to_list(
            cif_dict.get("_pdbx_struct_mod_residue.parent_comp_id", [])
        )

        modified_list = []
        for i in range(len(mod_names)):
            modified_list.append(
                {
                    "residue_name": mod_names[i],
                    "chain": mod_chains[i],
                    "residue_number": mod_seqs[i],
                    "parent_residue": parent_res[i] if i < len(parent_res) else "UNK",
                }
            )

        return modified_list if modified_list else []

    except Exception as e:
        print(f"Error extracting modified residues: {e}")
        return []


def extract_sequences_and_hetatms(structure) -> Dict[str, Any]:
    """Extract sequences and hetatom information per chain."""
    chain_data = {}

    for chain in structure.get_chains():
        chain_id = chain.get_id()
        seq = []
        hetatms = []

        for residue in chain:
            res_name = residue.get_resname().strip()

            # Try to convert 3-letter code to 1-letter code
            amino_acid = protein_letters_3to1.get(res_name.capitalize(), "X")

            if is_aa(residue):
                seq.append(amino_acid)
            else:
                hetatms.append(res_name)

        chain_data[f"chain_{chain_id}"] = "".join(seq)
        if hetatms:
            chain_data[f"hetatms_ch_{chain_id}"] = hetatms

    return chain_data


def parse_cif_file(cif_path: Path) -> tuple[Dict[str, Any], Any]:
    """Parse CIF file and extract structure information."""

    mmcif_dict = MMCIF2Dict(str(cif_path))
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", str(cif_path))

    return mmcif_dict, structure


def write_metadata(
    metadata: Dict[str, Any],
    missing_residues: List[Dict[str, Any]],
    modified_residues: List[Dict[str, Any]],
    chain_data: Dict[str, Any],
    output_path: Path,
):
    """Write metadata to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    full_metadata = {
        **metadata,
        "structural_issues": {
            "missing_residues": missing_residues if missing_residues else None,
            "modified_residues": modified_residues if modified_residues else None,
        },
        **chain_data,
    }
    with open(output_path, "w") as f:
        json.dump(full_metadata, f, indent=2)


def main():
    args = parse_args()
    input_path = Path(args.input)
    metadata_path = Path(args.metadata)

    if not input_path.exists():
        raise SystemExit(f"Input CIF file not found: {input_path}")

    # Parse the CIF file
    cif_dict, structure = parse_cif_file(input_path)

    # Extract all components
    metadata = extract_metadata(cif_dict)
    missing_residues = extract_missing_residues(cif_dict)
    modified_residues = extract_modified_residues(cif_dict)
    chain_data = extract_sequences_and_hetatms(structure)

    # Write outputs
    write_metadata(
        metadata,
        missing_residues,
        modified_residues,
        chain_data,
        metadata_path,
    )

    print(f"Metadata written to: {metadata_path}")
    if missing_residues:
        print(f"Missing residues: {len(missing_residues)} found")
    if modified_residues:
        print(f"Modified residues: {len(modified_residues)} found")


if __name__ == "__main__":
    main()
