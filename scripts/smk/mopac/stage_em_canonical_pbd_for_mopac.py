def convert_pdb_to_mop(pdb_path, mop_path):
    """
    Converts a PDB structure to a MOPAC .mop input file.
    
    Header structure required by MOPAC:
      Line 1: Keywords
      Line 2: Title / Descriptor
      Line 3: Comments / Blanks
      Line 4+: Atomic coordinates (Symbol X 0 Y 0 Z 0)
    """
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Input PDB file not found: {pdb_path}")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(mop_path)), exist_ok=True)

    # 1. Prepare Header Block
    header_keywords = "PM7 1SCF EPS=78.3 MOZYME\n"
    title_line = f"Canonical EM baseline for {os.path.basename(pdb_path)}\n"
    comment_line = "Generated automatically via Snakemake stage_em_pdb_for_mopac\n"

    lines = [header_keywords, title_line, comment_line]
    atom_count = 0

    # 2. Extract Atom Coordinates
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # Extract Element Symbol (PDB columns 77-78)
                element = line[76:78].strip()
                
                # Fallback to atom name if element column is blank
                if not element:
                    atom_name = line[12:16].strip()
                    element = atom_name[0]  # Take first character (e.g., C, N, O, H, S)

                # Extract Coordinates (PDB columns 31-54)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Format for MOPAC: Element X 0 Y 0 Z 0 (0s freeze coords for 1SCF)
                mop_line = f"{element:<2} {x:12.6f} 0 {y:12.6f} 0 {z:12.6f} 0\n"
                lines.append(mop_line)
                atom_count += 1

    if atom_count == 0:
        raise ValueError(f"No valid ATOM or HETATM records found in PDB file: {pdb_path}")

    # 3. Write Output .mop File
    with open(mop_path, "w") as f:
        f.writelines(lines)

    print(f"✅ Successfully staged MOPAC input file with {atom_count} atoms: {mop_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python stage_em_canonical_pdb_for_mopac.py <input_pdb> <output_mop>")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_mop = sys.argv[2]

    try:
        convert_pdb_to_mop(input_pdb, output_mop)
    except Exception as e:
        print(f"❌ Error staging MOPAC input: {str(e)}", file=sys.stderr)
        sys.exit(1)