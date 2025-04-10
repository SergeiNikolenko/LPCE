import os
import re
from pathlib import Path

from Bio import Align
from Bio.PDB import PDBParser, Superimposer
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def get_pdb_files(input_dir):
    """
    Finds all PDB files in the input directory and subdirectories.
    """
    pdb_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))
    return pdb_files


def parse_biomolecule_chains(pdb_lines):
    """
    Parses REMARK 350 sections to extract biomolecule chain assignments.
    """
    biomolecules = {}
    current_biomolecule = None
    for line in pdb_lines:
        if line.startswith("REMARK 350 BIOMOLECULE:"):
            current_biomolecule = int(line.strip().split(":")[1])
            biomolecules[current_biomolecule] = []
        elif line.startswith("REMARK 350 APPLY THE FOLLOWING TO CHAINS:"):
            chains = re.findall(r"[A-Za-z0-9]+", line.split(":")[1])
            biomolecules[current_biomolecule].extend(chains)
    return biomolecules


def get_header_footer(pdb_lines):
    """
    Extracts the header and footer sections from the PDB file.
    """
    header_lines = []
    footer_lines = []
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            break
        header_lines.append(line)
    for line in pdb_lines[::-1]:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            break
        footer_lines.insert(0, line)
    return header_lines, footer_lines


def extract_chains(pdb_lines, chains):
    """
    Extracts the ATOM/HETATM lines for the specified chains.
    """
    biomol_atom_lines = []
    atom_section_started = False
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_section_started = True
            chain_id = line[21]
            if chain_id in chains:
                biomol_atom_lines.append(line)
        elif atom_section_started and line.startswith("TER"):
            chain_id = line[21]
            if chain_id in chains:
                biomol_atom_lines.append(line)
        elif atom_section_started and line.startswith("ENDMDL"):
            biomol_atom_lines.append(line)
            break
    return biomol_atom_lines


def split_biomolecule_pdb(pdb_file, output_dir):
    """
    Splits a PDB file into separate files for each biomolecule based on REMARK 350.
    """
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    with open(pdb_file) as file:
        pdb_lines = file.readlines()
    biomolecules = parse_biomolecule_chains(pdb_lines)
    header_lines, footer_lines = get_header_footer(pdb_lines)
    biomolecule_files = []
    for biomol_id, chains in biomolecules.items():
        biomol_atom_lines = extract_chains(pdb_lines, chains)
        biomol_file = os.path.join(output_dir, f"{pdb_name}_bioml_{biomol_id}.pdb")
        with open(biomol_file, "w") as file:
            file.writelines(header_lines + biomol_atom_lines)
        biomolecule_files.append(biomol_file)
    return biomolecule_files


def get_chain_sequences(structure):
    """
    Extracts the sequences of amino acids for each chain in a structure.
    """
    sequences = {}
    for model in structure:
        for chain in model:
            seq = "".join(
                seq1(residue.resname) for residue in chain if residue.id[0] == " "
            )
            sequences[chain.id] = Seq(seq)
    return sequences


def compute_identity(seq1, seq2):
    """
    Computes the sequence identity percentage between two sequences using Bio.Align.PairwiseAligner.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    aligned_seq1 = alignment.aligned[0]
    aligned_seq2 = alignment.aligned[1]
    matches = 0
    total = 0
    for (start1, end1), (start2, end2) in zip(aligned_seq1, aligned_seq2):
        length = end1 - start1
        matches += length
        total += length
    identity = (matches / total) * 100 if total > 0 else 0
    return identity


def are_structures_identical(structure1, structure2, cutoff_identity=95.0):
    """
    Checks if two structures are identical based on sequence identity and RMSD.
    """
    chains1 = list(structure1.get_chains())
    chains2 = list(structure2.get_chains())
    if len(chains1) != len(chains2):
        return False
    sequences1 = get_chain_sequences(structure1)
    sequences2 = get_chain_sequences(structure2)
    for chain_id1, seq_chain1 in sequences1.items():
        found = False
        for chain_id2, seq_chain2 in sequences2.items():
            identity = compute_identity(seq_chain1, seq_chain2)
            if identity >= cutoff_identity:
                ca_atoms1 = [
                    atom
                    for atom in structure1[0][chain_id1].get_atoms()
                    if atom.get_id() == "CA"
                ]
                ca_atoms2 = [
                    atom
                    for atom in structure2[0][chain_id2].get_atoms()
                    if atom.get_id() == "CA"
                ]
                if len(ca_atoms1) != len(ca_atoms2):
                    continue
                sup = Superimposer()
                sup.set_atoms(ca_atoms1, ca_atoms2)
                rms = sup.rms
                if rms < 0.5:
                    found = True
                    break
        if not found:
            return False
    return True


def check_inclusion(seq1, seq2):
    """
    Checks if seq1 includes seq2 or seq2 includes seq1.
    """
    return str(seq1).find(str(seq2)) != -1 or str(seq2).find(str(seq1)) != -1


def remove_duplicate_and_included_biomolecules(biomolecule_files):
    """
    Removes duplicate biomolecules and those where one includes the other.
    """
    unique_files = []
    duplicates = []
    inclusions = []
    parser = PDBParser(QUIET=True)
    for i, file1 in enumerate(biomolecule_files):
        structure1 = parser.get_structure(f"biomol_{i+1}", file1)
        is_duplicate_or_included = False
        for j, file2 in enumerate(unique_files):
            structure2 = parser.get_structure(f"unique_{j+1}", file2)
            if are_structures_identical(structure1, structure2):
                duplicates.append(file1)
                is_duplicate_or_included = True
                break
            sequences1 = get_chain_sequences(structure1)
            sequences2 = get_chain_sequences(structure2)
            for chain_id1, seq_chain1 in sequences1.items():
                for chain_id2, seq_chain2 in sequences2.items():
                    if check_inclusion(seq_chain1, seq_chain2):
                        inclusions.append(file1)
                        is_duplicate_or_included = True
                        break
                if is_duplicate_or_included:
                    break
        if not is_duplicate_or_included:
            unique_files.append(file1)
    for dup_file in duplicates + inclusions:
        os.remove(dup_file)
    return unique_files, duplicates, inclusions


def process_pdb_file_with_inclusion_check(pdb_file, output_dir):
    """
    Processes a single PDB file by splitting it into biomolecules,
    checking for duplicates and inclusions, and removing any duplicates or included biomolecules.
    """
    try:
        biomolecule_files = split_biomolecule_pdb(pdb_file, output_dir)
        unique_files, duplicates, inclusions = (
            remove_duplicate_and_included_biomolecules(biomolecule_files)
        )
        return unique_files, duplicates, inclusions
    except Exception:
        # logger.error(f"Error processing {pdb_file}: {e}")
        return [], [], []


def bioml_split(cfg) -> dict:
    """
    Splits PDB files into biomolecules based on REMARK 350 information,
    removes duplicate biomolecules, and saves them into output_directory.
    Returns a dictionary containing lists of created, removed (duplicates), and included files.

    Args:
        cfg (DictConfig): Configuration object containing paths and logging settings.

    Returns:
        dict: Contains lists of unique files (kept), duplicates removed, and inclusions.
    """
    input_dir = Path(cfg.paths.processed_dir)
    output_dir = Path(cfg.paths.bioml_dir)

    logger.info("\n========== Split Bioml ==========")

    pdb_files = get_pdb_files(input_dir)
    total_files = len(pdb_files)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process files in parallel
    all_unique_files = []
    all_duplicates = []
    all_inclusions = []
    total_biomolecules_created = 0

    with Parallel(n_jobs=cfg.n_jobs) as parallel:
        results = parallel(
            delayed(process_pdb_file_with_inclusion_check)(pdb_file, output_dir)
            for pdb_file in tqdm(pdb_files, desc="Processing PDB files")
        )

    # Collect results
    for unique_files, duplicates, inclusions in results:
        total_biomolecules_created += (
            len(unique_files) + len(duplicates) + len(inclusions)
        )
        all_unique_files.extend(unique_files)
        all_duplicates.extend(duplicates)
        all_inclusions.extend(inclusions)

    logger.info(f"Total PDB files: {total_files}")
    logger.info(f"Total biomolecules created: {total_biomolecules_created:,}")
    logger.info(f"Total duplicates removed: {len(all_duplicates):,}")
    logger.info(f"Total inclusions found: {len(all_inclusions):,}")
    logger.info(f"Biomolecules remaining after filtering: {len(all_unique_files):,}")

    return {
        "unique_files": all_unique_files,
        "duplicates_removed": all_duplicates,
        "inclusions": all_inclusions,
    }
