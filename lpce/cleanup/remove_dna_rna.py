import os
from pathlib import Path

import joblib
from loguru import logger
from tqdm import tqdm


def contains_dna_rna_sequence(content: str) -> bool:
    """
    Checks if the content of a PDB file contains DNA or RNA sequences.

    Args:
        content (str): The content of the PDB file.

    Returns:
        bool: True if the file contains DNA/RNA sequences, otherwise False.
    """
    nucleotides = {
        "A",
        "T",
        "G",
        "C",
        "U",
        "DA",
        "DT",
        "DG",
        "DC",
        "DU",
        "RA",
        "RT",
        "RG",
        "RC",
        "RU",
    }
    seqres_lines = [line for line in content.splitlines() if line.startswith("SEQRES")]
    atom_lines = [line for line in content.splitlines() if line.startswith("ATOM")]

    # Check SEQRES lines
    for line in seqres_lines:
        sequence = line[19:].split()
        if any(n in nucleotides for n in sequence):
            return True

    # Check ATOM lines for specific nucleotides
    for line in atom_lines:
        if line[17:22].strip() in nucleotides:
            return True

    return False


def process_file(file: Path) -> str:
    """
    Processes a single PDB file and removes it if it contains DNA or RNA sequences.

    Args:
        file (Path): The path to the PDB file.

    Returns:
        str: 'removed' if the file was removed, or 'retained' if the file was kept.
    """
    try:
        with open(file) as f:
            content = f.read()
            if contains_dna_rna_sequence(content):
                os.remove(file)
                return "removed"
        return "retained"
    except Exception:
        # If an error occurs, the file is not removed
        return "error"


def remove_dna_rna_from_directory(cfg) -> dict:
    """
    Removes PDB files containing DNA or RNA sequences from the specified directory.

    Args:
        cfg (object): Configuration containing directory paths.

    Returns:
        dict: A dictionary with results: lists of removed and retained files.
    """
    input_dir = Path(cfg.paths.processed_dir)
    files = list(input_dir.glob("*.pdb"))

    # Use parallel processing for files
    results = joblib.Parallel(n_jobs=cfg.n_jobs)(
        joblib.delayed(process_file)(file)
        for file in tqdm(files, desc="Removing DNA/RNA")
    )

    # Form lists of removed and retained files
    retained_files = [
        files[i].stem[3:] if files[i].stem[:3] == "pdb" else files[i].stem
        for i, result in enumerate(results)
        if result == "retained"
    ]
    removed_files = [
        files[i].stem[3:] if files[i].stem[:3] == "pdb" else files[i].stem
        for i, result in enumerate(results)
        if result == "removed"
    ]
    errors = [
        files[i].stem[3:] if files[i].stem[:3] == "pdb" else files[i].stem
        for i, result in enumerate(results)
        if result == "error"
    ]

    # Return results
    logger.info(f"\n======== Removing DNA/RNA from PDB files ========")
    logger.info(f"Total files analyzed: {len(files):,}")
    logger.info(f"Files removed: {len(removed_files):,}")
    logger.info(f"Files retained: {len(retained_files):,}")
    logger.info(
        f"Percentage of files retained: {len(retained_files) / len(files) * 100:.2f}%"
    )

    return {"removed_files": removed_files}
