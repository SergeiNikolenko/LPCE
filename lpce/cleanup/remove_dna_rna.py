import os
from pathlib import Path
import joblib
from tqdm import tqdm
from loguru import logger


def contains_dna_rna_sequence(content: str) -> bool:
    """
    Checks if the content of a PDB file contains DNA or RNA sequences in SEQRES or ATOM lines.

    Args:
        content (str): The content of the PDB file.

    Returns:
        bool: True if the file contains DNA or RNA sequences, False otherwise.
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

    # Check ATOM lines for specific nucleotide atoms
    for line in atom_lines:
        if line[17:20].strip() in nucleotides:
            return True

    return False


def process_file(file: Path) -> bool:
    """
    Processes a single PDB file and removes it if it contains DNA or RNA sequences.

    Args:
        file (Path): The path to the PDB file.

    Returns:
        bool: True if the file is retained, False if it is removed.
    """
    try:
        with open(file) as f:
            content = f.read()
            if contains_dna_rna_sequence(content):
                os.remove(file)
                # logger.info(f"Removed file {file} containing DNA/RNA sequences.")
                return False
        return True
    except Exception as e:
        logger.error(f"Error processing file {file}: {e}")
        return False


def remove_dna_rna_from_directory(input_dir: Path, log_file: str) -> None:
    """
    Removes PDB files containing DNA or RNA sequences from the specified directory.

    Args:
        input_dir (Path): The directory to scan for PDB files.
        log_file (str): The path to the log file where the process is logged.

    Returns:
        None
    """
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )
    logger.info("========== Removing DNA/RNA ==========")
    files = list(input_dir.glob("*.pdb"))
    total_files = len(files)

    logger.info(f"Total PDB files to analyze: {total_files}")

    results = joblib.Parallel(n_jobs=-1)(
        joblib.delayed(process_file)(file)
        for file in tqdm(files, desc="Removing DNA/RNA")
    )
    retained_files = sum(results)

    remaining_percentage = (retained_files / total_files) * 100

    logger.info(f"Total files analyzed: {total_files:,}")
    logger.info(f"Files retained after removal: {retained_files:,}")
    logger.info(f"Percentage of files retained: {remaining_percentage:.2f}%")
