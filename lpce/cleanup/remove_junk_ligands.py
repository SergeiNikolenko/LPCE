import json
import multiprocessing
from collections import Counter
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm
import sys


def remove_junk_ligands_from_file(input_file_path: Path, junk_ligands: set) -> dict:
    """
    Removes junk ligands from a PDB file and writes the cleaned content back to the same file.

    Args:
        input_file_path (Path): The path to the PDB file to be processed.
        junk_ligands (set): A set of ligand IDs to be removed.

    Returns:
        dict: A dictionary with the count of removed ligands.
    """
    try:
        ligand_counts = Counter()
        temp_file_path = input_file_path.with_suffix(".tmp")
        changes_made = False

        with open(input_file_path) as f_in, open(temp_file_path, "w") as f_out:
            for line in f_in:
                if line.startswith("HETATM"):
                    ligand_id = line[17:20].strip()
                    if ligand_id in junk_ligands:
                        ligand_counts[ligand_id] += 1
                        changes_made = True
                        continue
                f_out.write(line)

        if changes_made:
            temp_file_path.replace(input_file_path)
        else:
            temp_file_path.unlink()

        return dict(ligand_counts)
    except Exception as e:
        logger.error(f"Failed to process {input_file_path}: {e}")
        return {"error": f"Error processing {input_file_path}: {e}"}


def remove_junk_ligands_from_directory(cfg) -> None:
    """
    Processes all PDB files in the specified directory, removing junk ligands from each file.

    Args:
        cfg (DictConfig): Configuration dictionary with paths and logging settings.

    Returns:
        None
    """
    input_directory = Path(cfg.paths.processed_dir)
    junk_ligands_file = Path(cfg.output_files.trash_ligands_json)
    summary_file_path = Path(cfg.output_files.removed_ligands_summary_json)
    log_file = cfg.logging.remove_junk_ligands_log_file

    # Setup logging to file
    logger.remove()
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )
    logger.info("========== Removing Junk Ligands ==========")
    # Load junk ligands from JSON file
    with open(junk_ligands_file) as file:
        junk_ligands = set(json.load(file))

    # Get all PDB files in the processed directory
    pdb_files = list(input_directory.rglob("*.pdb"))
    total_files = len(pdb_files)

    logger.info(f"Found {total_files} PDB files in {input_directory}")

    num_cores = multiprocessing.cpu_count() - 1

    # Process each file in parallel
    results = Parallel(n_jobs=num_cores)(
        delayed(remove_junk_ligands_from_file)(pdb_file, junk_ligands)
        for pdb_file in tqdm(
            pdb_files, desc="Removing junk ligands", unit="file", total=total_files
        )
    )

    # Summarize results
    total_ligand_counts = Counter()
    failed_files = 0

    for result in results:
        if isinstance(result, dict) and "error" not in result:
            total_ligand_counts.update(result)
        else:
            failed_files += 1

    successful_files = total_files - failed_files

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully processed: {successful_files}")
    logger.info(f"Failed to process: {failed_files}")
    logger.info(f"Total ligands removed: {sum(total_ligand_counts.values())}")

    # Save the summary to a JSON file
    with open(summary_file_path, "w") as summary_file:
        json.dump(dict(total_ligand_counts), summary_file, indent=4)

    logger.info(f"Summary of removed ligands saved to {summary_file_path}")
