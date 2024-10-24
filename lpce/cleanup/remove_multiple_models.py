import os
import sys
from pathlib import Path

import joblib
from loguru import logger
from tqdm import tqdm


def count_models_in_file(file_path: Path) -> int:
    """
    Counts the number of models in a PDB file.

    Args:
        file_path (Path): Path to the PDB file.

    Returns:
        int: Number of models in the PDB file.
    """
    models = 0
    try:
        with open(file_path) as file:
            for line in file:
                if line.startswith("MODEL"):
                    models += 1
        return models
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return 0


def process_file(file: Path) -> bool:
    """
    Processes a single PDB file and removes it if it contains multiple models.

    Args:
        file (Path): The path to the PDB file.

    Returns:
        bool: True if the file is retained, False if it is removed.
    """
    models = count_models_in_file(file)
    if models > 1:
        try:
            os.remove(file)
            # logger.debug(f"Removed file {file} containing {models} models.")
            return False
        except Exception as e:
            logger.error(f"Error removing file {file}: {e}")
            return False
    return True


def remove_multiple_models_from_directory(cfg: object) -> None:
    """
    Removes PDB files containing multiple models from the specified directory.

    Args:
        input_dir (Path): The directory to scan for PDB files.
        log_file (str): The path to the log file where the process is logged.

    Returns:
        None
    """
    input_dir = Path(cfg.paths.processed_dir)

    logger.info("========== Removing PDB files with multiple models ==========")
    files = list(input_dir.glob("*.pdb"))
    total_files = len(files)

    logger.info(f"Total PDB files to analyze: {total_files}")

    results = joblib.Parallel(n_jobs=-1)(
        joblib.delayed(process_file)(file)
        for file in tqdm(files, desc="Removing multiple models")
    )
    retained_files = sum(results)
    removed_files = total_files - retained_files
    remaining_percentage = (
        (retained_files / total_files * 100) if total_files > 0 else 0
    )

    logger.info(f"Total files analyzed: {total_files:,}")
    logger.info(f"Files retained after removal: {retained_files:,}")
    logger.info(f"Files removed: {removed_files:,}")
    logger.info(f"Percentage of files retained: {remaining_percentage:.2f}%")
