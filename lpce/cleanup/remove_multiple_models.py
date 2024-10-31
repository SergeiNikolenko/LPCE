import os
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


def process_file(file: Path) -> tuple[bool, str]:
    """
    Processes a single PDB file and removes it if it contains multiple models.

    Args:
        file (Path): The path to the PDB file.

    Returns:
        tuple: (bool, str) where bool indicates if the file is retained, 
               and str is the filename (used for tracking removed files).
    """
    models = count_models_in_file(file)
    if models > 1:
        try:
            os.remove(file)
            return False, file.stem[3:]  # Extract file name without "pdb" prefix and ".pdb" suffix
        except Exception as e:
            logger.error(f"Error removing file {file}: {e}")
            return False, file.stem[3:]
    return True, file.stem[3:]


def remove_multiple_models_from_directory(cfg: object) -> dict:
    """
    Removes PDB files containing multiple models from the specified directory.

    Args:
        cfg (object): Configuration object with paths.

    Returns:
        dict: {"removed_files": list} - list of removed file names.
    """
    input_dir = Path(cfg.paths.processed_dir)

    logger.info("========== Removing PDB files with multiple models ==========")
    files = list(input_dir.glob("*.pdb"))
    total_files = len(files)

    logger.info(f"Total PDB files to analyze: {total_files}")

    results = joblib.Parallel(n_jobs=cfg.n_jobs)(
        joblib.delayed(process_file)(file)
        for file in tqdm(files, desc="Removing multiple models")
    )

    removed_files = [file for retained, file in results if not retained]
    retained_files = total_files - len(removed_files)
    remaining_percentage = (retained_files / total_files * 100) if total_files > 0 else 0

    logger.info(f"Total files analyzed: {total_files:,}")
    logger.info(f"Files retained after removal: {retained_files:,}")
    logger.info(f"Files removed: {len(removed_files):,}")
    logger.info(f"Percentage of files retained: {remaining_percentage:.2f}%")

    return {"removed_files": removed_files}
