import subprocess
from pathlib import Path

from loguru import logger
from tqdm import tqdm


def remove_water_from_directory(cfg: object) -> dict:
    """
    Processes all PDB files in the specified directory, removing water molecules from each file.

    This function uses a compiled C program (`remove_water`) to perform the water removal.
    Logs the process and statistics using the existing logger.

    Args:
        input_dir (Path): The directory containing the PDB files to process.

    Returns:
        dict: Dictionary with lists of successfully processed and failed files.
    """
    logger.info("\n========== Removing Water from PDB files ==========")
    pdb_files = list(Path(cfg.paths.processed_dir).rglob("*.pdb"))
    total_files = len(pdb_files)

    logger.info(f"Found {total_files} PDB files in {Path(cfg.paths.processed_dir)}")

    executable_path = "lpce/cleanup/remove_water"

    processed_files = []
    failed_files = []

    for pdb_file in tqdm(
        pdb_files, desc="Removing water", unit="file", total=total_files
    ):
        result = subprocess.run(
            [executable_path, str(pdb_file)],
            capture_output=True,
            text=True,
            shell=False,
        )

        if result.returncode == 0:
            processed_files.append(str(pdb_file))
        else:
            logger.error(
                f"Failed to remove water from {pdb_file}. Error: {result.stderr}"
            )
            failed_files.append(str(pdb_file))

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully processed: {len(processed_files)}")
    logger.info(f"Failed to process: {len(failed_files)}")

    # Returning the results as a dictionary
    return {"processed": processed_files, "failed": failed_files}
