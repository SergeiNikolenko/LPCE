import subprocess
from pathlib import Path
from tqdm import tqdm
from loguru import logger


def remove_water_from_directory(input_dir: Path, log_file: str) -> None:
    """
    Processes all PDB files in the specified directory, removing water molecules from each file.

    This function uses a compiled C program (`remove_water`) to perform the water removal.
    Logs the process and statistics into the specified log file.

    Args:
        processed_dir (Path): The directory containing the PDB files to process.
        log_file (str): Path to the log file for logging the water removal process.

    Returns:
        None
    """
    # Add a separate log file for the water removal process
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )
    logger.info("========== Removing Water ==========")
    pdb_files = list(input_dir.rglob("*.pdb"))
    total_files = len(pdb_files)

    logger.info(f"Found {total_files} PDB files in {input_dir}")

    executable_path = "lpce/cleanup/remove_water"

    for pdb_file in tqdm(
        pdb_files, desc="Removing water", unit="file", total=total_files
    ):
        result = subprocess.run(
            [executable_path, str(pdb_file)], capture_output=True, text=True
        )

        if result.returncode == 0:
            pass
        else:
            logger.error(
                f"Failed to remove water from {pdb_file}. Error: {result.stderr}"
            )

    logger.info(f"Total structures processed: {total_files}")
