import subprocess
import sys
from pathlib import Path
from tqdm import tqdm
from loguru import logger


def count_structures(directory: Path) -> int:
    """
    Counts the number of .ent.gz files in the specified directory.

    Args:
        directory (Path): The directory to search for .ent.gz files.

    Returns:
        int: The number of .ent.gz files found.
    """
    directory_path = Path(directory)
    return sum(1 for _ in directory_path.rglob("*.ent.gz"))


def extract_complexes(raw_dir: Path, rsync_port: int = 33444, rsync_host: str = "rsync.rcsb.org") -> int:
    """
    Synchronizes PDB structures from the RCSB PDB FTP server to the local directory specified by RAW_DIR.

    This function performs a dry-run to estimate the total number of files to be synced, then proceeds with the actual
    synchronization using rsync. The function prints progress and statistics, including the number of new structures
    added.

    Returns:
        int: The number of new structures added during the sync process. Returns 0 if an error occurred.
    """
    output_path = Path(raw_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    initial_count = count_structures(output_path)
    logger.info(f"Initial count of structures: {initial_count}")

    rsync_command = [
        "rsync",
        "-rlPt",
        "--delete",
        f"--port={rsync_port}",
        f"{rsync_host}::ftp_data/structures/divided/pdb/",
        str(output_path),
    ]

    try:
        dry_run_command = rsync_command + ["--dry-run"]
        logger.info("Running dry-run to estimate total files...")

        dry_run_process = subprocess.Popen(
            dry_run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        dry_run_output, dry_run_error = dry_run_process.communicate()

        if dry_run_process.returncode != 0:
            logger.error(f"Dry-run failed with error: {dry_run_error}")
            return 0

        estimated_total_files = len(dry_run_output.strip().split("\n"))
        logger.info(f"Estimated total files to sync: {estimated_total_files}")

        with tqdm(
            total=estimated_total_files,
            desc="Syncing files",
            unit="file",
            file=sys.stdout,
        ) as pbar:
            process = subprocess.Popen(
                rsync_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )

            for _ in process.stdout:
                pbar.update(1)

            process.wait()

            if process.returncode == 0:
                final_count = count_structures(output_path)
                new_structures = final_count - initial_count
                logger.info(f"Complexes successfully extracted and saved to {output_path}")
                logger.info(f"Total structures: {final_count}")
                logger.info(f"New structures added: {new_structures}")
                return new_structures
            else:
                logger.error(f"rsync finished with errors, return code: {process.returncode}")
                logger.error(process.stderr.read())
                return 0

    except Exception as e:
        logger.error(f"Error during rsync: {e}")
        raise e
