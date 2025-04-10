import subprocess
import sys
from pathlib import Path

from loguru import logger
from tqdm import tqdm


def count_structures(directory: Path) -> int:
    """
    Counts the number of .ent.gz files in the specified directory.

    Args:
        directory (Path): Directory to search for .ent.gz files.

    Returns:
        int: Number of found .ent.gz files.
    """
    directory_path = Path(directory)
    return sum(1 for _ in directory_path.rglob("*.ent.gz"))


def extract_complexes(cfg) -> dict:
    """
    Synchronizes PDB structures from the RCSB PDB FTP server to the local directory specified by cfg.paths.raw_dir.

    Args:
        cfg: Configuration object containing paths and connection settings.

    Returns:
        dict: Dictionary with the number of new structures and total statistics.
    """
    raw_dir = Path(cfg.paths.raw_dir)
    rsync_port = cfg.rsync.port
    rsync_host = cfg.rsync.host

    raw_dir.mkdir(parents=True, exist_ok=True)
    logger.info("\n========== Extracting Complexes ==========")

    initial_count = count_structures(raw_dir)
    logger.info(f"Initial count of structures: {initial_count}")

    rsync_command = [
        "rsync",
        "-rlPt",
        "--delete",
        f"--port={rsync_port}",
        "--out-format='%n'",
        f"{rsync_host}::ftp_data/structures/divided/pdb/",
        str(raw_dir),
    ]

    try:
        # Dry-run to estimate total number of files
        dry_run_command = rsync_command + ["--dry-run"]
        logger.info("Running dry-run to estimate total files...")

        dry_run_process = subprocess.Popen(
            dry_run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        dry_run_output, dry_run_error = dry_run_process.communicate()

        if dry_run_process.returncode != 0:
            logger.error(f"Dry-run failed with error: {dry_run_error}")
            return {"new_structures": 0, "status": "error"}

        file_list = [line for line in dry_run_output.strip().split("\n") if line]
        estimated_total_files = len(file_list)
        logger.info(f"Estimated total files to sync: {estimated_total_files}")

        # Run the actual rsync process with progress bar
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
                final_count = count_structures(raw_dir)
                new_structures = final_count - initial_count
                logger.info(f"Complexes successfully extracted and saved to {raw_dir}")
                logger.info(f"Total structures: {final_count}")
                logger.info(f"New structures added: {new_structures}")
                return {
                    "new_structures": new_structures,
                    "total_structures": final_count,
                    "status": "success",
                }
            else:
                logger.error(
                    f"Rsync finished with errors, return code: {process.returncode}"
                )
                logger.error(process.stderr.read())
                return {"new_structures": 0, "status": "error"}

    except Exception as e:
        logger.error(f"Error during rsync: {e}")
        return {"new_structures": 0, "status": "error"}
