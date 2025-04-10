import json
import subprocess
import sys
import tempfile
from pathlib import Path

from loguru import logger
from tqdm import tqdm

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


def run_foldseek(
    cfg, input_dir, tmscore_threshold=0.5, fident_threshold=0.9, n_jobs=56
):
    """
    Runs Foldseek with specified parameters.

    :param cfg: Configuration object with paths and settings.
    :param input_dir: Directory containing input PDB files.
    :param tmscore_threshold: TM-score threshold for Foldseek.
    :param fident_threshold: Identity threshold for grouping results.
    :param n_jobs: Number of threads for Foldseek.
    """
    identical_pairs = {}

    with tempfile.TemporaryDirectory() as tmp_dir:
        output_dir = Path(tmp_dir)
        db_path = output_dir / "foldseek_db"
        aln_output = output_dir / "alignment_results.m8"
        foldseek_log = cfg.foldseek.foldseek_log
        foldseek_bin = cfg.foldseek.foldseek_path
        prostt5_path = cfg.foldseek.prostt5_path

        logger.info("Starting run_foldseek")
        logger.info(f"Input directory: {input_dir}")
        logger.info(f"Foldseek binary: {foldseek_bin}")

        try:
            # Handle the 'createdb' step
            with open(foldseek_log, "w") as log_file:
                logger.info("Creating Foldseek database...")
                createdb_command = [
                    foldseek_bin,
                    "createdb",
                    str(input_dir),
                    str(db_path),
                    "--threads",
                    str(n_jobs),
                ]
                process = subprocess.Popen(
                    createdb_command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )
                logger.info("createdb process started")

                # Initialize progress bar for 'createdb'
                total_equals_createdb = (
                    65  # Assuming progress bar has 65 '=' characters
                )
                pbar_createdb = tqdm(
                    total=total_equals_createdb, desc="Creating Database", unit="chars"
                )

                for line in process.stdout:
                    log_file.write(line)
                    if line.startswith("[") and "]" in line:
                        # Extract the progress bar part
                        progress_bar = line[line.find("[") : line.find("]") + 1]
                        num_equals = progress_bar.count("=")
                        pbar_createdb.n = num_equals
                        pbar_createdb.refresh()
                process.wait()
                pbar_createdb.close()
                if process.returncode != 0:
                    logger.error("Foldseek exited with an error during 'createdb'")
                    raise RuntimeError("Foldseek 'createdb' exited with an error.")
                logger.info(f"Database created at {db_path}")

            # Handle the 'easy-search' step
            with open(foldseek_log, "a") as log_file:
                logger.info("Starting easy-search...")
                command = [
                    foldseek_bin,
                    "easy-search",
                    str(db_path),
                    str(db_path),
                    str(aln_output),
                    str(output_dir),
                    "--format-output",
                    "query,target,fident,evalue",
                    "--max-seqs",
                    "10000",
                    "--alignment-type",
                    "2",
                    "--tmscore-threshold",
                    str(tmscore_threshold),
                    "--threads",
                    str(n_jobs),
                    "--prostt5-model",
                    prostt5_path,
                ]

                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )

                logger.info("easy-search process started")

                # Initialize progress bar for 'easy-search'
                total_equals_easysearch = 5 * 65
                pbar_easysearch = tqdm(
                    total=total_equals_easysearch,
                    desc="Foldseek Progress",
                    unit="chars",
                )

                # Read the output character by character
                while True:
                    char = process.stdout.read(1)
                    if not char:
                        break
                    log_file.write(char)
                    if char == "=":
                        pbar_easysearch.update(1)

                process.wait()
                pbar_easysearch.close()

                if process.returncode != 0:
                    logger.error("Foldseek exited with an error in easy-search")
                    raise RuntimeError("Foldseek exited with an error.")

        except subprocess.CalledProcessError as e:
            logger.error(f"Foldseek failed to create the database: {e}")
            sys.exit(1)

        # Reading and processing Foldseek results
        logger.info("Reading Foldseek alignment results")
        with open(aln_output) as f:
            for line in f:
                query, target, fident, evalue = line.strip().split()[:4]
                fident = float(fident)
                if fident >= fident_threshold:
                    if query not in identical_pairs:
                        identical_pairs[query] = {query}
                    identical_pairs[query].add(target)

        logger.info("Grouping identical structures")
        groups = []
        visited = set()
        for key, group in identical_pairs.items():
            if key not in visited:
                groups.append(group)
                visited.update(group)

        logger.info(
            f"run_foldseek completed, found {len(groups)} groups of identical structures"
        )
        return groups


def find_duplicates_foldseek(cfg):
    """
    Find duplicates in the input directory using Foldseek.
    """
    input_dir = Path(cfg.paths.separated_dir)
    tmscore_threshold = cfg.foldseek.tmscore_threshold
    fident_threshold = cfg.foldseek.fident_threshold
    n_jobs = cfg.foldseek.n_jobs
    logger.info("\n========== Foldseek analysis completed ==========")

    pdb_files = list(input_dir.glob("*.pdb"))
    total_input_files = len(pdb_files)
    logger.info(f"Total input files for Foldseek: {total_input_files}")

    logger.info(f"Searching for identical structures in {input_dir} with Foldseek")
    logger.info(
        f"tmscore_threshold = {tmscore_threshold}, fident_threshold = {fident_threshold}, n_jobs = {n_jobs}"
    )

    identical_groups = run_foldseek(
        cfg, input_dir, tmscore_threshold, fident_threshold, n_jobs
    )

    identical_groups_json_compatible = [list(group) for group in identical_groups]
    with open(cfg.foldseek.identical_groups, "w") as f:
        json.dump(identical_groups_json_compatible, f)

    logger.info(
        f"Search completed. Found {len(identical_groups)} groups of identical structures."
    )
    return identical_groups
