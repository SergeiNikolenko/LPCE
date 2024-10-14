import gzip
import shutil
import sys
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def decompress_file(input_file_path: Path, output_file_path: Path) -> str:
    """
    Decompress a gzip-compressed PDB file and writes it to the specified output path.

    Args:
        input_file_path (Path): The path to the gzip-compressed PDB file.
        output_file_path (Path): The path where the decompressed PDB file will be saved.

    Returns:
        str: "File already exists" if the output file already exists, "True" for success, or an error message if an error occurs.
    """
    if output_file_path.exists():
        return "File already exists"

    try:
        with gzip.open(input_file_path, "rb") as f_in:
            with open(output_file_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True
    except EOFError as e:
        logger.error(f"Error decompressing {input_file_path}: {e}")
        return f"Error decompressing {input_file_path}: {e}"
    except Exception as e:
        logger.error(f"Unexpected error decompressing {input_file_path}: {e}")
        return f"Unexpected error decompressing {input_file_path}: {e}"


def get_file_size_in_gb(file_path: Path) -> float:
    """
    Get the size of the specified file in gigabytes.

    Args:
        file_path (Path): The path to the file.

    Returns:
        float: The size of the file in gigabytes.
    """
    return file_path.stat().st_size / (1024**3)


def decompress_pdb_files(input_dir: Path, output_dir: Path, log_file: str) -> None:
    """
    Decompress all gzip-compressed PDB files in the input directory and save them
    in the specified output directory. Logs the progress and statistics to the specified log file.

    Args:
        input_dir (Path): Directory containing the gzip-compressed files.
        output_dir (Path): Directory to store the decompressed files.
        log_file (str): Path to the log file for logging the decompression process.

    Returns:
        None
    """
    # Add a separate log file for the decompression process
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )
    logger.info("========== Decompressing Files ==========")

    output_dir.mkdir(parents=True, exist_ok=True)

    input_output_paths = [
        (input_file, output_dir / input_file.name.replace(".ent.gz", ".pdb"))
        for input_file in input_dir.rglob("*.ent.gz")
    ]
    total_files = len(input_output_paths)

    logger.info(f"Total .ent.gz files to decompress: {total_files}")

    results = Parallel(n_jobs=-1)(
        delayed(decompress_file)(input_path, output_path)
        for input_path, output_path in tqdm(
            input_output_paths,
            desc="Decompressing files",
            unit="file",
            total=total_files,
        )
    )

    successful_files = sum(1 for result in results if result is True)
    skipped_files = sum(1 for result in results if result == "File already exists")
    failed_files = sum(
        1 for result in results if result not in [True, "File already exists"]
    )

    compressed_size = sum(
        get_file_size_in_gb(input_path) for input_path, _ in input_output_paths
    )
    decompressed_size = sum(
        get_file_size_in_gb(output_path)
        for _, output_path in input_output_paths
        if output_path.exists()
    )

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully decompressed: {successful_files}")
    logger.info(f"Skipped (already exists): {skipped_files}")
    logger.info(f"Failed to decompress: {failed_files}")
    logger.info(f"Total size of compressed files: {compressed_size:.2f} GB")
    logger.info(f"Total size of decompressed files: {decompressed_size:.2f} GB")
