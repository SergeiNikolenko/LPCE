import gzip
import shutil
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def decompress_file(input_file_path: Path, output_file_path: Path) -> str:
    if output_file_path.exists():
        return "File already exists"
    try:
        with input_file_path.open("rb") as f_in:
            with output_file_path.open("wb") as f_out:
                shutil.copyfileobj(gzip.GzipFile(fileobj=f_in), f_out)
        return True
    except (EOFError, OSError) as e:
        logger.error(f"Error decompressing {input_file_path}: {e}")
        return f"Error decompressing {input_file_path}: {e}"


def get_file_size_in_gb(file_path: Path) -> float:
    return file_path.stat().st_size / (1024**3)


def decompress_pdb_files(cfg: object) -> None:
    logger.info("\n========== Decompressing Files ==========")
    output_dir = Path(cfg.paths.processed_dir)
    n_jobs = cfg.n_jobs
    output_dir.mkdir(parents=True, exist_ok=True)

    input_output_paths = [
        (input_file, (output_dir / input_file.stem).with_suffix(".pdb"))
        for input_file in Path(cfg.paths.raw_dir).rglob("*.ent.gz")
    ]
    total_files = len(input_output_paths)

    logger.info(f"Total .ent.gz files to decompress: {total_files}")

    results = Parallel(n_jobs=n_jobs, backend="threading")(
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
    failed_files = total_files - successful_files - skipped_files

    compressed_size = sum(
        get_file_size_in_gb(input_path) for input_path, _ in input_output_paths
    )
    decompressed_size = sum(
        get_file_size_in_gb(output_path)
        for _, output_path in input_output_paths
        if output_path.exists()
    )

    logger.info(f"Total files processed: {total_files}")
    logger.info(f"Successfully decompressed: {successful_files}")
    logger.info(f"Skipped (already exists): {skipped_files}")
    logger.info(f"Failed to decompress: {failed_files}")
    logger.info(f"Total size of compressed files: {compressed_size:.2f} GB")
    logger.info(f"Total size of decompressed files: {decompressed_size:.2f} GB")
