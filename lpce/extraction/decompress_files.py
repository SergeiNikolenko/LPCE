import os
from pathlib import Path
import gzip
import shutil
from tqdm import tqdm
from joblib import Parallel, delayed
from config.settings import RAW_DIR, PROCESSED_DIR

def decompress_file(input_file_path: Path, output_file_path: Path) -> str:
    """
    Decompresses a gzip-compressed PDB file and writes it to the specified output path.

    Args:
        input_file_path (Path): The path to the gzip-compressed PDB file.
        output_file_path (Path): The path where the decompressed PDB file will be saved.

    Returns:
        str: A message indicating the outcome of the decompression ("File already exists", 
        "True" for success, or an error message).
    """
    if os.path.exists(output_file_path):
        return "File already exists"

    try:
        with gzip.open(input_file_path, 'rb') as f_in:
            with open(output_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True
    except EOFError as e:
        return f"Error decompressing {input_file_path}: {e}"
    except Exception as e:
        return f"Unexpected error decompressing {input_file_path}: {e}"

def get_file_size_in_gb(file_path: Path) -> float:
    """
    Calculates the size of a file in gigabytes.

    Args:
        file_path (Path): The path to the file.

    Returns:
        float: The size of the file in gigabytes.
    """
    return os.path.getsize(file_path) / (1024 ** 3)

def decompress_pdb_files() -> None:
    """
    Decompresses all gzip-compressed PDB files in the RAW_DIR directory and saves them 
    in the PROCESSED_DIR directory. 

    This function also prints out statistics about the number of files processed, the 
    number of successful decompressions, and the total size of compressed and decompressed files.

    Returns:
        None
    """
    os.makedirs(PROCESSED_DIR, exist_ok=True)

    input_output_paths = []
    for root, _, files in os.walk(RAW_DIR):
        for file in files:
            if file.endswith(".ent.gz"):
                input_file_path = os.path.join(root, file)
                output_file_path = os.path.join(PROCESSED_DIR, file.replace(".ent.gz", ".pdb"))
                input_output_paths.append((input_file_path, output_file_path))

    total_files = len(input_output_paths)

    results = Parallel(n_jobs=-1)(
        delayed(decompress_file)(input_path, output_path) 
        for input_path, output_path in tqdm(input_output_paths, desc="Decompressing files", unit="file", total=total_files)
    )

    successful_files = sum(1 for result in results if result is True)
    skipped_files = sum(1 for result in results if result == "File already exists")
    failed_files = sum(1 for result in results if result not in [True, "File already exists"])

    compressed_size = sum(get_file_size_in_gb(input_path) for input_path, _ in input_output_paths)
    decompressed_size = sum(get_file_size_in_gb(output_path) for _, output_path in input_output_paths if os.path.exists(output_path))

    print(f"Total structures processed: {total_files}")
    print(f"Successfully decompressed: {successful_files}")
    print(f"Skipped (already exists): {skipped_files}")
    print(f"Failed to decompress: {failed_files}")
    print(f"Total size of compressed files: {compressed_size:.2f} GB")
    print(f"Total size of decompressed files: {decompressed_size:.2f} GB")
