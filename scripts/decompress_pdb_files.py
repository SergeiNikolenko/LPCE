import gzip
import shutil
import os
from tqdm import tqdm
from joblib import Parallel, delayed

def decompress_file(input_file_path, output_file_path):
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

def get_file_size_in_gb(file_path):
    return os.path.getsize(file_path) / (1024 ** 3)

def decompress_pdb_files(input_directory, output_directory, n_jobs=-1):
    os.makedirs(output_directory, exist_ok=True)

    input_output_paths = []
    for root, _, files in os.walk(input_directory):
        for file in files:
            if file.endswith(".ent.gz"):
                input_file_path = os.path.join(root, file)
                output_file_path = os.path.join(output_directory, file.replace(".ent.gz", ".pdb"))
                input_output_paths.append((input_file_path, output_file_path))

    total_files = len(input_output_paths)

    results = Parallel(n_jobs=n_jobs)(
        delayed(decompress_file)(input_path, output_path) for input_path, output_path in tqdm(input_output_paths, desc="Decompressing files", unit="file", total=total_files)
    )

    successful_files = sum(1 for result in results if result is True)
    skipped_files = sum(1 for result in results if result == "File already exists")
    failed_files = sum(1 for result in results if result not in [True, "File already exists"])

    for result in results:
        if result not in [True, "File already exists"]:
            print(result)

    compressed_size = sum(get_file_size_in_gb(input_path) for input_path, _ in input_output_paths)
    decompressed_size = sum(get_file_size_in_gb(output_path) for _, output_path in input_output_paths if os.path.exists(output_path))

    print(f"Total structures processed: {total_files}")
    print(f"Successfully decompressed: {successful_files}")
    print(f"Skipped (already exists): {skipped_files}")
    print(f"Failed to decompress: {failed_files}")
    print(f"Total size of compressed files: {compressed_size:.2f} GB")
    print(f"Total size of decompressed files: {decompressed_size:.2f} GB")

if __name__ == "__main__":
    input_dir = "/mnt/ligandpro/db/PDB/pdb2/raw/pdb"
    output_dir = "/mnt/ligandpro/db/PDB/pdb2/processed/decompressed"
    decompress_pdb_files(input_dir, output_dir)
