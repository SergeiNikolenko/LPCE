from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
from config.settings import PROCESSED_DIR

def remove_water_from_file(input_file_path):
    try:
        with open(input_file_path, 'r') as f_in:
            filtered_lines = [line for line in f_in if not (line.startswith("HETATM") and "HOH" in line)]

        with open(input_file_path, 'w') as f_out:
            f_out.writelines(filtered_lines)

        print(f"Processed and overwrote {input_file_path}")
        return True
    except Exception as e:
        print(f"Failed to process {input_file_path}: {e}")
        return f"Error processing {input_file_path}: {e}"

def remove_water_from_directory():
    input_directory = Path(PROCESSED_DIR)
    pdb_files = list(input_directory.rglob("*.pdb"))

    total_files = len(pdb_files)
    print(f"Found {total_files} PDB files in {input_directory}")

    results = Parallel(n_jobs=-1)(
        delayed(remove_water_from_file)(pdb_file) for pdb_file in tqdm(pdb_files, desc="Removing water", unit="file", total=total_files)
    )

    successful_files = sum(1 for result in results if result is True)
    failed_files = sum(1 for result in results if result is not True)

    print(f"Total structures processed: {total_files}")
    print(f"Successfully processed: {successful_files}")
    print(f"Failed to process: {failed_files}")
