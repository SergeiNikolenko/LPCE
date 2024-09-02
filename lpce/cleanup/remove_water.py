from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
import subprocess
from config.settings import PROCESSED_DIR

def remove_water_from_directory() -> None:
    """
    Processes all PDB files in the specified directory, removing water molecules from each file.

    This function utilizes a compiled C program (`remove_water`) to perform the water removal.
    The results are printed to the console, summarizing the number of files processed.

    Returns:
        None
    """
    input_directory = Path(PROCESSED_DIR)
    pdb_files = list(input_directory.rglob("*.pdb"))

    total_files = len(pdb_files)
    print(f"Found {total_files} PDB files in {input_directory}")

    executable_path = "cleanup/remove_water"

    for pdb_file in tqdm(pdb_files, desc="Removing water", unit="file", total=total_files):
        subprocess.run([executable_path, str(pdb_file)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print(f"Total structures processed: {total_files}")
