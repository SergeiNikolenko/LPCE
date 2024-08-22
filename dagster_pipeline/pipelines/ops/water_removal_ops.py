from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def remove_water_from_file(input_file_path):
    try:
        with open(input_file_path, 'r') as f_in:
            filtered_lines = [line for line in f_in if not (line.startswith("HETATM") and "HOH" in line)]
        
        with open(input_file_path, 'w') as f_out:
            f_out.writelines(filtered_lines)

        logger.info(f"Processed and overwrote {input_file_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to process {input_file_path}: {e}")
        return f"Error processing {input_file_path}: {e}"

def remove_water_from_directory(input_directory, n_jobs=-1):
    input_directory = Path(input_directory)
    pdb_files = list(input_directory.rglob("*.pdb"))

    total_files = len(pdb_files)
    logger.info(f"Found {total_files} PDB files in {input_directory}")

    results = Parallel(n_jobs=n_jobs)(
        delayed(remove_water_from_file)(pdb_file) for pdb_file in tqdm(pdb_files, desc="Removing water", unit="file", total=total_files)
    )

    successful_files = sum(1 for result in results if result is True)
    failed_files = sum(1 for result in results if result is not True)

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully processed: {successful_files}")
    logger.info(f"Failed to process: {failed_files}")


remove_water_from_directory("/mnt/ligandpro/db/PDB/pdb2/processed")
