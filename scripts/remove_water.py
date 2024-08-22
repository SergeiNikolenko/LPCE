import os
from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
import logging

from dagster import repository


# Настройка логирования
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def remove_water_from_file(input_file_path, output_directory):
    output_file_path = Path(output_directory) / (Path(input_file_path).stem + "_nowater.pdb")
    
    if output_file_path.exists():
        logger.info(f"Output file {output_file_path} already exists. Skipping.")
        return "File already exists"

    try:
        with open(input_file_path, 'r') as f_in, open(output_file_path, 'w') as f_out:
            for line in f_in:
                # Удаляем строки, относящиеся к воде
                if not line.startswith("HETATM") or "HOH" not in line:
                    f_out.write(line)
        logger.info(f"Processed {input_file_path}, saved to {output_file_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to process {input_file_path}: {e}")
        return f"Error processing {input_file_path}: {e}"

def remove_water_from_directory(input_directory, output_directory, n_jobs=-1):
    input_directory = Path(input_directory)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    
    pdb_files = list(input_directory.rglob("*.pdb"))

    total_files = len(pdb_files)
    logger.info(f"Found {total_files} PDB files in {input_directory}")

    results = Parallel(n_jobs=n_jobs)(
        delayed(remove_water_from_file)(pdb_file, output_directory) for pdb_file in tqdm(pdb_files, desc="Removing water", unit="file", total=total_files)
    )

    successful_files = sum(1 for result in results if result is True)
    skipped_files = sum(1 for result in results if result == "File already exists")
    failed_files = sum(1 for result in results if result not in [True, "File already exists"])

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully processed: {successful_files}")
    logger.info(f"Skipped (already exists): {skipped_files}")
    logger.info(f"Failed to process: {failed_files}")

if __name__ == "__main__":
    input_dir = "/mnt/ligandpro/db/PDB/pdb2/processed/decompressed"  # Укажите путь к директории с PDB файлами
    output_dir = "/mnt/ligandpro/db/PDB/pdb2/processed/no_water"  # Директория для сохранения файлов без воды
    remove_water_from_directory(input_dir, output_dir)



import os
from dagster import job, op, get_dagster_logger
import subprocess
from pathlib import Path

@op
def run_remove_water():
    logger = get_dagster_logger()
    script_path = Path(__file__).resolve().parent.parent / "scripts" / "remove_water.py"
    
    input_dir = "/mnt/ligandpro/db/PDB/pdb2/processed/decompressed"
    output_dir = "/mnt/ligandpro/db/PDB/pdb2/processed/no_water" 

    process = subprocess.Popen(
        ["python", str(script_path), input_dir, output_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    for line in process.stdout:
        logger.info(line.strip())
    
    process.wait()
    
    if process.returncode != 0:
        for line in process.stderr:
            logger.error(line.strip())
        raise Exception(f"Script exited with return code {process.returncode}")
    
    logger.info("Water removal completed successfully")

@job
def water_removal_job():
    run_remove_water()


@repository
def my_pipeline_repository():
    return [pdb_extraction_job, water_removal_job]
