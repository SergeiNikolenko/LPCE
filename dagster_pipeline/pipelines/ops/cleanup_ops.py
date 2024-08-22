from dagster import op, get_dagster_logger, In
from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
import json

@op(ins={"input_dir": In(dagster_type=str, default_value="/mnt/ligandpro/db/PDB/pdb2/processed")})
def run_remove_junk_ligands(context, input_dir: str = None):
    logger = get_dagster_logger()

    input_directory = Path(input_dir)

    with open('./data/trash_ligands.json', 'r') as file:
        junk_ligands = json.load(file)
    
    def remove_junk_ligands_from_file(input_file_path):
        try:
            # Чтение и фильтрация строк в памяти
            filtered_lines = []
            with open(input_file_path, 'r') as f_in:
                for line in f_in:
                    if not (line.startswith("HETATM") and any(ligand in line for ligand in junk_ligands)):
                        filtered_lines.append(line)
            
            with open(input_file_path, 'w') as f_out:
                f_out.writelines(filtered_lines)
            
            return True
        except Exception as e:
            logger.error(f"Failed to process {input_file_path}: {e}")
            return f"Error processing {input_file_path}: {e}"

    def remove_junk_ligands_from_directory(input_directory, n_jobs=-1):
        pdb_files = list(input_directory.rglob("*.pdb"))

        total_files = len(pdb_files)
        logger.info(f"Found {total_files} PDB files in {input_directory}")

        results = Parallel(n_jobs=n_jobs)(
            delayed(remove_junk_ligands_from_file)(pdb_file) for pdb_file in tqdm(pdb_files, desc="Removing junk ligands", unit="file", total=total_files)
        )

        successful_files = sum(1 for result in results if result is True)
        failed_files = sum(1 for result in results if result is not True)

        logger.info(f"Total structures processed: {total_files}")
        logger.info(f"Successfully processed: {successful_files}")
        logger.info(f"Failed to process: {failed_files}")

    remove_junk_ligands_from_directory(input_directory)
