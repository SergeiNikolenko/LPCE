

import shutil
import tempfile
from pathlib import Path
import sys
from hydra import initialize, compose
from loguru import logger

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from cleanup.remove_dna_rna import remove_dna_rna_from_directory
from cleanup.remove_water import remove_water_from_directory
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from extraction.convert_pdb_to_smiles_sdf import convert_pdb_to_smiles_sdf
from extraction.parse_dict import extract_and_save_complexes_with_ligands
from cleanup.filter_ligands import filter_ligands
# from extraction.separate_complexes import separate_complexes

from utils.send_email import send_email_notification


def test_run_pipeline():
    with initialize(config_path="../config", version_base=None):
        cfg = compose(config_name="config")

    logs_dir = Path("../logs")
    if logs_dir.exists():
        shutil.rmtree(logs_dir)

    logger.remove()
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(cfg.logging.log_file, format="{time} | {level} | {message}", level="INFO")

    test_data_dir = Path("lpce/tests/test_data")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        processed_dir = temp_path / "processed"
        processed_dir.mkdir(parents=True, exist_ok=True)

        ligands_dir = temp_path / "ligands"
        ligands_dir.mkdir(parents=True, exist_ok=True)

        test_cfg = cfg.copy()
        test_cfg.paths.processed_dir = str(processed_dir)
        test_cfg.output_files.complexes_json = str(temp_path / "complexes.json")
        test_cfg.output_files.grouped_complexes_json = str(temp_path / "grouped_complexes.json")
        test_cfg.output_files.site_info_json = str(temp_path / "site_info.json")
        test_cfg.paths.ligands_dir = str(temp_path / "ligands")
        test_cfg.output_files.removed_ligands_summary_json = str(temp_path / "removed_ligands_summary.json")
        test_cfg.output_files.filtered_ligands_json = str(temp_path / "filtered_ligands.json")
        shutil.copytree(test_data_dir, processed_dir, dirs_exist_ok=True)

        logger.info(f"Running tests in temporary directory: {processed_dir}")

        # Проверка наличия PDB файлов в директории
        pdb_files = list(processed_dir.glob("*.pdb"))
        if not pdb_files:
            logger.warning("No PDB files found in the processed directory for testing.")

        remove_dna_rna_from_directory(input_dir=processed_dir, log_file=cfg.logging.dna_rna_removal_log_file)
        remove_water_from_directory(input_dir=processed_dir, log_file=cfg.logging.water_removal_log_file)
        remove_junk_ligands_from_directory(test_cfg)
        convert_pdb_to_smiles_sdf(input_dir=processed_dir, output_dir=ligands_dir)
        extract_and_save_complexes_with_ligands(test_cfg)
        filter_ligands(test_cfg)

        send_email_notification(new_structures='test', email_user=cfg.email.user, email_password=cfg.email.password, receiver_email=cfg.email.recipient, log_file=cfg.logging.email_log_file)
        logger.info(f"Test pipeline completed successfully in {temp_path}")

if __name__ == "__main__":
    test_run_pipeline()
