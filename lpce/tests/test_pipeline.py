import shutil
import sys
import tempfile
from pathlib import Path
import argparse

from hydra import compose, initialize
from loguru import logger

sys.path.append(str(Path(__file__).resolve().parent.parent))

import warnings
from Bio import BiopythonDeprecationWarning
from cleanup.filter_ligands import filter_ligands
from cleanup.remove_dna_rna import remove_dna_rna_from_directory
from cleanup.remove_empty_structures import remove_unused_pdb_files
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from cleanup.remove_multiple_models import remove_multiple_models_from_directory
from cleanup.remove_water import remove_water_from_directory
from extraction.convert_pdb_to_smiles_sdf import convert_pdb_to_smiles_sdf
from extraction.parse_dict import extract_and_save_complexes_with_ligands
from pdb_manipulations.add_h_to_ligands import add_h_to_ligands
from pdb_manipulations.foldseek import find_duplicates_foldseek
from pdb_manipulations.protein_ligand_separator import protein_ligand_separator
from pdb_manipulations.remove_not_buried_ligands import remove_not_buried_ligands
from pdb_manipulations.remove_similar_structures import remove_similar_structures
from pdb_manipulations.split_bioml import bioml_split
from utils.clean_names import clean_multiple_paths
from utils.utils import save_removed_files_to_json

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


def test_run_pipeline(config_name):
    with initialize(config_path="../config", version_base=None):
        cfg = compose(config_name=config_name)

    log_file_path = Path(cfg.logging.log_file_tests)
    if log_file_path.exists():
        log_file_path.unlink()

    logger.remove()
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(
        cfg.logging.log_file_tests,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )

    test_data_dir = Path("lpce/tests/test_data")

    processed_path = Path("./lpce/tests/processed")
    if processed_path.exists():
        shutil.rmtree(processed_path)

    bioml_path = Path("./lpce/tests/bioml")
    if bioml_path.exists():
        shutil.rmtree(bioml_path)

    separated_path = Path("./lpce/tests/separated")
    if separated_path.exists():
        shutil.rmtree(separated_path)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        processed_dir = temp_path / "processed"
        processed_dir.mkdir(parents=True, exist_ok=True)

        ligands_dir = temp_path / "ligands"
        ligands_dir.mkdir(parents=True, exist_ok=True)

        logs_dir = temp_path / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

        test_cfg = cfg.copy()
        test_cfg.paths.processed_dir = str(processed_dir)
        test_cfg.paths.bioml_dir = str(temp_path / "bioml")
        test_cfg.paths.ligands_dir = str(temp_path / "ligands")
        test_cfg.paths.separated_dir = str(temp_path / "separated")
        test_cfg.paths.separated_fixed = str(temp_path / "separated_fixed")
        test_cfg.logging.log_file = str(logs_dir / "foldseek.log")

        shutil.copytree(test_data_dir, processed_dir, dirs_exist_ok=True)

        logger.info(f"Running tests in temporary directory: {processed_dir}")

        pdb_files = list(processed_dir.glob("*.pdb"))
        logger.info(f"PDB files found: {len(pdb_files)}")
        if not pdb_files:
            logger.warning("No PDB files found in the processed directory for testing.")

        dna_rna = remove_dna_rna_from_directory(test_cfg)
        models = remove_multiple_models_from_directory(test_cfg)
        remove_water_from_directory(test_cfg)
        remove_junk_ligands_from_directory(test_cfg)
        convert_pdb_to_smiles_sdf(test_cfg)
        extract_and_save_complexes_with_ligands(test_cfg)
        filter_ligands(test_cfg)
        unused = remove_unused_pdb_files(test_cfg)
        bioml_split(test_cfg)
        protein_ligand_separator(test_cfg)
        clean_multiple_paths(test_cfg)
        find_duplicates_foldseek(test_cfg)
        remove_similar_structures(test_cfg)
        not_buried = remove_not_buried_ligands(test_cfg)
        add_h_to_ligands(test_cfg)

        json_output_path = Path("data/removed_files_tests.json")
        save_removed_files_to_json(
            dna_rna, models, unused, not_buried, json_output_path
        )

        shutil.copytree(processed_dir, processed_path, dirs_exist_ok=True)
        shutil.copytree(test_cfg.paths.bioml_dir, bioml_path, dirs_exist_ok=True)
        shutil.copytree(
            test_cfg.paths.separated_dir, separated_path, dirs_exist_ok=True
        )

        logger.info("DONE! Test pipeline completed successfully")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the test pipeline with a specified config name.")
    parser.add_argument(
        "config_name",
        type=str,
        help="Name of the configuration file (without .yaml extension)"
    )
    args = parser.parse_args()
    test_run_pipeline(args.config_name)
