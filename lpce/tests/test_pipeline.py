import shutil
import tempfile
import warnings
from pathlib import Path
from typing import Any

from Bio import BiopythonDeprecationWarning
from hydra import compose, initialize
from loguru import logger

from lpce.cleanup import (
    filter_ligands,
    remove_dna_rna_from_directory,
    remove_junk_ligands_from_directory,
    remove_multiple_models_from_directory,
    remove_unused_pdb_files,
    remove_water_from_directory,
)
from lpce.extraction import (
    convert_pdb_to_smiles_sdf,
    extract_and_save_complexes_with_ligands,
)
from lpce.pdb_manipulations import (
    add_h_to_ligands,
    bioml_split,
    create_final_files,
    find_duplicates_foldseek,
    protein_ligand_separator,
    remove_not_buried_ligands,
    remove_similar_structures,
    split_overlapping_ligands,
)
from lpce.utils import (
    clean_multiple_paths,
    cleanup_directories,
    copy_results_to_final,
    save_removed_files_to_json,
    setup_logging,
    setup_test_directories,
    update_test_config,
)

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


def run_pipeline_steps(test_cfg: dict[str, Any]) -> list[str]:
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
    split_overlapping_ligands(test_cfg)
    add_h_to_ligands(test_cfg)
    create_final_files(test_cfg)

    return [dna_rna, models, unused, not_buried]


def test_run_pipeline(config_name: str) -> None:
    with initialize(config_path="../config", version_base=None):
        cfg = compose(config_name=config_name)

    setup_logging(cfg)

    test_data_dir = Path("lpce/tests/test_data")
    final_dirs = {
        "processed": Path("./lpce/tests/processed"),
        "bioml": Path("./lpce/tests/bioml"),
        "separated": Path("./lpce/tests/separated"),
        "final": Path("./lpce/tests/final"),
    }

    cleanup_directories(final_dirs.values())

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        directories = setup_test_directories(temp_path)
        test_cfg = update_test_config(cfg, directories)

        shutil.copytree(test_data_dir, directories["processed"], dirs_exist_ok=True)
        logger.info(f"Running tests in temporary directory: {directories['processed']}")

        pdb_files = list(directories["processed"].glob("*.pdb"))
        logger.info(f"PDB files found: {len(pdb_files)}")
        if not pdb_files:
            logger.warning("No PDB files found in the processed directory for testing.")

        removed_files = run_pipeline_steps(test_cfg)

        json_output_path = Path("data/removed_files_tests.json")
        save_removed_files_to_json(
            removed_files[0],
            removed_files[1],
            removed_files[2],
            removed_files[3],
            json_output_path,
        )

        copy_results_to_final(directories["processed"], test_cfg, final_dirs)
        logger.info("DONE! Test pipeline completed successfully")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the test pipeline with a specified config name."
    )
    parser.add_argument(
        "config_name",
        type=str,
        help="Name of the configuration file (without .yaml extension)",
    )
    args = parser.parse_args()
    test_run_pipeline(args.config_name)
