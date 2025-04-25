import argparse
import shutil
import sys
from pathlib import Path

from cleanup import (
    remove_dna_rna_from_directory,
    remove_junk_ligands_from_directory,
    remove_multiple_models_from_directory,
    remove_water_from_directory,
)
from extraction import decompress_pdb_files, extract_complexes
from hydra import compose, initialize
from loguru import logger
from pdb_manipulations import (
    add_h_to_ligands,
    bioml_split,
    create_final_files,
    protein_ligand_separator,
    remove_not_buried_ligands,
    split_overlapping_ligands,
)
from utils import clean_multiple_paths, save_removed_files_to_json


def main(config_name):
    config_path = "config"
    with initialize(config_path=config_path, version_base=None):
        cfg = compose(config_name=config_name)
    
    logs_dir = Path("logs")
    if logs_dir.exists():
        shutil.rmtree(logs_dir)

    logger.remove()
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(
        cfg.logging.log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )

    logger.debug(f"Config: {cfg}")

    extract_complexes(cfg)
    decompress_pdb_files(cfg)

    dna_rna = remove_dna_rna_from_directory(cfg)
    models = remove_multiple_models_from_directory(cfg)

    remove_water_from_directory(cfg)

    remove_junk_ligands_from_directory(cfg)
    
    # convert_pdb_to_smiles_sdf(cfg)
    # extract_and_save_complexes_with_ligands(cfg)
    # filter_ligands(cfg)
    # unused = remove_unused_pdb_files(cfg)
    unused = {"removed_files": []}
    bioml_split(cfg)
    protein_ligand_separator(cfg)

    clean_multiple_paths(cfg)
    find_duplicates_foldseek(cfg)
    # remove_similar_structures(cfg)
    not_buried = remove_not_buried_ligands(cfg)
    split_overlapping_ligands(cfg)
    add_h_to_ligands(cfg)
    create_final_files(cfg)

    json_output_path = Path(cfg.output_files.removed_files_json)
    save_removed_files_to_json(dna_rna, models, unused, not_buried, json_output_path)

    # send_email_notification(cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the full pipeline with a specified config name."
    )
    parser.add_argument(
        "config_name",
        type=str,
        help="Name of the configuration file (without .yaml extension)",
    )
    args = parser.parse_args()
    main(args.config_name)
