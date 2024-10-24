import shutil
import sys
from pathlib import Path

from cleanup.filter_ligands import filter_ligands
from cleanup.remove_dna_rna import remove_dna_rna_from_directory
from cleanup.remove_empty_structures import remove_unused_pdb_files
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from cleanup.remove_multiple_models import remove_multiple_models_from_directory
from cleanup.remove_water import remove_water_from_directory
from extraction.convert_pdb_to_smiles_sdf import convert_pdb_to_smiles_sdf
from extraction.decompress_files import decompress_pdb_files
from extraction.extract_complexes import extract_complexes
from extraction.parse_dict import extract_and_save_complexes_with_ligands
from hydra import compose, initialize
from loguru import logger
from pdb_manipulations.protein_ligand_separator import protein_ligand_separator
from pdb_manipulations.split_bioml import bioml_split
from utils.send_email import send_email_notification


def main():
    with initialize(config_path="config", version_base=None):
        cfg = compose(config_name="config")

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
    (decompress_pdb_files(cfg),)
    remove_dna_rna_from_directory(cfg)
    remove_multiple_models_from_directory(cfg)
    remove_water_from_directory(cfg)
    remove_junk_ligands_from_directory(cfg)
    convert_pdb_to_smiles_sdf(cfg)
    extract_and_save_complexes_with_ligands(cfg)
    filter_ligands(cfg)
    remove_unused_pdb_files(cfg)
    bioml_split(cfg)
    protein_ligand_separator(cfg)
    send_email_notification(cfg)


if __name__ == "__main__":
    main()
