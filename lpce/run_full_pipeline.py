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
from pdb_manipulations.foldseek import find_duplicates_foldseek
from pdb_manipulations.remove_similar_structures import remove_similar_structures

from utils.clean_names import clean_multiple_paths
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

    _ = extract_complexes(cfg)
    _ = decompress_pdb_files(cfg)
    _ = remove_dna_rna_from_directory(cfg)
    _ = remove_multiple_models_from_directory(cfg)
    _ = remove_water_from_directory(cfg)
    _ = remove_junk_ligands_from_directory(cfg)
    _ = convert_pdb_to_smiles_sdf(cfg)
    _ = extract_and_save_complexes_with_ligands(cfg)
    _ = filter_ligands(cfg)
    _ = remove_unused_pdb_files(cfg)
    _ = bioml_split(cfg)
    _ = protein_ligand_separator(cfg)
    _ = send_email_notification(cfg)
    _ = clean_multiple_paths(cfg)
    _ = find_duplicates_foldseek(cfg)
    _ = remove_similar_structures(cfg)

if __name__ == "__main__":
    main()
