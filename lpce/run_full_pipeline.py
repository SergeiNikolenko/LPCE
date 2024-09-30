from hydra import initialize, compose
from loguru import logger
from pathlib import Path
import sys
import shutil

from cleanup.filter_ligands import filter_ligands
from cleanup.remove_dna_rna import remove_dna_rna_from_directory
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from cleanup.remove_water import remove_water_from_directory
from extraction.convert_pdb_to_smiles_sdf import convert_pdb_to_smiles_sdf
from extraction.decompress_files import decompress_pdb_files
from extraction.extract_complexes import extract_complexes
from extraction.parse_dict import extract_and_save_complexes_with_ligands
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
        cfg.logging.log_file, format="{time} | {level} | {message}", level="INFO"
    )

    raw_dir = Path(cfg.paths.raw_dir)
    processed_dir = Path(cfg.paths.processed_dir)

    new_structures = extract_complexes(
        raw_dir=raw_dir, rsync_port=cfg.rsync.port, rsync_host=cfg.rsync.host
    )

    decompress_pdb_files(input_dir=raw_dir, output_dir=processed_dir)
    remove_dna_rna_from_directory(
        input_dir=processed_dir, log_file=cfg.logging.dna_rna_removal_log_file
    )
    remove_water_from_directory(
        input_dir=processed_dir, log_file=cfg.logging.water_removal_log_file
    )
    remove_junk_ligands_from_directory(cfg)
    convert_pdb_to_smiles_sdf(input_dir=processed_dir, output_dir=cfg.paths.ligands_dir)
    extract_and_save_complexes_with_ligands(cfg)
    filter_ligands(cfg)

    send_email_notification(
        new_structures=new_structures,
        email_user=cfg.email.user,
        email_password=cfg.email.password,
        receiver_email=cfg.email.recipient,
        log_file=cfg.logging.email_log_file,
    )


if __name__ == "__main__":
    main()
