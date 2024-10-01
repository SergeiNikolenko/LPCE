import json
from pathlib import Path
from loguru import logger


def remove_unused_pdb_files(cfg):
    input_filtered_complexes = Path(cfg.output_files.filtered_ligands_json)
    processed_dir = Path(cfg.paths.processed_dir)
    log_file = cfg.logging.pdb_cleanup_log_file

    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )

    with open(input_filtered_complexes) as f:
        filtered_complexes = json.load(f)

    proteins_to_keep = set(filtered_complexes.keys())

    pdb_files = list(processed_dir.glob("*.pdb"))
    total_files = len(pdb_files)
    total_kept = 0
    total_removed = 0

    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem[3:].lower()

        if pdb_id not in proteins_to_keep:
            pdb_file.unlink()
            total_removed += 1
        else:
            total_kept += 1

    percent_removed = (total_removed / total_files) * 100 if total_files else 0

    logger.info("========== Removing Unused PDB Files ==========")
    logger.info(f"Total PDB files in directory: {total_files:,}")
    logger.info(f"Filtered PDB files to keep: {total_kept:,}")
    logger.info(f"PDB files removed: {total_removed:,}")
    logger.info(f"Percentage of PDB files removed: {percent_removed:.1f}%")
