import json
from pathlib import Path

from loguru import logger


def remove_unused_pdb_files(cfg) -> dict:
    """
    Removes unused PDB files based on the filtered complexes and keeps only relevant PDB files.

    Args:
        cfg (DictConfig): Configuration object containing paths and log file locations.

    Returns:
        dict: A summary with a list of removed PDB files without prefixes or extensions.
    """
    input_filtered_complexes = Path(cfg.output_files.filtered_ligands_json)
    processed_dir = Path(cfg.paths.processed_dir)

    with open(input_filtered_complexes) as f:
        filtered_complexes = json.load(f)

    if not filtered_complexes:
        raise ValueError("Filtered complexes JSON is empty. Check the input file.")

    proteins_to_keep = set(filtered_complexes.keys())

    pdb_files = list(processed_dir.glob("*.pdb"))
    kept_files = []
    removed_files = []

    for pdb_file in pdb_files:
        file_stem = pdb_file.stem.lower()
        if file_stem.startswith("pdb") and len(file_stem) > 3:
            pdb_id = file_stem[3:]  # Убираем префикс "pdb"
        elif len(file_stem) == 4:
            pdb_id = file_stem  # Это валидный PDB ID без префикса
        else:
            logger.warning(f"Skipping unexpected file format: {pdb_file.name}")
            continue

        if pdb_id not in proteins_to_keep:
            logger.debug(f"Removing PDB file: {pdb_file.name}")
            pdb_file.unlink()
            removed_files.append(pdb_id)
        else:
            kept_files.append(pdb_id)

    logger.info("\n========== Removing Unused PDB Files ==========")
    logger.info(f"Total PDB files in directory: {len(pdb_files):,}")
    logger.info(f"Filtered PDB files to keep: {len(kept_files):,}")
    logger.info(f"PDB files removed: {len(removed_files):,}")

    return {"removed_files": removed_files}
