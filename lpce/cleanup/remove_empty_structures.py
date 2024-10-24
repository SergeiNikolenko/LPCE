import json
from pathlib import Path
from loguru import logger


def remove_unused_pdb_files(cfg) -> dict:
    """
    Removes unused PDB files based on the filtered complexes and keeps only relevant PDB files.

    Args:
        cfg (DictConfig): Configuration object containing paths and log file locations.

    Returns:
        dict: A summary with lists of kept and removed PDB files.
    """
    input_filtered_complexes = Path(cfg.output_files.filtered_ligands_json)
    processed_dir = Path(cfg.paths.processed_dir)

    with open(input_filtered_complexes) as f:
        filtered_complexes = json.load(f)

    proteins_to_keep = set(filtered_complexes.keys())

    pdb_files = list(processed_dir.glob("*.pdb"))
    kept_files = []
    removed_files = []

    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem[3:].lower()

        if pdb_id not in proteins_to_keep:
            pdb_file.unlink()
            removed_files.append(pdb_file.name)
        else:
            kept_files.append(pdb_file.name)

    # Логируем как и раньше
    logger.info("========== Removing Unused PDB Files ==========")
    logger.info(f"Total PDB files in directory: {len(pdb_files):,}")
    logger.info(f"Filtered PDB files to keep: {len(kept_files):,}")
    logger.info(f"PDB files removed: {len(removed_files):,}")

    # Возвращаем только списки оставшихся и удалённых файлов
    return {
        'kept_files': kept_files,
        'removed_files': removed_files
    }
