import json
from pathlib import Path

from loguru import logger


def filter_ligands(cfg):
    """
    Filters ligands based on their presence in site information and saves updated complexes.

    Args:
        cfg (DictConfig): Configuration object containing paths and log file locations.

    Returns:
        dict: A summary of the filtering process, including counts of analyzed proteins,
              ligands, ligands deleted, and remaining ligands.
    """
    input_complexes = Path(cfg.output_files.grouped_complexes_json)
    input_site_info = Path(cfg.output_files.site_info_json)
    output_filtered = Path(cfg.output_files.filtered_ligands_json)

    with open(input_complexes) as f:
        grouped_complexes = json.load(f)

    with open(input_site_info) as f:
        site_info = json.load(f)

    total_proteins_grouped = len(grouped_complexes)
    total_ligands_grouped = sum(
        len(ligands)
        for chains in grouped_complexes.values()
        for ligands in chains.values()
    )

    total_proteins_site_info = len(site_info)
    ligand_intersections = 0
    ligand_deletions = 0
    remaining_ligands = 0

    for protein, chains in grouped_complexes.items():
        if protein not in site_info:
            remaining_ligands += sum(len(ligands) for ligands in chains.values())
            continue  # Skip protein if it's not in site_info.json
        for chain, ligands in chains.items():
            ligands_to_keep = []
            for ligand_info in ligands:
                ligand_name = ligand_info["ligand"]
                if (
                    ligand_name in site_info[protein]
                    and chain in site_info[protein][ligand_name]
                ):
                    ligands_to_keep.append(ligand_info)
                    ligand_intersections += 1
                else:
                    ligand_deletions += 1
            grouped_complexes[protein][chain] = ligands_to_keep
            remaining_ligands += len(ligands_to_keep)

    with open(output_filtered, "w") as f:
        json.dump(grouped_complexes, f, indent=4)

    percent_deleted = (
        (ligand_deletions / total_ligands_grouped) * 100 if total_ligands_grouped else 0
    )

    logger.info("\n========== Filtering Ligands ==========")
    logger.info(f"Total proteins analyzed: {total_proteins_grouped:,}")
    logger.info(f"Total ligands analyzed: {total_ligands_grouped:,}")
    logger.info(f"Proteins with site info available: {total_proteins_site_info:,}")
    logger.info(f"Relevant ligands found in sites: {ligand_intersections:,}")
    logger.info(f"Ligands removed during filtering: {ligand_deletions:,}")
    logger.info(f"Percentage of ligands removed: {percent_deleted:.1f}%")
    logger.info(f"Ligands remaining after filtering: {remaining_ligands:,}")

    return {
        "total_proteins_grouped": total_proteins_grouped,
        "total_ligands_grouped": total_ligands_grouped,
        "total_proteins_site_info": total_proteins_site_info,
        "ligand_intersections": ligand_intersections,
        "ligand_deletions": ligand_deletions,
        "percent_deleted": percent_deleted,
        "remaining_ligands": remaining_ligands,
    }
