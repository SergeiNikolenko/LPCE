import json
import os
from collections import defaultdict

from loguru import logger


def extract_het_and_chain_identifiers(filename):
    try:
        start_het = filename.index("bioml_1_") + len("bioml_1_")
        end_het = filename.index("_chains_")
        het_identifier = filename[start_het:end_het]

        start_chain = end_het + len("_chains_")
        end_chain = filename.index("_processed", start_chain)
        chain_identifier = filename[start_chain:end_chain]

        return het_identifier, chain_identifier
    except ValueError:
        return None, None


def get_resolution_from_pdb(pdb_file_path):
    try:
        with open(pdb_file_path) as f:
            for line in f:
                if line.startswith("REMARK   2 RESOLUTION."):
                    parts = line.strip().split()
                    for i, part in enumerate(parts):
                        if "ANGSTROM" in part.upper():
                            try:
                                resolution = float(parts[i - 1])
                                return resolution
                            except ValueError:
                                # logger.warning(f"Unable to convert resolution to float for {pdb_file_path}")
                                return None
                    resolution_str = line.strip().split("RESOLUTION.")[-1].split()[0]
                    try:
                        resolution = float(resolution_str)
                        return resolution
                    except ValueError:
                        # logger.warning(f"Unable to convert resolution to float for {pdb_file_path}")
                        return None
        return None
    except Exception as e:
        logger.error(f"Error reading resolution from {pdb_file_path}: {e}")
        return None


def remove_duplicate_groups(groups):
    unique_groups = []
    groups_sorted = sorted(
        [sorted(group) for group in groups], key=lambda x: (-len(x), x)
    )

    for group in groups_sorted:
        group_set = set(group)
        is_subset = False
        for unique_group in unique_groups:
            if group_set <= unique_group:
                is_subset = True
                break
        if not is_subset:
            unique_groups.append(group_set)
    return unique_groups


def process_groups_with_resolution(identical_groups, pdb_directory):
    processed_groups = []
    for group in identical_groups:
        processed_group = {
            file_name.split("_processed")[0] + "_processed" for file_name in group
        }
        processed_groups.append(processed_group)

    unique_groups = remove_duplicate_groups(processed_groups)
    final_groups = []

    for group in unique_groups:
        het_chain_groups = defaultdict(list)
        for file_name in group:
            het_identifier, chain_identifier = extract_het_and_chain_identifiers(
                file_name
            )
            key = (het_identifier, chain_identifier)
            het_chain_groups[key].append(file_name)

        selected_structures = []
        for key, structures in het_chain_groups.items():
            best_structure = None
            best_resolution = None
            for structure in structures:
                pdb_file_path = os.path.join(pdb_directory, f"{structure}.pdb")
                if not os.path.exists(pdb_file_path):
                    logger.warning(f"File {pdb_file_path} not found.")
                    continue
                resolution = get_resolution_from_pdb(pdb_file_path)
                if resolution is None:
                    logger.warning(f"Resolution not found for {structure}")
                    continue
                if best_resolution is None or resolution < best_resolution:
                    best_resolution = resolution
                    best_structure = structure
            if best_structure:
                selected_structures.append(best_structure)
            else:
                selected_structures.append(structures[0])
        final_groups.append(set(selected_structures))
    return final_groups


def remove_similar_structures(cfg):
    def load_identical_groups(cfg):
        with open(cfg.foldseek.identical_groups) as f:
            identical_groups_json_compatible = json.load(f)
        identical_groups = [set(group) for group in identical_groups_json_compatible]

        return identical_groups

    identical_groups = load_identical_groups(cfg)
    initial_count = sum(len(group) for group in identical_groups)
    final_groups = process_groups_with_resolution(
        identical_groups, cfg.paths.separated_dir
    )

    total_files = set()
    for group in final_groups:
        # Add names with .pdb extension for correct file deletion
        total_files.update({file + ".pdb" for file in group})

    final_count = len(total_files)
    deleted_files = 0

    for filename in os.listdir(cfg.paths.separated_dir):
        file_path = os.path.join(cfg.paths.separated_dir, filename)
        if filename.endswith(".pdb") and filename not in total_files:
            os.remove(file_path)
            deleted_files += 1

    deleted_percentage = (deleted_files / (deleted_files + final_count)) * 100
    logger.info("========== Remove similar structures completed ==========")
    logger.info(f"Initially there were {initial_count} structures.")
    logger.info(f"{final_count} unique structures remain.")
    logger.info(
        f"Deleted {deleted_files} files, which is {deleted_percentage:.2f}% of the total."
    )