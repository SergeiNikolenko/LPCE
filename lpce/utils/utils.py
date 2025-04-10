import json

from loguru import logger


def save_removed_files_to_json(dna_rna, models, unused, not_buried, output_path):
    """
    Combines and saves removed files information from three functions into a single JSON file.

    Args:
        dna_rna (dict): Result from remove_dna_rna_from_directory function.
        models (dict): Result from remove_multiple_models_from_directory function.
        unused (dict): Result from remove_unused_pdb_files function.
        output_path (Path): Path to save the combined JSON file.
    """
    combined_data = {
        "dna_rna_removed_files": dna_rna["removed_files"],
        "multiple_models_removed_files": models["modified_files"],
        "unused_pdb_removed_files": unused["removed_files"],
        "not_buried_removed_files": not_buried["removed_files"],
    }
    with open(output_path, "w") as f:
        json.dump(combined_data, f, indent=4)
    logger.info(f"Removed files saved to {output_path}")
