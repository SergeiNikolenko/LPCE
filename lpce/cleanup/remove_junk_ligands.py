import json
from collections import defaultdict
from pathlib import Path

from config.settings import PROCESSED_DIR
from joblib import Parallel, delayed
from tqdm import tqdm


def remove_junk_ligands_from_file(input_file_path: Path, junk_ligands: set) -> dict:
    """
    Removes junk ligands from a PDB file and writes the cleaned content back to the same file.

    Args:
        input_file_path (Path): The path to the PDB file to be processed.
        junk_ligands (set): A set of ligand IDs to be removed.

    Returns:
        dict: A dictionary with the count of removed ligands.
    """
    try:
        ligand_counts = defaultdict(int)
        filtered_lines = []

        with open(input_file_path) as f_in:
            for line in f_in:
                if line.startswith("HETATM"):
                    ligand_id = line[17:20].strip()
                    if ligand_id in junk_ligands or len(ligand_id) < 3:
                        ligand_counts[ligand_id] += 1
                        continue
                filtered_lines.append(line)

        with open(input_file_path, "w") as f_out:
            f_out.writelines(filtered_lines)

        return ligand_counts
    except Exception as e:
        print(f"Failed to process {input_file_path}: {e}")
        return {"error": f"Error processing {input_file_path}: {e}"}


def remove_junk_ligands_from_directory() -> None:
    """
    Processes all PDB files in the specified directory, removing junk ligands from each file.

    The results are summarized and saved to a JSON file. The function also prints out statistics
    about the number of files processed and the number of ligands removed.

    Returns:
        None
    """
    input_directory = Path(PROCESSED_DIR)
    with open("data/trash_ligands.json") as file:
        junk_ligands = set(json.load(file))

    pdb_files = list(input_directory.rglob("*.pdb"))

    total_files = len(pdb_files)
    print(f"Found {total_files} PDB files in {input_directory}")

    results = Parallel(n_jobs=-1)(
        delayed(remove_junk_ligands_from_file)(pdb_file, junk_ligands)
        for pdb_file in tqdm(
            pdb_files, desc="Removing junk ligands", unit="file", total=total_files
        )
    )

    total_ligand_counts = defaultdict(int)
    failed_files = 0

    for result in results:
        if isinstance(result, dict) and "error" not in result:
            for ligand, count in result.items():
                total_ligand_counts[ligand] += count
        else:
            failed_files += 1

    successful_files = total_files - failed_files

    print(f"Total structures processed: {total_files}")
    print(f"Successfully processed: {successful_files}")
    print(f"Failed to process: {failed_files}")

    print("Removed ligands summary:")
    for ligand, count in total_ligand_counts.items():
        print(f"{ligand}: {count} occurrences removed")

    with open("../data/removed_ligands_summary.json", "w") as summary_file:
        json.dump(total_ligand_counts, summary_file, indent=4)
