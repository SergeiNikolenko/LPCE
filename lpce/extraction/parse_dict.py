from pathlib import Path
from Bio import PDB
from joblib import Parallel, delayed
import json
from tqdm import tqdm
import re
import string
from config.settings import PROCESSED_DIR, OUTPUT_COMPLEXES_JSON, OUTPUT_GROUPED_COMPLEXES_JSON, OUTPUT_SITE_INFO_JSON, TRASH_LIGANDS_JSON

warnings = []

def parse_pdb_file(file_path: Path) -> list:
    """
    Parses a PDB file to extract ligand information.

    Args:
        file_path (Path): Path to the PDB file to be parsed.

    Returns:
        list: A list of tuples containing ligand information (ligand_id, res_number, chain_id).
        Returns None if the file is empty or if an error occurs during parsing.
    """
    if file_path.stat().st_size == 0:
        warnings.append(f"Warning: {file_path} is empty and will be skipped.")
        return None

    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('', file_path)
    except (ValueError, AttributeError) as e:
        warnings.append(f"Error parsing {file_path}: {e}")
        return None

    has_protein = False
    ligands = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    has_protein = True
                elif residue.id[0] != " ":
                    ligand_id = residue.resname
                    res_number = residue.id[1]
                    chain_id = chain.id
                    ligands.append((ligand_id, res_number, chain_id))

    if has_protein and ligands:
        return ligands
    else:
        return None

def process_single_file(pdb_file: Path, pdb_directory: Path) -> tuple:
    """
    Processes a single PDB file to extract ligand information.

    Args:
        pdb_file (Path): The PDB file to process.
        pdb_directory (Path): The directory containing the PDB files.

    Returns:
        tuple: A tuple containing the PDB ID and the extracted ligand information.
    """
    pdb_id = pdb_file.stem[3:]
    file_path = pdb_directory / pdb_file
    return pdb_id, parse_pdb_file(file_path)

def extract_complexes_with_ligands(pdb_directory: Path, max_files: int = None) -> dict:
    """
    Extracts complexes with ligands from PDB files in a directory.

    Args:
        pdb_directory (Path): The directory containing the PDB files.
        max_files (int, optional): Maximum number of files to process. Defaults to None.

    Returns:
        dict: A dictionary with PDB IDs as keys and ligand information as values.
    """
    pdb_files = list(pdb_directory.glob("*.pdb"))

    if max_files is not None:
        pdb_files = pdb_files[:max_files]

    print(f"Starting to process {len(pdb_files)} PDB files...")

    results = []
    for result in Parallel(n_jobs=56)(
            delayed(process_single_file)(pdb_file, pdb_directory) 
            for pdb_file in tqdm(pdb_files, desc="Processing files", unit="file")):
        results.append(result)

    pdb_complexes = {pdb_id: ligands for pdb_id, ligands in results if ligands is not None}
    
    print(f"Completed processing of {len(pdb_files)} PDB files.")
    return pdb_complexes

def save_complexes_to_json(complexes: dict, output_path: Path) -> None:
    """
    Saves complexes to a JSON file.

    Args:
        complexes (dict): The complexes to save.
        output_path (Path): The path to the output JSON file.
    """
    with open(output_path, 'w') as f:
        json.dump(complexes, f, indent=4)
    print(f"Complexes saved to {output_path}")

def create_grouped_complexes_dict(complexes: dict) -> dict:
    """
    Groups complexes by chain.

    Args:
        complexes (dict): The complexes to group.

    Returns:
        dict: A dictionary where complexes are grouped by chain.
    """
    new_complexes_dict = {}
    for pdb_id, ligands in complexes.items():
        grouped_by_chain = {}

        for ligand in ligands:
            ligand_id, res_number, chain_id = ligand

            if chain_id not in grouped_by_chain:
                grouped_by_chain[chain_id] = []

            ligand_entry = {"ligand": ligand_id, "residue": res_number}
            if ligand_entry not in grouped_by_chain[chain_id]:
                grouped_by_chain[chain_id].append(ligand_entry)

        new_complexes_dict[pdb_id] = grouped_by_chain

    return new_complexes_dict

def extract_site_info(pdb_file: Path) -> list:
    """
    Extracts site information from a PDB file.

    Args:
        pdb_file (Path): The PDB file from which to extract site information.

    Returns:
        list: A list of site information lines.
    """
    site_info = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("SITE"):
                site_info.append(line.strip())
    return site_info

def get_pdb_id(pdb_file_name: str) -> str:
    """
    Extracts the PDB ID from a PDB file name.

    Args:
        pdb_file_name (str): The name of the PDB file.

    Returns:
        str: The extracted PDB ID.
    """
    return pdb_file_name.replace(".pdb", "")[3:]

def load_trash_ligands(file_path: Path) -> set:
    """
    Loads trash ligands from a JSON file.

    Args:
        file_path (Path): The path to the JSON file containing trash ligands.

    Returns:
        set: A set of trash ligands.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        return set(json.load(file).keys())

def process_site_info_from_pdb_files(pdb_directory: Path) -> dict:
    """
    Processes site information from PDB files in a directory.

    Args:
        pdb_directory (Path): The directory containing the PDB files.

    Returns:
        dict: A dictionary with PDB IDs as keys and site information as values.
    """
    pdb_files = list(pdb_directory.glob('*.pdb'))
    site_info_list = Parallel(n_jobs=-1)(delayed(extract_site_info)(pdb_file) for pdb_file in tqdm(pdb_files))
    return {get_pdb_id(pdb_file.name): site_info for pdb_file, site_info in zip(pdb_files, site_info_list) if site_info}

def clean_site_info(site_info: list) -> list:
    """
    Cleans site information by removing unwanted characters and formatting it.

    Args:
        site_info (list): The site information to clean.

    Returns:
        list: The cleaned site information.
    """
    cleaned_site_info = []
    
    for site in site_info:
        site = re.sub(r'SITE|\b\d+\b', '', site).strip()
        site = site.replace(',', '')
        site = re.sub(r'\b([A-Z]{2})\b(?!\s[A-Z])', '', site).strip()
        site = re.sub(r'\s+', ' ', site)
        cleaned_site_info.append(site)
    
    return cleaned_site_info

def transform_site_info(site_info_dict: dict, trash_ligands: set) -> dict:
    """
    Transforms site information by excluding unwanted ligands and amino acids.

    Args:
        site_info_dict (dict): The dictionary containing site information.
        trash_ligands (set): The set of trash ligands to exclude.

    Returns:
        dict: The transformed site information.
    """
    amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}

    to_exclude = amino_acids.union(trash_ligands)
    
    remove_digits = str.maketrans('', '', string.digits)

    transformed_dict = {}

    for key, entries in site_info_dict.items():
        ligand_dict = {}
        for entry in entries:
            parts = entry.split()
            residues = parts[1:]

            for i in range(0, len(residues) - 1, 2):
                residue = residues[i]
                chain = residues[i+1]

                if residue in to_exclude or len(residue) != 3:
                    continue

                chain = chain.translate(remove_digits)

                if residue not in ligand_dict:
                    ligand_dict[residue] = []

                if chain not in ligand_dict[residue]:
                    ligand_dict[residue].append(chain)

        if ligand_dict:
            transformed_dict[key] = ligand_dict

    return transformed_dict

def print_warnings() -> None:
    """
    Prints any warnings that were collected during the execution of the program.
    """
    if warnings:
        print("\nWarnings:")
        for warning in warnings:
            print(warning)
    else:
        print("\nNo warnings.")

def extract_and_save_complexes_with_ligands(max_files: int = None) -> None:
    """
    Extracts complexes with ligands from PDB files, processes and saves them to JSON files.

    Args:
        max_files (int, optional): Maximum number of files to process. Defaults to None.
    """
    pdb_directory = Path(PROCESSED_DIR)

    complexes = extract_complexes_with_ligands(pdb_directory, max_files)
    save_complexes_to_json(complexes, Path(OUTPUT_COMPLEXES_JSON))
    
    grouped_complexes = create_grouped_complexes_dict(complexes)
    save_complexes_to_json(grouped_complexes, Path(OUTPUT_GROUPED_COMPLEXES_JSON))

    site_info_dict = process_site_info_from_pdb_files(pdb_directory)
    final_cleaned_site_info_dict = {pdb_id: clean_site_info(site_info) for pdb_id, site_info in site_info_dict.items()}
    trash_ligands = load_trash_ligands(Path(TRASH_LIGANDS_JSON))
    final_transformed_site_info = transform_site_info(final_cleaned_site_info_dict, trash_ligands)
    save_complexes_to_json(final_transformed_site_info, Path(OUTPUT_SITE_INFO_JSON))

    print_warnings()

if __name__ == "__main__":
    extract_and_save_complexes_with_ligands()
