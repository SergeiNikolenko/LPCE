import logging
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from tqdm import tqdm

logging.basicConfig(level=logging.ERROR)

sdf_directory = Path("/mnt/ligandpro/db/lpce/ligands/sdf")


def clean_ligand_info(ligand_info):
    return ligand_info


def extract_ligand_and_chain(ligand_info):
    ligand_info = ligand_info.split("/")[-1]
    ligand_info = ligand_info.replace(".pdb", "")
    het, chain = ligand_info.split("_")
    return het, chain


def process_sdf_file(filepath):
    local_data = []
    local_processed = 0
    local_successful = 0

    try:
        supplier = SDMolSupplier(str(filepath), sanitize=False, removeHs=False)
        pdb_id = filepath.stem.replace("pdb", "")
    except Exception:
        return local_data, local_processed, local_successful

    for mol in supplier:
        local_processed += 1

        if mol is None:
            continue

        try:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            local_successful += 1
        except Exception as e:
            smiles = None
            logging.error(f"Error converting molecule in {filepath}: {e}")

        try:
            ligand_info = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        except Exception:
            ligand_info = "unknown"

        if "/" in ligand_info:
            het, chain = extract_ligand_and_chain(ligand_info)
        elif len(ligand_info.split("_")) == 2:
            het, chain = ligand_info.split("_")
            het = clean_ligand_info(het)
            chain = clean_ligand_info(chain)
        else:
            het, chain = "unknown", "unknown"

        local_data.append({"SMILES": smiles, "PDB": pdb_id, "HET": het, "Chain": chain})

    return local_data, local_processed, local_successful


filepaths = list(sdf_directory.glob("*.sdf"))

results = Parallel(n_jobs=-1)(
    delayed(process_sdf_file)(filepath)
    for filepath in tqdm(filepaths, desc="Processing SDF files")
)

data = []
total_processed = 0
total_successful = 0

for result in results:
    data.extend(result[0])
    total_processed += result[1]
    total_successful += result[2]

df = pd.DataFrame(data)
df.to_csv("data/standardized_ligands_sdf.csv", index=False)

success_rate = (total_successful / total_processed) * 100 if total_processed > 0 else 0
print(f"Total molecules processed: {total_processed}")
print(f"Successfully processed molecules: {total_successful}")
print(f"Processing success rate: {success_rate:.2f}%")
