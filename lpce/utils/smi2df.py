# %%
import json
import warnings
from pathlib import Path

import datamol as dm
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")

smi_directory = Path("/mnt/ligandpro/db/PDB/pdb2/lig/smi_files")

data = []
for filepath in smi_directory.glob("*.smi"):
    pdb_id = filepath.stem.replace("pdb", "")
    with filepath.open("r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                smiles, _, ligand_chain = parts
                het, chain = ligand_chain.split("_")
                data.append(
                    {"SMILES": smiles, "PDB": pdb_id, "HET": het, "Chain": chain}
                )

df = pd.DataFrame(data)


with Path("data/trash_ligands.json").open("r") as f:
    trash_ligands = json.load(f)

initial_count = len(df)
df_filtered = df[~df["HET"].isin(trash_ligands)]
final_count = len(df_filtered)
removed_count = initial_count - final_count

print(f"Total records before filtering: {initial_count}")
print(f"Total records after filtering: {final_count}")
print(f"Total records removed: {removed_count}")

smiles_column = "SMILES"


def _preprocess(i, row):
    try:
        dm.disable_rdkit_log()
        mol = dm.to_mol(row[smiles_column], ordered=True)
        if mol is None:
            return None

        mol = dm.fix_mol(mol)
        mol = dm.sanitize_mol(mol, sanifix=True, charge_neutral=False)
        mol = dm.standardize_mol(
            mol,
            disconnect_metals=False,
            normalize=True,
            reionize=True,
            uncharge=True,
            stereo=True,
        )

        row["standard_smiles"] = dm.standardize_smiles(dm.to_smiles(mol))
        return row

    except Exception:
        return None


results = Parallel(n_jobs=-1)(
    delayed(_preprocess)(i, row)
    for i, row in tqdm(df_filtered.iterrows(), total=len(df_filtered))
)

df_standardized = pd.DataFrame([res for res in results if res is not None])

successful_conversions = len(df_standardized)
total_molecules = len(df_filtered)
conversion_rate = (successful_conversions / total_molecules) * 100

print(f"Total molecules: {total_molecules}")
print(f"Successfully standardized molecules: {successful_conversions}")
print(f"Conversion success rate: {conversion_rate:.2f}%")

df_standardized.to_csv("data/standardized_ligands.csv", index=False)
