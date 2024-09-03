import os
from dotenv import load_dotenv

load_dotenv()

BASE_DIR = "/mnt/ligandpro/db/PDB"

RAW_DIR = os.path.join(BASE_DIR, "pdb2/raw/pdb")
PROCESSED_DIR = os.path.join(BASE_DIR, "pdb2/processed")

EMAIL_USER = os.getenv("EMAIL_USER")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
RECEIVER_EMAIL = "Nikolenko.Sergei@icloud.com"

OUTPUT_COMPLEXES_JSON = "../data/complexes.json"
OUTPUT_GROUPED_COMPLEXES_JSON = "../data/grouped_complexes.json"
OUTPUT_SITE_INFO_JSON = "../data/site_info.json"
TRASH_LIGANDS_JSON = "../data/trash_ligands.json"

OUTPUT_LIG = os.path.join(BASE_DIR, "pdb2/lig")
