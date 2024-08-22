import os
from dotenv import load_dotenv

load_dotenv()

BASE_DIR = "/mnt/ligandpro/db/PDB"

RAW_DIR = os.path.join(BASE_DIR, "pdb2/raw/pdb")
PROCESSED_DIR = os.path.join(BASE_DIR, "pdb2/processed")

EMAIL_USER = os.getenv("EMAIL_USER")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
RECEIVER_EMAIL = "Nikolenko.Sergei@icloud.com"
