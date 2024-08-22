from dagster import op, get_dagster_logger
import subprocess
from pathlib import Path
import yagmail
import os
from dotenv import load_dotenv
from tqdm import tqdm
import sys
from joblib import Parallel, delayed

load_dotenv()


@op
def run_extract_complexes():
    logger = get_dagster_logger()
    
    def count_structures(directory):
        directory_path = Path(directory)
        return sum(1 for _ in directory_path.rglob("*.ent.gz"))

    def extract_complexes(output_directory):
        logger.info("Starting extraction of ligand-protein complexes from PDB using rsync...")
        
        output_path = Path(output_directory)
        output_path.mkdir(parents=True, exist_ok=True)

        initial_count = count_structures(output_path)
        logger.info(f"Initial count of structures: {initial_count}")
        
        rsync_command = [
            "rsync",
            "-rlPt",
            "--delete",
            "--port=33444",
            "rsync.rcsb.org::ftp_data/structures/divided/pdb/",
            str(output_path)
        ]

        try:
            dry_run_command = rsync_command + ["--dry-run"]

            logger.info("Running dry-run to estimate total files...")
        
            dry_run_process = subprocess.Popen(dry_run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            dry_run_output, dry_run_error = dry_run_process.communicate()

            if dry_run_process.returncode != 0:
                logger.error(f"Dry-run failed with error: {dry_run_error}")
                return 0

            estimated_total_files = len(dry_run_output.strip().split('\n'))

            logger.info(f"Estimated total files to sync: {estimated_total_files}")

            with tqdm(total=estimated_total_files, desc="Syncing files", unit="file", file=sys.stdout) as pbar:
                process = subprocess.Popen(rsync_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                for line in process.stdout:
                    pbar.update(1)
                    logger.debug(line.strip())

                process.wait()

                if process.returncode == 0:
                    final_count = count_structures(output_path)
                    new_structures = final_count - initial_count
                    logger.info(f"Complexes successfully extracted and saved to {output_path}")
                    logger.info(f"Total structures: {final_count}")
                    logger.info(f"New structures added: {new_structures}")
                    return new_structures
                else:
                    logger.error(f"rsync finished with errors, return code: {process.returncode}")
                    logger.error(process.stderr.read())
                    return 0
            
        except Exception as e:
            logger.exception(f"Error during rsync: {e}")
            raise e

    output_directory = "/mnt/ligandpro/db/PDB/pdb2/raw/pdb"
    new_structures = extract_complexes(output_directory)
    return new_structures


@op
def send_email_notification(new_structures: int):
    logger = get_dagster_logger()

    email_user = os.getenv("EMAIL_USER")
    email_password = os.getenv("EMAIL_PASSWORD")

    receiver = "Nikolenko.Sergei@icloud.com"
    subject = "PDB Extraction Job Completed"
    body = f"The job has completed successfully. {new_structures} new ligands were added."

    yag = yagmail.SMTP(email_user, email_password)
    yag.send(to=receiver, subject=subject, contents=body)

    logger.info(f"Email sent to {receiver} from {email_user} with subject '{subject}' and body '{body}'")