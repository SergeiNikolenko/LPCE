import subprocess
from pathlib import Path
from tqdm import tqdm
import sys
from config.settings import RAW_DIR

def count_structures(directory):
    directory_path = Path(directory)
    return sum(1 for _ in directory_path.rglob("*.ent.gz"))

def extract_complexes():
    output_path = Path(RAW_DIR)
    output_path.mkdir(parents=True, exist_ok=True)

    initial_count = count_structures(output_path)
    print(f"Initial count of structures: {initial_count}")
    
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

        print("Running dry-run to estimate total files...")
    
        dry_run_process = subprocess.Popen(dry_run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        dry_run_output, dry_run_error = dry_run_process.communicate()

        if dry_run_process.returncode != 0:
            print(f"Dry-run failed with error: {dry_run_error}")
            return 0

        estimated_total_files = len(dry_run_output.strip().split('\n'))

        print(f"Estimated total files to sync: {estimated_total_files}")

        with tqdm(total=estimated_total_files, desc="Syncing files", unit="file", file=sys.stdout) as pbar:
            process = subprocess.Popen(rsync_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            for line in process.stdout:
                pbar.update(1)
                print(line.strip())

            process.wait()

            if process.returncode == 0:
                final_count = count_structures(output_path)
                new_structures = final_count - initial_count
                print(f"Complexes successfully extracted and saved to {output_path}")
                print(f"Total structures: {final_count}")
                print(f"New structures added: {new_structures}")
                return new_structures
            else:
                print(f"rsync finished with errors, return code: {process.returncode}")
                print(process.stderr.read())
                return 0
        
    except Exception as e:
        print(f"Error during rsync: {e}")
        raise e
