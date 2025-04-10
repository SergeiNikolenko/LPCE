import subprocess
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def add_hydrogens_with_babel(input_file, output_file):
    command = ["obabel", str(input_file), "-O", str(output_file), "-h"]
    subprocess.run(
        command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def add_h_to_ligands(cfg) -> dict:
    if not cfg.add_h_to_ligands:
        logger.info("Skipping hydrogen addition as per configuration")
        return {"status": "skipped", "message": "Hydrogen addition disabled in config"}

    input_path = Path(cfg.paths.separated_dir)
    output_path = Path(cfg.paths.separated_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("\n========== Add hydrogens to ligands ==========")
    pdb_files = list(input_path.glob("*.pdb"))
    logger.info(f"Starting to add hydrogens to {len(pdb_files)} files")

    Parallel(n_jobs=-1)(
        delayed(add_hydrogens_with_babel)(pdb_file, output_path / pdb_file.name)
        for pdb_file in tqdm(pdb_files, desc="Adding hydrogens")
    )

    logger.info("Done!")
    return {"status": "completed", "processed_files": len(pdb_files)}
