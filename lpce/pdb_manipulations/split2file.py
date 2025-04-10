import subprocess
import sys
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


def process_pdb_file(input_file, output_dir, obabel_path):
    """
    Processes a PDB file to separate protein and ligands, convert ligands to SDF and MOL2 formats,
    and save the results into a structured directory.

    Args:
        input_file (Path): Path to the input PDB file.
        output_dir (Path): Path to the output directory.
        obabel_path (str): Path to the obabel executable.
    """
    # Create output folder for the structure
    folder_name = input_file.stem.replace("_processed", "")
    folder_path = output_dir / folder_name
    folder_path.mkdir(parents=True, exist_ok=True)

    # Read the input PDB file
    with input_file.open("r") as f:
        lines = f.readlines()

    header = []
    atom_lines = []
    hetatm_lines = {}

    # Separate header, protein ATOM lines, and ligand HETATM lines
    is_header = True
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            is_header = False
        if is_header:
            header.append(line)
        elif line.startswith("ATOM"):
            atom_lines.append(line)
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            ligand_id = f"{input_file.stem[:4]}_{resname}"
            if ligand_id not in hetatm_lines:
                hetatm_lines[ligand_id] = []
            hetatm_lines[ligand_id].append(line)

    # Save protein without ligands
    protein_file = folder_path / input_file.name.replace("_processed", "")
    with protein_file.open("w") as f:
        f.writelines(header)
        f.writelines(atom_lines)
        f.write("END\n")

    # Process each ligand
    for ligand_id, lines in hetatm_lines.items():
        ligand_pdb_file = folder_path / f"{ligand_id}.pdb"
        # Save ligand in PDB format
        with ligand_pdb_file.open("w") as f:
            f.writelines(header)
            f.writelines(lines)
            f.write("END\n")

        # Define output files
        ligand_sdf_file = folder_path / f"{ligand_id}.sdf"
        ligand_mol2_file = folder_path / f"{ligand_id}.mol2"

        # Convert to SDF and MOL2 using obabel
        subprocess.run(
            [
                obabel_path,
                "-ipdb",
                str(ligand_pdb_file),
                "-osdf",
                "-O",
                str(ligand_sdf_file),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [
                obabel_path,
                "-ipdb",
                str(ligand_pdb_file),
                "-omol2",
                "-O",
                str(ligand_mol2_file),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Remove temporary PDB file
        ligand_pdb_file.unlink()


def create_final_files(cfg):
    """
    Processes all PDB files in a directory, separating proteins and ligands,
    and converting ligands to SDF and MOL2 formats in parallel.

    Args:
        cfg (object): Configuration object with input and output paths.
    """
    input_path = Path(cfg.paths.separated_dir)
    output_path = Path(cfg.paths.final_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("\n========== Create final files ==========")
    pdb_files = list(input_path.glob("*.pdb"))
    logger.info(f"Found {len(pdb_files)} PDB files to process.")

    # Parallel processing of PDB files
    Parallel(n_jobs=cfg.n_jobs)(
        delayed(process_pdb_file)(pdb_file, output_path, cfg.paths.obabel_path)
        for pdb_file in tqdm(pdb_files, desc="Processing PDB files")
    )

    logger.info("Processing completed.")
