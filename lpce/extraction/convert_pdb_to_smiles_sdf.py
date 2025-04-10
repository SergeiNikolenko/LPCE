import subprocess
import tempfile
from pathlib import Path

from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def process_pdb_file(pdb_file_path, smiles_output_dir, sdf_output_dir):
    """
    Processes a single PDB file to extract HETATM blocks, converts them to SMILES and SDF formats,
    and saves the results to the specified output directories.

    Args:
        pdb_file_path (Path): Path to the input PDB file.
        smiles_output_dir (Path): Directory where SMILES files will be saved.
        sdf_output_dir (Path): Directory where SDF files will be saved.

    Returns:
        int: The number of ligands found in the PDB file.
    """
    try:
        with pdb_file_path.open("r") as file:
            lines = file.readlines()

        hetatm_blocks = {}
        ligands_found = 0

        for line in lines:
            if line.startswith("HETATM"):
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                key = f"{res_name}_{chain_id}"
                if key not in hetatm_blocks:
                    hetatm_blocks[key] = []
                hetatm_blocks[key].append(line)
            elif line.startswith("END"):
                break

        if not hetatm_blocks:
            return 0

        ligands_found = len(hetatm_blocks)

        with tempfile.TemporaryDirectory() as hetatm_dir:
            hetatm_dir = Path(hetatm_dir)
            for key, hetatm_lines in hetatm_blocks.items():
                output_file = hetatm_dir / f"{key}.pdb"
                with output_file.open("w") as out_file:
                    out_file.writelines(hetatm_lines)

            smiles_file = smiles_output_dir / pdb_file_path.with_suffix(".smi").name
            sdf_file = sdf_output_dir / pdb_file_path.with_suffix(".sdf").name

            with smiles_file.open("w") as out_smiles:
                for pdb_file in hetatm_dir.glob("*.pdb"):
                    temp_smiles_file = pdb_file.with_suffix(".smi")
                    temp_sdf_file = pdb_file.with_suffix(".sdf")

                    command_smiles = [
                        "obabel",
                        str(pdb_file),
                        "-O",
                        str(temp_smiles_file),
                        "-osmi",
                    ]
                    subprocess.run(command_smiles, capture_output=True, shell=False)

                    command_sdf = [
                        "obabel",
                        str(pdb_file),
                        "-O",
                        str(temp_sdf_file),
                        "-osdf",
                    ]
                    subprocess.run(command_sdf, capture_output=True, shell=False)

                    if temp_smiles_file.exists():
                        with temp_smiles_file.open("r") as temp_file:
                            for line in temp_file:
                                smiles_line = f"{line.strip()}\t{pdb_file.stem}\n"
                                out_smiles.write(smiles_line)

                    if temp_sdf_file.exists():
                        with sdf_file.open("a") as out_sdf:
                            with temp_sdf_file.open("r") as temp_file:
                                out_sdf.write(temp_file.read())

        return ligands_found
    except Exception as e:
        logger.error(f"Failed to process {pdb_file_path}: {e}")
        return 0


def convert_pdb_to_smiles_sdf(cfg: object) -> None:
    """
    Processes all PDB files in the input directory by extracting ligands and converting
    them to SMILES and SDF formats. The results are saved in the output directories.

    Args:
        cfg (object): Configuration object containing paths and other parameters.

    Returns:
        None
    """

    smiles_output_dir = Path(cfg.paths.ligands_dir) / "smiles"
    sdf_output_dir = Path(cfg.paths.ligands_dir) / "sdf"

    smiles_output_dir.mkdir(parents=True, exist_ok=True)
    sdf_output_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = list(Path(cfg.paths.processed_dir).glob("*.pdb"))
    total_files = len(pdb_files)

    logger.info("\n========== Converting PDB to SMILES and SDF ==========")
    logger.info(
        f"Processing {total_files} PDB files from {Path(cfg.paths.processed_dir)}"
    )

    results = Parallel(n_jobs=cfg.n_jobs)(
        delayed(process_pdb_file)(pdb_file, smiles_output_dir, sdf_output_dir)
        for pdb_file in tqdm(pdb_files, desc="Processing PDB files")
    )

    total_ligands = sum(results)

    logger.info(f"Total proteins processed: {total_files}")
    logger.info(f"Total ligands found: {total_ligands}")
