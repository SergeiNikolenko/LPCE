import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


def calculate_distances(ligand_atoms, protein_atoms):
    distances = []
    for lig_atom in ligand_atoms:
        lig_coord = lig_atom.get_coord()
        min_distance = np.min(
            [
                np.linalg.norm(lig_coord - prot_atom.get_coord())
                for prot_atom in protein_atoms
            ]
        )
        distances.append(min_distance)
    return distances


def is_ligand_buried(pdb_file, threshold, distance_cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)

    ligand_atoms = []
    protein_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0].strip() == "" or residue.id[0].strip() == "W":
                    protein_atoms.extend(residue.get_atoms())
                else:
                    ligand_atoms.extend(residue.get_atoms())

    if not ligand_atoms or not protein_atoms:
        return False, pdb_file.name

    distances = calculate_distances(ligand_atoms, protein_atoms)
    buried_atoms = sum(1 for d in distances if d <= distance_cutoff)
    fraction_buried = buried_atoms / len(ligand_atoms)

    return fraction_buried >= threshold, pdb_file.name


def process_structures(output_dir, threshold=0.3, distance_cutoff=5.0, n_jobs=112):
    pdb_files = list(Path(output_dir).glob("*.pdb"))
    results = Parallel(n_jobs=n_jobs)(
        delayed(is_ligand_buried)(pdb_file, threshold, distance_cutoff)
        for pdb_file in tqdm(pdb_files)
    )

    buried = [name for buried, name in results if buried]
    not_buried = [name for buried, name in results if not buried]

    removed_files = []
    for name in not_buried:
        file_path = Path(output_dir) / name
        file_path.unlink()
        removed_files.append(name)

    logger.info("\n======== Removing not buried ligands ========")
    logger.info(f"Total buried ligands: {len(buried)}")
    logger.info(f"Total not buried ligands: {len(not_buried)}")
    logger.info(
        f"Deleted {len(removed_files)} files, which is {len(removed_files) / len(pdb_files) * 100:.2f}% of the total."
    )

    return {"removed_files": removed_files}


def remove_not_buried_ligands(cfg):
    return process_structures(
        output_dir=cfg.paths.separated_dir,
        threshold=cfg.buried.buried_threshold,
        distance_cutoff=cfg.buried.buried_distance_cutoff,
        n_jobs=cfg.n_jobs,
    )
