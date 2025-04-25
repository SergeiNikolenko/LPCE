import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from loguru import logger


def parse_ligands_from_pdb(pdb_file: Path):
    ligands = {}
    with pdb_file.open() as f:
        for line in f:
            if line.startswith("HETATM"):
                resn = line[17:20].strip()
                if resn == "HOH":
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligands.setdefault(resn, []).append((x, y, z))
    return ligands


def compute_rmsd(coords1, coords2):
    if len(coords1) != len(coords2):
        return float("inf")
    arr1 = np.array(coords1)
    arr2 = np.array(coords2)
    diff = arr1 - arr2
    return math.sqrt(np.mean(np.sum(diff * diff, axis=1)))


def are_ligands_same(pdbA: Path, pdbB: Path, ligand_rmsd_thr=1.0):
    a = parse_ligands_from_pdb(pdbA)
    b = parse_ligands_from_pdb(pdbB)
    if set(a.keys()) != set(b.keys()):
        return False
    for lig_name in a:
        coordsA = a[lig_name]
        coordsB = b[lig_name]
        rmsd_val = compute_rmsd(coordsA, coordsB)
        if rmsd_val > ligand_rmsd_thr:
            return False
    return True


def refine_groups_by_ligands(identical_groups, pdb_dir: Path, ligand_rmsd_thr=1.0):
    refined = []
    for group in identical_groups:
        items = sorted(group)
        subgroups = []
        for pdb_name in items:
            placed = False
            pdb_file = pdb_dir / f"{pdb_name}.pdb"
            if not pdb_file.exists():
                continue
            for sg in subgroups:
                rep = next(iter(sg))
                rep_file = pdb_dir / f"{rep}.pdb"
                if are_ligands_same(pdb_file, rep_file, ligand_rmsd_thr):
                    sg.add(pdb_name)
                    placed = True
                    break
            if not placed:
                subgroups.append({pdb_name})
        refined.extend(subgroups)
    return refined


def extract_het_and_chain_identifiers(filename: str):
    chain_tag = "_chains_"
    processed_tag = "_processed"
    try:
        start_bioml = filename.index("bioml_") + len("bioml_")
        pos_underscore = filename.index("_", start_bioml)
        start_ligand = pos_underscore + 1
        chain_start = filename.index(chain_tag, start_ligand)
        ligand_part = filename[start_ligand:chain_start]
        chain_identifier_start = chain_start + len(chain_tag)
        processed_pos = filename.index(processed_tag, chain_identifier_start)
        chain_identifier = filename[chain_identifier_start:processed_pos]
        return (ligand_part, chain_identifier)
    except ValueError:
        return None, None


def get_resolution_from_pdb(pdb_file: Path):
    try:
        with pdb_file.open() as f:
            for line in f:
                if line.startswith("REMARK   2 RESOLUTION."):
                    parts = line.strip().split()
                    for i, part in enumerate(parts):
                        if "ANGSTROM" in part.upper():
                            try:
                                return float(parts[i - 1])
                            except ValueError:
                                return None
                    val = line.strip().split("RESOLUTION.")[-1].split()[0]
                    try:
                        return float(val)
                    except ValueError:
                        return None
        return None
    except Exception as e:
        logger.error(f"Error reading resolution from {pdb_file}: {e}")
        return None


def remove_duplicate_groups(groups):
    unique_groups = []
    groups_sorted = sorted([sorted(g) for g in groups], key=lambda x: (-len(x), x))
    for group in groups_sorted:
        group_set = set(group)
        is_subset = False
        for uniq in unique_groups:
            if group_set <= uniq:
                is_subset = True
                break
        if not is_subset:
            unique_groups.append(group_set)
    return unique_groups


def process_groups_with_resolution(identical_groups, pdb_dir: Path):
    processed_groups = []
    for group in identical_groups:
        short = {f.split("_processed")[0] + "_processed" for f in group}
        processed_groups.append(short)
    unique_groups = remove_duplicate_groups(processed_groups)
    final_groups = []
    for group in unique_groups:
        hc_groups = defaultdict(list)
        for file_name in group:
            hetID, chainID = extract_het_and_chain_identifiers(file_name)
            hc_groups[(hetID, chainID)].append(file_name)
        chosen = []
        for _, items in hc_groups.items():
            best_struc = None
            best_resol = None
            for structure in items:
                pdb_file = pdb_dir / f"{structure}.pdb"
                if not pdb_file.exists():
                    logger.debug(f"File {pdb_file} not found.")
                    continue
                reso = get_resolution_from_pdb(pdb_file)
                if reso is None:
                    logger.debug(f"Resolution not found for {structure}")
                    continue
                if best_resol is None or reso < best_resol:
                    best_resol = reso
                    best_struc = structure
            if best_struc:
                chosen.append(best_struc)
            else:
                chosen.append(items[0])
        final_groups.append(set(chosen))
    return final_groups


def remove_similar_structures(cfg):
    logger.info("\n========== Remove similar structures ==========")

    with open(cfg.foldseek.identical_groups) as f:
        data = json.load(f)
    identical_groups = [set(g) for g in data]

    pdb_dir = Path(cfg.paths.separated_dir)
    all_files = {p.stem for p in pdb_dir.glob("*.pdb")}

    refined = refine_groups_by_ligands(identical_groups, pdb_dir, ligand_rmsd_thr=1.0)
    final_groups = process_groups_with_resolution(refined, pdb_dir)

    in_groups = set()
    for grp in final_groups:
        in_groups |= {f"{x}.pdb" for x in grp}

    # Anything that never appeared in Foldseek is also kept
    never_foldseek = {f"{x}.pdb" for x in (all_files - set().union(*identical_groups))}
    keep_files = in_groups | never_foldseek
    initial_count = len(all_files)
    deleted_files = 0

    for item in pdb_dir.glob("*.pdb"):
        if item.name not in keep_files:
            item.unlink()
            deleted_files += 1

    final_count = initial_count - deleted_files
    perc = (deleted_files / initial_count) * 100
    logger.info(f"Initially {initial_count} PDB.")
    logger.info(f"{final_count} remain.")
    logger.info(f"Deleted {deleted_files} ({perc:.2f}%).")
