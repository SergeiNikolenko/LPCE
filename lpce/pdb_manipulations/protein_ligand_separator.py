import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import NeighborSearch, PDBParser, Select, Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm

warnings.filterwarnings(
    "ignore", category=UserWarning, module="MDAnalysis.core.universe"
)

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


class LigandSelect(Select):
    def __init__(self, ligand_residues, interacting_chains):
        self.ligand_residues = ligand_residues
        self.interacting_chains = interacting_chains
        self.ligand_ids = {
            (residue.get_parent().id, residue.id) for residue in ligand_residues
        }

    def accept_chain(self, chain):
        return chain.id in self.interacting_chains or any(
            chain.id == residue.get_parent().id for residue in self.ligand_residues
        )

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        if (chain_id, residue.id) in self.ligand_ids:
            return True
        elif residue.id[0] == " ":  # Protein residues
            return True
        else:
            return False


def calculate_aligned_rmsd(
    structure1_file, structure2_file, input_file_path, identity_threshold=0.9
):
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("struct1", structure1_file)
    structure2 = parser.get_structure("struct2", structure2_file)

    residues1 = [res for res in structure1.get_residues() if res.id[0] == " "]
    residues2 = [res for res in structure2.get_residues() if res.id[0] == " "]

    seq1_str = "".join([seq1(res.get_resname()) for res in residues1])
    seq2_str = "".join([seq1(res.get_resname()) for res in residues2])

    if not seq1_str or not seq2_str:
        logger.warning(f"Empty sequence in file: {input_file_path}")
        return np.inf

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq1_str, seq2_str)

    if len(alignments) == 0:
        logger.warning(f"Sequence alignment failed. File: {input_file_path}")
        return np.inf

    aligned_seq1, aligned_seq2 = get_aligned_sequences(
        alignments[0], seq1_str, seq2_str
    )

    if aligned_seq1[0] == "-" or aligned_seq2[0] == "-":
        logger.debug(f"Detected shift in aligned sequences for file: {input_file_path}")
        aligned_seq1 = aligned_seq1.lstrip("-")
        aligned_seq2 = aligned_seq2.lstrip("-")

    logger.debug(f"Aligned sequences for file: {input_file_path}")
    logger.debug(f"Sequence 1: {aligned_seq1}")
    logger.debug(f"Sequence 2: {aligned_seq2}")

    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    identity_ratio = matches / max(len(aligned_seq1), len(aligned_seq2))

    if identity_ratio >= identity_threshold:
        logger.debug(
            f"Sequences are {identity_ratio:.2%} identical. Skipping RMSD calculation."
        )
        return 0.0

    idx1, idx2 = 0, 0
    atom_pairs = []
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != "-" and b != "-":
            res1 = residues1[idx1]
            res2 = residues2[idx2]
            for atom_name in ["N", "CA", "C"]:
                if atom_name in res1 and atom_name in res2:
                    atom_pairs.append((res1[atom_name], res2[atom_name]))
        if a != "-":
            idx1 += 1
        if b != "-":
            idx2 += 1

    if len(atom_pairs) == 0:
        logger.warning(
            f"No common atoms for RMSD comparison in file: {input_file_path}"
        )
        return np.inf

    sup = Superimposer()
    sup.set_atoms([pair[0] for pair in atom_pairs], [pair[1] for pair in atom_pairs])
    rmsd_value = sup.rms
    logger.debug(f"RMSD value: {rmsd_value:.2f} for file: {input_file_path}")
    return rmsd_value


def get_aligned_sequences(alignment, seq1, seq2):
    aligned_seq1, aligned_seq2 = [], []
    last_end1, last_end2 = 0, 0

    for (start1, end1), (start2, end2) in zip(
        alignment.aligned[0], alignment.aligned[1]
    ):
        if start1 > last_end1:
            aligned_seq1.append("-" * (start1 - last_end1))
            aligned_seq2.append(seq2[last_end2:start2])
        if start2 > last_end2:
            aligned_seq1.append(seq1[last_end1:start1])
            aligned_seq2.append("-" * (start2 - last_end2))
        aligned_seq1.append(seq1[start1:end1])
        aligned_seq2.append(seq2[start2:end2])
        last_end1 = end1
        last_end2 = end2

    if last_end1 < len(seq1):
        aligned_seq1.append(seq1[last_end1:])
        aligned_seq2.append("-" * (len(seq1) - last_end1))
    if last_end2 < len(seq2):
        aligned_seq1.append("-" * (len(seq2) - last_end2))
        aligned_seq2.append(seq2[last_end2:])

    return "".join(aligned_seq1), "".join(aligned_seq2)


def get_interacting_chains(ligand_atoms, protein_atoms, distance):
    if not protein_atoms:
        logger.debug("No protein atoms found in the structure")
        return set()

    ns = NeighborSearch(protein_atoms)
    interacting_atoms = []
    for atom in ligand_atoms:
        nearby_atoms = ns.search(atom.coord, distance, level="A")
        interacting_atoms.extend(nearby_atoms)
    interacting_chains = {
        atom.get_parent().get_parent().id for atom in interacting_atoms
    }
    logger.debug(f"Interacting chains in file: {interacting_chains}")

    if interacting_chains:
        ligand_residue = ligand_atoms[0].get_parent().get_resname()
        logger.debug(
            f"Ligand {ligand_residue} interacts with chains: {', '.join(interacting_chains)}"
        )
    else:
        logger.debug("No interacting chains found for this ligand group.")

    return interacting_chains


def find_close_ligands(ligand_atoms, all_ligand_atoms, distance):
    ns = NeighborSearch(all_ligand_atoms)
    close_residues = set()
    for atom in ligand_atoms:
        nearby_atoms = ns.search(atom.coord, distance, level="A")
        for nearby_atom in nearby_atoms:
            residue = nearby_atom.get_parent()
            if residue not in ligand_atoms and residue.id[0] != " ":
                close_residues.add(residue)
    return list(close_residues)


def save_structure_to_tempfile(structure, select):
    io = PDBIO()
    temp_file = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    io.set_structure(structure)
    io.save(temp_file.name, select=select)
    temp_file.close()
    return temp_file.name


def process_ligands(structure, interact_distance=4.5, ligand_ligand_distance=3.0):
    all_ligand_atoms = [
        atom for atom in structure.get_atoms() if atom.get_parent().id[0] != " "
    ]
    protein_atoms = [
        atom for atom in structure.get_atoms() if atom.get_parent().id[0] == " "
    ]
    ligand_residues = [
        residue for residue in structure.get_residues() if residue.id[0] != " "
    ]

    if not ligand_residues:
        logger.debug("No ligands found in the structure")
        return []

    processed_ligands = set()
    ligand_groups = []

    for ligand in ligand_residues:
        ligand_id = (ligand.get_resname(), ligand.get_parent().id, ligand.id[1])
        if ligand_id in processed_ligands:
            continue
        ligand_atoms = list(ligand.get_atoms())
        interacting_chains = get_interacting_chains(
            ligand_atoms, protein_atoms, interact_distance
        )
        if not interacting_chains:
            continue
        close_ligands = find_close_ligands(
            ligand_atoms, all_ligand_atoms, ligand_ligand_distance
        )

        all_ligands_in_group = [ligand] + close_ligands

        for close_ligand in close_ligands:
            processed_ligands.add(
                (
                    close_ligand.get_resname(),
                    close_ligand.get_parent().id,
                    close_ligand.id[1],
                )
            )
        processed_ligands.add(ligand_id)
        ligand_groups.append(
            {"ligands": all_ligands_in_group, "interacting_chains": interacting_chains}
        )

    return ligand_groups


def fix_conect_format(conect_lines):
    fixed_conect_lines = []
    for line in conect_lines:
        if line.startswith("CONECT"):
            atom_numbers = line[6:].strip()
            if len(atom_numbers.replace(" ", "")) % 5 == 0:
                atom_numbers_fixed = " ".join(
                    [atom_numbers[i : i + 5] for i in range(0, len(atom_numbers), 5)]
                )
                fixed_line = f"CONECT {atom_numbers_fixed}\n"
            else:
                fixed_line = line
            fixed_conect_lines.append(fixed_line)
        else:
            fixed_conect_lines.append(line)
    return fixed_conect_lines


def filter_conect_lines(conect_lines, saved_atoms):
    """Filter out incorrect atom numbers in CONECT lines, keeping valid ones."""
    filtered_conect_lines = []
    for line in conect_lines:
        if line.startswith("CONECT"):
            atom_numbers = [line[i : i + 5].strip() for i in range(6, len(line), 5)]
            valid_atom_numbers = [
                f"{int(atom_num):5d}"
                for atom_num in atom_numbers
                if atom_num.isdigit() and int(atom_num) in saved_atoms
            ]

            if valid_atom_numbers:
                filtered_conect_lines.append(f"CONECT{''.join(valid_atom_numbers)}\n")
    return filtered_conect_lines


def save_pocket_structure(
    structure,
    pocket_info,
    output_dir,
    saved_structures,
    original_lines,
    conect_lines,
    input_file_path,
    rmsd_threshold=2.0,
    identity_threshold=0.98,
):
    ligands = pocket_info["ligands"]
    interacting_chains = pocket_info["interacting_chains"]
    input_filename = Path(input_file_path).stem

    io = PDBIO()
    io.set_structure(structure)

    select = LigandSelect(ligands, interacting_chains)
    ligand_names = "_".join(sorted({ligand.get_resname() for ligand in ligands}))
    chains_str = "_".join(sorted(interacting_chains))

    output_file = (
        output_dir
        / f"{input_filename}_{ligand_names}_chains_{chains_str}_processed.pdb"
    )

    old_atom_serials = {
        atom.get_serial_number(): atom
        for atom in structure.get_atoms()
        if select.accept_atom(atom)
    }

    temp_structure_file = save_structure_to_tempfile(structure, select)

    similar_found = False
    for saved in saved_structures:
        rmsd_value = calculate_aligned_rmsd(
            saved["temp_file"],
            temp_structure_file,
            input_file_path,
            identity_threshold=identity_threshold,
        )
        if rmsd_value < rmsd_threshold:
            logger.debug(
                f"Found similar structure (RMSD: {rmsd_value:.2f}), skipping file: {output_file}"
            )
            similar_found = True
            break

    if similar_found:
        saved_structures.append(
            {
                "output_file": output_file,
                "temp_file": temp_structure_file,
                "skipped_similar": True,
            }
        )
        return

    parser = PDBParser(QUIET=True)
    new_structure = parser.get_structure("new_struct", temp_structure_file)

    new_atom_serials = {
        atom.get_serial_number(): atom
        for atom in new_structure.get_atoms()
        if select.accept_atom(atom)
    }
    serial_map = {
        old_serial: new_serial
        for (old_serial, old_atom), (new_serial, new_atom) in zip(
            old_atom_serials.items(), new_atom_serials.items()
        )
    }

    updated_conect_lines = []
    for line in conect_lines:
        if line.startswith("CONECT"):
            atom_numbers = [line[i : i + 5].strip() for i in range(6, len(line), 5)]
            new_atom_numbers = []
            for atom_num in atom_numbers:
                if atom_num.isdigit():
                    old_atom_num = int(atom_num)
                    new_atom_num = serial_map.get(old_atom_num, old_atom_num)
                    new_atom_numbers.append(f"{new_atom_num:5d}")
                else:
                    new_atom_numbers.append(atom_num)
            updated_conect_lines.append(f"CONECT{''.join(new_atom_numbers)}\n")

    filter_conect_lines(updated_conect_lines, serial_map.keys())

    with open(output_file, "w") as f_out:
        f_out.writelines(original_lines)
        io.save(f_out, select=select)

    with open(output_file, "r+") as f_out:
        lines = f_out.readlines()
        if lines[-1].strip() == "END":
            lines = lines[:-1]
        f_out.seek(0)
        f_out.writelines(lines)
        f_out.write("END\n")

    saved_structures.append(
        {
            "output_file": output_file,
            "temp_file": temp_structure_file,
            "skipped_similar": False,
        }
    )

    logger.debug(f"Saved pocket structure to file: {output_file}")


def load_original_lines(file_path):
    conect_lines = []
    other_lines = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("CONECT"):
                conect_lines.append(line)
            elif not line.startswith(("ATOM", "HETATM", "END", "MASTER", "TER")):
                other_lines.append(line)
    return other_lines, conect_lines


def process_pdb_file(
    input_file_path,
    output_dir_path,
    interact_distance=4.5,
    ligand_ligand_distance=3.0,
    rmsd_threshold=2.0,
    identity_threshold=0.98,
):
    output_dir = Path(output_dir_path)
    output_dir.mkdir(exist_ok=True)

    original_lines, conect_lines = load_original_lines(input_file_path)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", input_file_path)

    ligand_groups = process_ligands(
        structure, interact_distance, ligand_ligand_distance
    )

    saved_structures = []
    logger.debug(
        f"Processing {len(ligand_groups)} ligand groups in file: {input_file_path}"
    )

    for pocket_info in ligand_groups:
        save_pocket_structure(
            structure,
            pocket_info,
            output_dir,
            saved_structures,
            original_lines,
            conect_lines,
            input_file_path,
            rmsd_threshold,
            identity_threshold,
        )

    return saved_structures


def analyze_protein(
    input_file_path,
    output_dir_path="separated_complexes",
    interact_distance=4.5,
    ligand_ligand_distance=3.5,
    rmsd_threshold=2.0,
    identity_threshold=0.98,
):
    logger.debug(f"Starting analysis of protein: {input_file_path}")

    saved_structures = process_pdb_file(
        input_file_path,
        output_dir_path,
        interact_distance,
        ligand_ligand_distance,
        rmsd_threshold,
    )

    saved_count = sum(
        1 for s in saved_structures if not s.get("skipped_similar", False)
    )
    skipped_similar = sum(
        1 for s in saved_structures if s.get("skipped_similar", False)
    )

    logger.debug("Protein analysis completed")

    return {
        "structures_saved": saved_count,
        "structures_skipped_similar": skipped_similar,
    }


def get_pdb_files(input_dir):
    pdb_files = list(Path(input_dir).glob("*.pdb"))
    return pdb_files


def protein_ligand_separator(cfg):
    input_dir = Path(cfg.paths.bioml_dir)
    output_dir = Path(cfg.paths.separated_dir)

    interact_distance = cfg.separator_params.interact_distance
    ligand_ligand_distance = cfg.separator_params.ligand_ligand_distance
    rmsd_threshold = cfg.separator_params.rmsd_threshold
    identity_threshold = cfg.separator_params.identity_threshold
    pdb_files = get_pdb_files(input_dir)
    total_files = len(pdb_files)

    os.makedirs(output_dir, exist_ok=True)
    logger.info("\n========== Protein ligand separator completed ==========")
    results = Parallel(n_jobs=cfg.n_jobs)(
        delayed(analyze_protein)(
            pdb_file,
            output_dir,
            interact_distance,
            ligand_ligand_distance,
            rmsd_threshold,
            identity_threshold,
        )
        for pdb_file in tqdm(pdb_files, desc="Separating ligand pockets in PDB files")
    )

    total_structures_saved = sum(result["structures_saved"] for result in results)
    total_structures_skipped_similar = sum(
        result["structures_skipped_similar"] for result in results
    )
    total_structures_left = total_structures_saved

    logger.info(f"Total PDB files found: {total_files}")

    logger.info(f"Total similar structures skipped: {total_structures_skipped_similar}")
    logger.info(f"Structures remaining after filtering: {total_structures_left}")
    logger.info(f"Total structures SAVED: {total_structures_saved}")
