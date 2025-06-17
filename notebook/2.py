import sys
import re
import os
from pathlib import Path
import tempfile
import threading
import uuid
import warnings
from dataclasses import dataclass

import logging
import numpy as np
import pymol
from Bio.PDB import PDBIO, NeighborSearch, PDBParser, Select
from loguru import logger
from pymol import cmd
from spyrmsd import graph, molecule
import datamol as dm
from rdkit import RDLogger
from collections import Counter
from typing import FrozenSet
from rdkit import Chem

for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(var, "1")

class _SilentStream:
    _pat = re.compile(
        r"(Undo has been disabled|Reason: Memory exceeded|"
        r"Failed to find the pandas get_adjustment\(\) function to patch|"
        r"Failed to patch pandas|Matrix: Warning: no convergence|"
        r"A worker stopped while some jobs were given to the executor)"
    )
    def __init__(self, stream): self._s = stream
    def write(self, txt):
        if not self._pat.search(txt): self._s.write(txt)
    def flush(self): self._s.flush()

sys.stdout = _SilentStream(sys.stdout)
sys.stderr = _SilentStream(sys.stderr)

dm.disable_rdkit_log()

logging.getLogger("rdkit").setLevel(logging.ERROR)
logging.getLogger("rdkit.Chem.PandasTools").setLevel(logging.CRITICAL)
RDLogger.DisableLog("rdApp.*")

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.core.universe")
warnings.filterwarnings("ignore", message="Disconnected graph detected.*", module="spyrmsd.graphs.rx")

logger.remove()
logger.add(sys.stdout, format="{message}", level="DEBUG")


BOND_LENGTHS = {
    ("C", "C"): 1.54,
    ("C", "N"): 1.65,
    ("C", "O"): 1.63,
    ("C", "S"): 1.80,
    ("P", "O"): 1.64,
    ("N", "O"): 1.55,
    ("N", "S"): 1.70,
    ("Fe", "S"): 2.33,
    ("Fe", "O"): 2.20,
    ("Fe", "N"): 2.15,
    ("Zn", "S"): 2.39,
    ("Zn", "N"): 2.26,
    ("Zn", "O"): 2.42,
    ("Cu", "S"): 2.26,
    ("Cu", "N"): 1.98,
    ("Cu", "O"): 1.95,
    ("Ni", "S"): 2.28,
    ("Ni", "N"): 2.18,
    ("Ni", "O"): 2.20,
    ("Mn", "O"): 2.18,
    ("Mn", "N"): 2.20,
    ("Co", "S"): 2.32,
    ("Co", "N"): 2.25,
    ("Co", "O"): 2.20,
    ("Mg", "O"): 2.18,
    ("Ca", "O"): 2.45,
    ("Mo", "S"): 2.42,
    ("Mo", "O"): 2.00,
    ("Na", "O"): 2.35,
    ("K", "O"): 2.75,
}

METALS = {"Fe", "Zn", "Cu", "Ni", "Mn", "Co", "Mg", "Ca", "Mo", "Na", "K"}

def _get_bond_threshold(el1, el2, default_bond_distance):
    pair = (el1, el2)
    if pair in BOND_LENGTHS:
        return BOND_LENGTHS[pair]
    rev = (el2, el1)
    if rev in BOND_LENGTHS:
        return BOND_LENGTHS[rev]
    return default_bond_distance

def _elem(atom) -> str:
    e = getattr(atom, "element", "").strip()
    if e:
        return e.capitalize()
    name = atom.get_name().strip()
    m = re.match(r"([A-Za-z]{1,2})", name)
    if not m:
        return name[0].upper()
    s = m.group(1).upper()
    return s.capitalize()

class _PymolSession:
    _lock = threading.Lock()
    _started = False

    @classmethod
    def ensure(cls) -> None:
        with cls._lock:
            if not cls._started:
                pymol.finish_launching(["pymol", "-qc"])
                cmd.feedback("disable", "all", "everything")
                cls._started = True

    @classmethod
    def locked(cls):
        cls.ensure()
        return cls._lock

def _ligands_connected(r1, r2, default_bond_distance):
    for a1 in r1.get_atoms():
        for a2 in r2.get_atoms():
            el1 = _elem(a1)
            el2 = _elem(a2)
            thr = _get_bond_threshold(el1, el2, default_bond_distance)
            dist = np.linalg.norm(a1.coord - a2.coord)
            # logger.debug(f"_ligands_connected: Checking {r1.get_resname()}-{r2.get_resname()} : "
            #             f"{a1.get_name()}-{a2.get_name()}, dist={dist:.3f}, thr={thr}")
            if dist <= thr:
                # logger.debug(f"_ligands_connected: FOUND BOND between {r1.get_resname()} and {r2.get_resname()}")
                return True
    return False

def _rmsd_on_ca(pdb1: str, pdb2: str) -> float:
    with _PymolSession.locked():
        o1, o2 = f"o{uuid.uuid4().hex}", f"o{uuid.uuid4().hex}"
        cmd.delete("all")
        cmd.load(pdb1, o1)
        cmd.load(pdb2, o2)
        try:
            return float(cmd.align(f"{o1} and name CA", f"{o2} and name CA")[0])
        finally:
            cmd.delete("all")
            cmd.reinitialize()

def _rmsd_on_ligand(pdb1: str, pdb2: str) -> float:
    with _PymolSession.locked():
        o1, o2 = f"o{uuid.uuid4().hex}", f"o{uuid.uuid4().hex}"
        cmd.delete("all")
        cmd.load(pdb1, o1)
        cmd.load(pdb2, o2)

        cmd.align(f"{o1} and name CA", f"{o2} and name CA", cycles=0)

        def _dump(obj, fn):
            tmp_obj = f"{obj}_x"
            cmd.create(tmp_obj, obj, 1, 1)
            cmd.save(fn, f"{tmp_obj} and hetatm and not resn HOH", state=1)
            cmd.delete(tmp_obj)

        t1 = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        t2 = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        _dump(o1, t1.name)
        _dump(o2, t2.name)

        m1 = Chem.MolFromPDBFile(t1.name, removeHs=False)
        m2 = Chem.MolFromPDBFile(t2.name, removeHs=False)
        t1.close()
        t2.close()
        if not m1 or not m2:
            logger.debug("    ligand RMSD -> RDKit build failed (m1 or m2 is None)")
            cmd.delete("all")
            return float("inf")

        c1 = np.array([a.coord for a in cmd.get_model(o1 + " and hetatm and not resn HOH").atom])
        c2 = np.array([a.coord for a in cmd.get_model(o2 + " and hetatm and not resn HOH").atom])

        try:
            mol1 = molecule.Molecule.from_rdkit(m1)
            mol2 = molecule.Molecule.from_rdkit(m2)
            G1 = graph.graph_from_adjacency_matrix(mol1.adjacency_matrix, mol1.atomicnums)
            G2 = graph.graph_from_adjacency_matrix(mol2.adjacency_matrix, mol2.atomicnums)


            for idx1, idx2 in graph.match_graphs(G1, G2):
                r = np.sqrt(np.mean((c1[idx1] - c2[idx2]) ** 2))
                logger.debug("    ligand RMSD -> graph mode matched. RMSD=%.3f" % r)
                cmd.delete("all")
                return float(r)
        except Exception as e:
            logger.debug(f"    ligand RMSD -> graph error: {e}")

        logger.debug("    ligand RMSD -> return inf (no match)")
        cmd.delete("all")
        return float("inf")


@dataclass
class Pocket:
    ligands: List
    chains: Set[str]
    bond_type: str


class LigandPocketExtractor:

    def __init__(self, interact_d=4.5, ligand_cluster_d=1.5, default_bond_d=1.5, short_peptide=10):
        self.interact_d = interact_d
        self.default_bond_d = default_bond_d
        self.short_peptide = short_peptide
        self.ligand_cluster_d = max(ligand_cluster_d, default_bond_d)
        self._search_radius = max(self.ligand_cluster_d, max(BOND_LENGTHS.values(), default=0.0))

    @staticmethod

    def _chains_near(lig_atoms, prot_atoms, d):
        if not prot_atoms:
            return set()
        ns = NeighborSearch(prot_atoms)
        return {
            a.get_parent().get_parent().id
            for la in lig_atoms
            for a in ns.search(la.coord, d, level="A")
        }

    @staticmethod
    def _save_temp(structure, selector):
        io = PDBIO()
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        io.set_structure(structure)
        io.save(tmp.name, select=selector)
        tmp.close()
        return tmp.name

    def _check_bond_type(self, la, pa):
        el_l = _elem(la)
        el_p = _elem(pa)
        thr = _get_bond_threshold(el_l, el_p, self.default_bond_d)
        dist = np.linalg.norm(la.coord - pa.coord)
        if dist <= thr:
            if el_l in METALS or el_p in METALS:
                return "coord"
            return "cov"
        return None

    def _get_pocket_bond_type(self, lig_atoms, prot_atoms):
        if not prot_atoms:
            return ""

        p_ns = NeighborSearch(prot_atoms)
        all_ns = NeighborSearch(prot_atoms + lig_atoms)

        for la in lig_atoms:
            for pa in all_ns.search(la.coord, 3.0, level="A"):
                if pa is la:
                    continue
                if _elem(la) in METALS or _elem(pa) in METALS:
                    bt = self._check_bond_type(la, pa)
                    if bt == "coord":
                        return "coord"

        for la in lig_atoms:
            for pa in p_ns.search(la.coord, 3.0, level="A"):
                bt = self._check_bond_type(la, pa)
                if bt == "cov":
                    return "cov"

        return ""

    def extract(self, structure, mod_positions: FrozenSet[tuple[str, int, str]]):
        if len(structure) > 1:
            logger.debug("extract: multiple models, using first model.")
            structure = structure[0]

        chain_res = {}
        for r in structure.get_residues():
            chain_res.setdefault(r.get_parent().id, []).append(r)

        protein_chains = {
            cid
            for cid, res in chain_res.items()
            if sum(1 for r in res if r.id[0] == " ") > self.short_peptide
        }

        prot_atoms = []
        for a in structure.get_atoms():
            r = a.get_parent()
            ch = r.get_parent().id
            if ch in protein_chains and (
                r.id[0] == " " or (ch, r.id[1], r.id[2]) in mod_positions
            ):
                prot_atoms.append(a)

        for cid, reslist in chain_res.items():
            logger.debug(f"Chain {cid} has {len(reslist)} residues. Num_std={sum(1 for r in reslist if r.id[0]==' ')}")

        candidates = []
        for r in structure.get_residues():
            ch, num, icode = r.get_parent().id, r.id[1], r.id[2]
            if (ch, num, icode) in mod_positions:
                continue
            if (r.id[0] != " " and r.get_resname() != "HOH") or \
               len(chain_res[ch]) <= self.short_peptide:
                candidates.append(r)

        logger.debug(f"Candidate residues: {len(candidates)}")

        if not candidates:
            logger.debug("extract: no candidate residues found — skip")
            return []

        index = {r: i for i, r in enumerate(candidates)}
        parent = list(range(len(candidates)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            pi, pj = find(i), find(j)
            if pi != pj:
                parent[pj] = pi

        ns = NeighborSearch([a for r in candidates for a in r.get_atoms()])
        for r in candidates:
            i = index[r]
            for a in r.get_atoms():
                for n in ns.search(a.coord, self._search_radius, level="A"):
                    r2 = n.get_parent()
                    if r2 is r or r2 not in index:
                        continue
                    dist = np.linalg.norm(a.coord - n.coord)
                    thr = _get_bond_threshold(_elem(a), _elem(n), self.ligand_cluster_d)
                    if dist <= thr:
                        union(i, index[r2])

        comps = {}
        for r in candidates:
            comps.setdefault(find(index[r]), []).append(r)

        pockets = []
        for residues in comps.values():
            lig_atoms = [a for rr in residues for a in rr.get_atoms()]
            chains = self._chains_near(lig_atoms, prot_atoms, self.interact_d)
            if not chains:
                continue
            btype = self._get_pocket_bond_type(lig_atoms, prot_atoms)
            logger.debug(f"New pocket with residues {[r.get_resname() for r in residues]} -> bond={btype}, chains={chains}")
            pockets.append(Pocket(residues, chains, btype))

        logger.debug(f"Pockets found: {len(pockets)}")
        return pockets

import copy

def _unique_name(base, taken):
    import string
    base = base.strip().upper()[:2]
    pool = string.ascii_uppercase + string.digits
    for s in pool:
        cand = (base + s).ljust(4)[:4]
        if cand not in taken:
            return cand
    for s1 in pool:
        for s2 in pool:
            cand = (base + s1 + s2)[:4]
            if cand not in taken:
                return cand
    raise ValueError("too many duplicates")

class _MergedSelect(Select):
    def __init__(self, ligands, chains):
        self._ligands = set(ligands)
        self._chains = chains
    def accept_chain(self, chain):
        return chain.id in self._chains or any(r.get_parent() is chain for r in self._ligands)
    def accept_residue(self, residue):
        if residue in self._ligands:
            return True
        return residue.get_parent().id in self._chains and residue.id[0] == " "

class PocketWriter:
    def __init__(
        self,
        rmsd_thr=2.0,
        lig_rmsd_thr=0.5,
        ligand_cluster_distance=3.0,
        default_bond_distance=1.5,
    ):
        self.rmsd_thr = rmsd_thr
        self.lig_rmsd_thr = lig_rmsd_thr
        self.ligand_cluster_distance = ligand_cluster_distance
        self.default_bond_distance = default_bond_distance
        self._saved = []
        self._skipped = 0

    def _is_duplicate(self, new_tmp, lig_resnames, new_chains):
        new_names = set(lig_resnames)
        for s in self._saved:
            if s["chains"] != new_chains and _rmsd_on_ca(s["tmp"], new_tmp) >= self.rmsd_thr:
                continue
            if new_names != set(s["ligand_res"]):
                continue
            if _rmsd_on_ligand(s["tmp"], new_tmp) < self.lig_rmsd_thr:
                return True
        return False

    @staticmethod
    def _merge_residues(residues, new_name):
        main = residues[0]
        main.resname = "LIG"
        taken = {a.get_name() for a in main.get_atoms()}
        
        for r in residues[1:]:
            if r is main:
                continue
            
            for at in list(r.get_atoms()):
                c = copy.copy(at)
                nn = _unique_name(at.name.strip(), taken)
                c.id = c.name = nn
                c.fullname = f"{nn:>4}"
                taken.add(nn)
                main.add(c)
            
            parent = r.get_parent()
            if parent is not None:
                parent.detach_child(r.id)
        return main

    def save(self, structure, pocket, out_dir, extra_lines, pdb_basename):
        lig_raw = [r for r in pocket.ligands if r.get_resname() != "HOH"]
        if not lig_raw:
            return False

        parent = list(range(len(lig_raw)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            pi, pj = find(i), find(j)
            if pi != pj:
                parent[pj] = pi

        for i in range(len(lig_raw)):
            for j in range(i + 1, len(lig_raw)):
                if _ligands_connected(
                    lig_raw[i], lig_raw[j], self.default_bond_distance
                ):
                    union(i, j)

        comps = {}
        for idx, r in enumerate(lig_raw):
            comps.setdefault(find(idx), []).append(r)

        merged_ligs = []
        group_names = []
        for comp in comps.values():
            counts = Counter(r.get_resname() for r in comp)
            unique_names = sorted(counts)
            if len(comp) > 1:
                new_name = "".join(unique_names)[:3]
                main = self._merge_residues(comp, new_name)
                merged_ligs.append(main)
            else:
                merged_ligs.append(comp[0])

            group_names.append(
                "-".join(
                    f"{res}{counts[res] if counts[res] > 1 else ''}"
                    for res in unique_names
                )
            )

        pocket.ligands = merged_ligs

        file_lig_code = (
            "_".join(sorted(group_names)) if len(group_names) > 1 else group_names[0]
        )
        chains_code = "_".join(sorted(pocket.chains))
        parts = [pdb_basename, file_lig_code, "chains", chains_code]
        if pocket.bond_type:
            parts.append(pocket.bond_type)
        root = "_".join(parts)

        out_path = out_dir / f"{root}.pdb"
        v = 2
        while out_path.exists():
            out_path = out_dir / f"{root}_v{v}.pdb"
            v += 1

        selector = _MergedSelect(pocket.ligands, pocket.chains)
        tmp = LigandPocketExtractor._save_temp(structure, selector)
        if self._is_duplicate(tmp, group_names, pocket.chains):
            Path(tmp).unlink(missing_ok=True)
            self._skipped += 1
            return False

        io = PDBIO()
        io.set_structure(structure)
        with open(out_path, "w") as fh:
            fh.writelines(extra_lines)
            io.save(fh, select=selector)
            fh.write("\n")

        self._saved.append(
            {
                "output": out_path,
                "tmp": tmp,
                "ligand_res": group_names,
                "chains": pocket.chains,
            }
        )
        logger.debug(f"Saved: {out_path}")
        return True

    @property
    def skipped(self):
        return self._skipped

class ProteinAnalyzer:
    def __init__(
        self,
        interaction_distance=4.5,
        ligand_cluster_distance=0.5,
        rmsd_threshold=2.0,
        default_bond_distance=1.5,
        short_peptide_length=10,
        ligand_rmsd_threshold=0.5,
        overlap_distance=0.5,
        limit_low=3,
        limit_high=5,
        protect_distance=3.0,
    ):
        self.extractor = LigandPocketExtractor(
            interaction_distance,
            ligand_cluster_distance,
            default_bond_distance,
            short_peptide_length,
        )
        self.writer = PocketWriter(
            rmsd_threshold,
            ligand_rmsd_threshold,
            ligand_cluster_distance,
            default_bond_distance,
        )
        self.overlap_distance = overlap_distance
        self.limit_low = limit_low
        self.limit_high = limit_high
        self.protect_distance = protect_distance

    @staticmethod
    def _read_extra_lines(pdb_filepath: Path):
        other = []
        with open(pdb_filepath) as fh:
            for line in fh:
                if not line.startswith(
                    ("ATOM", "HETATM", "END", "MASTER", "TER", "CONECT", "ANISOU")
                ):
                    other.append(line)
        return other

    @staticmethod
    def _parse_modres(extra_lines) -> FrozenSet[tuple[str, int, str]]:
        pos = set()
        for l in extra_lines:
            if l.startswith("MODRES"):
                chain = l[16].strip()
                seq   = l[18:22].strip()
                icode = l[22].strip() or " "
                if seq.isdigit():
                    pos.add((chain, int(seq), icode))
        return frozenset(pos)

    @staticmethod
    def _min_dist(r1, r2):
        a1 = np.array([a.coord for a in r1.get_atoms()])
        a2 = np.array([a.coord for a in r2.get_atoms()])
        return np.min(np.linalg.norm(a1[:, None, :] - a2[None, :, :], axis=-1))

    @staticmethod
    def _pocket_min_dist(p1, p2):
        a1 = np.array([a.coord for r in p1.ligands for a in r.get_atoms()])
        a2 = np.array([a.coord for r in p2.ligands for a in r.get_atoms()])
        return np.min(np.linalg.norm(a1[:, None, :] - a2[None, :, :], axis=-1))

    def _split_overlapping(self, pocket, prot_atoms):
        ligs = [r for r in pocket.ligands if r.get_resname() != "HOH"]
        if len(ligs) < 2:
            return [pocket]

        overlapped = {
            (i, j)
            for i in range(len(ligs))
            for j in range(i + 1, len(ligs))
            if self._min_dist(ligs[i], ligs[j]) < self.overlap_distance
        }
        if not overlapped:
            return [pocket]

        others = [
            r
            for r in ligs
            if all(
                (ligs.index(r), k) not in overlapped
                and (k, ligs.index(r)) not in overlapped
                for k in range(len(ligs))
            )
        ]

        pockets = []
        for lig in {ligs[i] for i, _ in overlapped}.union(
            {ligs[j] for _, j in overlapped}
        ):
            lig_atoms = [a for a in lig.get_atoms()]
            btype = self.extractor._get_pocket_bond_type(lig_atoms, prot_atoms)
            pockets.append(Pocket([lig] + others, pocket.chains, btype))
        return pockets

    def _filter_overabundant(self, pockets):
        counts = {}
        for p in pockets:
            for n in {r.get_resname() for r in p.ligands if r.get_resname() != "HOH"}:
                counts[n] = counts.get(n, 0) + 1
        if not counts:
            return pockets
        min_count = min(counts.values())
        limit = self.limit_low if min_count <= 2 else self.limit_high
        abundant = {n for n, c in counts.items() if c > limit}
        if not abundant:
            return pockets
        rare_pockets = [
            p
            for p in pockets
            if not {
                r.get_resname() for r in p.ligands if r.get_resname() != "HOH"
            }.issubset(abundant)
        ]
        kept, removed = [], 0
        for p in pockets:
            names = {r.get_resname() for r in p.ligands if r.get_resname() != "HOH"}
            if names and names.issubset(abundant):
                close_to_rare = any(
                    self._pocket_min_dist(p, rp) < self.protect_distance
                    for rp in rare_pockets
                )
                if not close_to_rare:
                    removed += 1
                    continue
            kept.append(p)
        logger.debug(
            f"Overabundant filter (min {min_count}, limit {limit}, protect {self.protect_distance} Å): removed {removed} pockets"
        )
        return kept

    def analyze(self, pdb_filepath, output_directory="separated_complexes"):
        pdb_filepath = Path(pdb_filepath)
        out_dir = Path(output_directory)
        out_dir.mkdir(exist_ok=True)
        logger.debug(f"Start: {pdb_filepath}")

        extra = self._read_extra_lines(pdb_filepath)
        mod_positions = self._parse_modres(extra)

        structure = PDBParser(QUIET=True).get_structure("pdb", pdb_filepath)

        raw = self.extractor.extract(structure, mod_positions)

        prot_atoms = [a for a in structure.get_atoms() if a.get_parent().id[0] == " "]

        pockets = []
        for p in raw:
            splitted = self._split_overlapping(p, prot_atoms)
            pockets.extend(splitted)
        pockets = self._filter_overabundant(pockets)
        logger.debug(f"Total pockets to export: {len(pockets)}")

        saved = 0
        try:
            for p in pockets:
                if self.writer.save(structure, p, out_dir, extra, pdb_filepath.stem):
                    saved += 1
        finally:
            for item in self.writer._saved:
                Path(item["tmp"]).unlink(missing_ok=True)
            self.writer._saved.clear()

        logger.debug("Done")
        return {"saved": saved, "skipped": self.writer.skipped}

def analyze_protein(
    pdb_filepath,
    output_directory="separated_complexes",
    interaction_distance=4.5,
    ligand_cluster_distance=4.5,
    rmsd_threshold=2.0,
    default_bond_distance=1.5,
    short_peptide_length=8,
    ligand_rmsd_threshold=1.0,
    overlap_distance=0.6,
):
    analyzer = ProteinAnalyzer(
        interaction_distance,
        ligand_cluster_distance,
        rmsd_threshold,
        default_bond_distance,
        short_peptide_length,
        ligand_rmsd_threshold,
        overlap_distance,
    )
    return analyzer.analyze(pdb_filepath, output_directory)

from pathlib import Path

def get_pdb_files(input_dir: Path):
    return list(input_dir.glob("*.pdb"))

from contextlib import contextmanager
from pathlib import Path

from joblib import Parallel, delayed, parallel
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

console = Console(force_terminal=True, width=100)
@contextmanager
def rich_joblib_progress(total: int):
    old_cb = parallel.BatchCompletionCallBack
    with Progress(
        SpinnerColumn(),
        TextColumn("Processing files:"),
        BarColumn(bar_width=None),
        TaskProgressColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=False,
        refresh_per_second=4,
    ) as progress:
        task_id = progress.add_task("processing", total=total)

        class RichCallBack(old_cb):
            def __call__(self, *args, **kwargs):
                res = old_cb.__call__(self, *args, **kwargs)
                progress.update(task_id, advance=self.batch_size)
                return res

        parallel.BatchCompletionCallBack = RichCallBack
        try:
            yield
        finally:
            parallel.BatchCompletionCallBack = old_cb


from joblib import parallel_backend
from pathlib import Path

def _run_single(pdb_file: Path, output_dir: Path, params: dict):
    result = {"saved": 0, "skipped": 0}
    try:
        result = analyze_protein(
            pdb_file,
            output_dir,
            interaction_distance=params["interaction_distance"],
            ligand_cluster_distance=params["ligand_cluster_distance"],
            rmsd_threshold=params["rmsd_threshold"],
            default_bond_distance=params["default_bond_distance"],
            short_peptide_length=params["short_peptide_length"],
            ligand_rmsd_threshold=params["ligand_rmsd_threshold"],
            overlap_distance=params["overlap_distance"],
        )
    except Exception as e:
        logger.exception(f"Failed to process {pdb_file}: {e}")
    finally:
        import gc
        gc.collect()
        with _PymolSession.locked():
            cmd.reinitialize()
    return result

def protein_ligand_separator(cfg):
    logger.info("\n========== Protein Ligand Separator ==========")

    input_dir = Path(cfg.paths.bioml_dir)
    output_dir = Path(cfg.paths.separated_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "interaction_distance": cfg.separator_params.interact_distance,
        "ligand_cluster_distance": cfg.separator_params.ligand_ligand_distance,
        "rmsd_threshold": cfg.separator_params.rmsd_threshold,
        "default_bond_distance": getattr(cfg.separator_params, "default_bond_distance", 1.5),
        "short_peptide_length": getattr(cfg.separator_params, "short_peptide_length", 8),
        "ligand_rmsd_threshold": getattr(cfg.separator_params, "ligand_rmsd_threshold", 1.0),
        "overlap_distance": getattr(cfg.separator_params, "overlap_distance", 0.6),
    }

    pdb_files = get_pdb_files(input_dir)
    total = len(pdb_files)
    logger.info(f"Total PDB files found: {total}")
    
    with parallel_backend("loky", n_jobs=cfg.n_jobs): 
        with rich_joblib_progress(total):
            results = Parallel(
                n_jobs=cfg.n_jobs,
                batch_size=1,
                pre_dispatch='n_jobs'
            )(delayed(_run_single)(p, output_dir, params) for p in pdb_files)
            
    saved = sum(r["saved"] for r in results)
    skipped = sum(r["skipped"] for r in results)
    logger.info(f"Total similar structures skipped: {skipped}")
    logger.info(f"Total structures SAVED: {saved}")