import sys
import warnings
import tempfile
import threading
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import List, Set

from Bio.PDB import PDBParser, NeighborSearch, Select, PDBIO
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm
import pymol
from pymol import cmd

# ───────────────────  basic setup  ───────────────────
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.core.universe")

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


# ───────────────────  PyMOL helpers  ───────────────────
class _PymolSession:
    _lock = threading.Lock()
    _started = False

    @classmethod
    def ensure(cls):
        with cls._lock:
            if not cls._started:
                pymol.finish_launching(["pymol", "-qc"])
                cls._started = True

    @classmethod
    def locked(cls):
        cls.ensure()
        return cls._lock


def _rmsd_on_ca(p1: str, p2: str) -> float:
    with _PymolSession.locked():
        o1, o2 = f"o{uuid.uuid4().hex}", f"o{uuid.uuid4().hex}"
        cmd.delete("all")
        cmd.load(p1, o1)
        cmd.load(p2, o2)
        try:
            return float(cmd.align(f"{o1} and name CA", f"{o2} and name CA")[0])
        finally:
            cmd.delete("all")


def _rmsd_on_ligand(p1: str, p2: str) -> float:
    with _PymolSession.locked():
        o1, o2 = f"o{uuid.uuid4().hex}", f"o{uuid.uuid4().hex}"
        cmd.delete("all")
        cmd.load(p1, o1)
        cmd.load(p2, o2)
        cmd.align(f"{o1} and name CA", f"{o2} and name CA", cycles=0)
        s1 = f"{o1} and hetatm and not resn HOH"
        s2 = f"{o2} and hetatm and not resn HOH"
        if cmd.count_atoms(s1) == 0 or cmd.count_atoms(s2) == 0:
            cmd.delete("all")
            return float("inf")
        rms = float(cmd.rms_cur(s1, s2))
        cmd.delete("all")
        return rms


# ───────────────────  selectors & DTO  ───────────────────
class _LigandPocketSelect(Select):
    def __init__(self, ligands, chains):
        self._ids = {(r.get_parent().id, r.id) for r in ligands}
        self._chains = chains

    def accept_chain(self, chain):
        return chain.id in self._chains or any(chain.id == cid for cid, _ in self._ids)

    def accept_residue(self, residue):
        cid = residue.get_parent().id
        return (cid, residue.id) in self._ids or residue.id[0] == " "


@dataclass
class Pocket:
    ligands: List
    chains: Set[str]
    covalent: bool


# ───────────────────  core extraction  ───────────────────
class LigandPocketExtractor:
    def __init__(self, interact_d=4.5, cluster_d=3.0, bond_d=2.0, short_pep=10):
        self.interact_d = interact_d
        self.cluster_d = cluster_d
        self.bond_d = bond_d
        self.short_pep = short_pep

    @staticmethod
    def _save_tmp(structure, sel):
        io = PDBIO()
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        io.set_structure(structure)
        io.save(tmp.name, select=sel)
        tmp.close()
        return tmp.name

    @staticmethod
    def _chains_near(lig_atoms, prot_atoms, d):
        ns = NeighborSearch(prot_atoms)
        chains = {
            a.get_parent().get_parent().id
            for lig in lig_atoms
            for a in ns.search(lig.coord, d, level="A")
        }
        return chains

    def extract(self, structure):
        prot_atoms = [a for a in structure.get_atoms() if a.get_parent().id[0] == " "]
        chain_res = {}
        for r in structure.get_residues():
            chain_res.setdefault(r.get_parent().id, []).append(r)

        candidates = [
            r for r in structure.get_residues()
            if r.id[0] != " " or len(chain_res[r.get_parent().id]) <= self.short_pep
        ]

        idx = {r: i for i, r in enumerate(candidates)}
        parent = list(range(len(candidates)))

        def f(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def u(i, j):
            pi, pj = f(i), f(j)
            if pi != pj:
                parent[pj] = pi

        ns = NeighborSearch([a for r in candidates for a in r.get_atoms()])
        for r in candidates:
            i = idx[r]
            for a in r.get_atoms():
                for n in ns.search(a.coord, self.cluster_d, level="A"):
                    r2 = n.get_parent()
                    if r2 is r or r2 not in idx:
                        continue
                    u(i, idx[r2])

        comps = {}
        for r in candidates:
            comps.setdefault(f(idx[r]), []).append(r)

        p_ns = NeighborSearch(prot_atoms)
        pockets = []
        for res in comps.values():
            lig_atoms = [a for rr in res for a in rr.get_atoms()]
            chains = self._chains_near(lig_atoms, prot_atoms, self.interact_d)
            if not chains:
                continue
            cov = any(p_ns.search(a.coord, self.bond_d, level="A") for a in lig_atoms)
            pockets.append(Pocket(res, chains, cov))
        return pockets


# ───────────────────  writer with statistics  ───────────────────
class PocketWriter:
    def __init__(self, rmsd_thr=2.0, lig_rmsd_thr=0.5):
        self.rmsd_thr = rmsd_thr
        self.lig_rmsd_thr = lig_rmsd_thr
        self.saved = 0
        self.skipped = 0
        self._seen = []

    def _is_dup(self, new_tmp, lig_res):
        for s in self._seen:
            if _rmsd_on_ca(s["tmp"], new_tmp) >= self.rmsd_thr:
                continue
            if s["lig_and_res"] != lig_res:
                continue
            if _rmsd_on_ligand(s["tmp"], new_tmp) < self.lig_rmsd_thr:
                return True
        return False

    def save(self, structure, pocket, out_dir, extra, base):
        lig_res = sorted({l.get_resname() for l in pocket.ligands})
        chains_s = "_".join(sorted(pocket.chains))
        lig_s = "_".join(lig_res) + ("_COV" if pocket.covalent else "")
        root = f"{base}_{lig_s}_chains_{chains_s}"
        out = out_dir / f"{root}.pdb"
        v = 2
        while out.exists():
            out = out_dir / f"{root}_v{v}.pdb"
            v += 1

        sel = _LigandPocketSelect(pocket.ligands, pocket.chains)
        tmp = LigandPocketExtractor._save_tmp(structure, sel)

        if self._is_dup(tmp, lig_res):
            Path(tmp).unlink(missing_ok=True)
            self.skipped += 1
            return

        io = PDBIO()
        io.set_structure(structure)
        with open(out, "w") as f:
            f.writelines(extra)
            io.save(f, select=sel)
            f.write("END\n")

        self._seen.append({"tmp": tmp, "lig_and_res": lig_res})
        self.saved += 1
        logger.debug(f"Saved -> {out}")


# ───────────────────  high‑level analyzer  ───────────────────
class ProteinAnalyzer:
    def __init__(
        self,
        interact_d=4.5,
        cluster_d=3.0,
        rmsd_thr=2.0,
        bond_d=2.0,
        short_pep=10,
        lig_rmsd_thr=0.5,
    ):
        self.extractor = LigandPocketExtractor(interact_d, cluster_d, bond_d, short_pep)
        self.writer = PocketWriter(rmsd_thr, lig_rmsd_thr)

    @staticmethod
    def _extra_lines(p):
        other = []
        with open(p) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM", "END", "MASTER", "TER", "CONECT")):
                    other.append(line)
        return other

    def analyze(self, pdb_path, out_dir="separated_complexes"):
        pdb_path = Path(pdb_path)
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        extra = self._extra_lines(pdb_path)
        structure = PDBParser(QUIET=True).get_structure("pdb", pdb_path)
        for p in self.extractor.extract(structure):
            self.writer.save(structure, p, out_dir, extra, pdb_path.stem)
        return {
            "structures_saved": self.writer.saved,
            "structures_skipped_similar": self.writer.skipped,
        }


# ───────────────────  functional wrapper  ───────────────────
def analyze_protein(
    pdb_path,
    out_dir="separated_complexes",
    interact_d=4.5,
    cluster_d=1.0,
    rmsd_thr=2.0,
    bond_d=1.0,
    short_pep=3,
    lig_rmsd_thr=0.1,
):
    return ProteinAnalyzer(
        interact_d,
        cluster_d,
        rmsd_thr,
        bond_d,
        short_pep,
        lig_rmsd_thr,
    ).analyze(pdb_path, out_dir)


# ───────────────────  batch‑level API  ───────────────────
def get_pdb_files(directory):
    return list(Path(directory).glob("*.pdb"))


def protein_ligand_separator(cfg):
    in_dir = Path(cfg.paths.bioml_dir)
    out_dir = Path(cfg.paths.separated_dir)

    p = cfg.separator_params
    pdb_files = get_pdb_files(in_dir)

    logger.info("\n========== Protein ligand separator started ==========")

    results = Parallel(n_jobs=cfg.n_jobs)(
        delayed(analyze_protein)(
            pdb,
            out_dir,
            p.interact_distance,
            p.ligand_ligand_distance,
            p.rmsd_threshold,
            p.bond_distance,
            p.short_pep,
            p.lig_rmsd_thr,
        )
        for pdb in tqdm(pdb_files, desc="Separating ligand pockets")
    )

    total_saved = sum(r["structures_saved"] for r in results)
    total_skipped = sum(r["structures_skipped_similar"] for r in results)

    logger.info(f"Total PDB files: {len(pdb_files)}")
    logger.info(f"Total similar skipped: {total_skipped}")
    logger.info(f"Total structures SAVED: {total_saved}")
