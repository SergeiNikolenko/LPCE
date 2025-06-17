# convert_pdb_to_smiles_sdf.py
import subprocess
import tempfile
from pathlib import Path

import requests
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def download_ligand_sdf(pdb_id: str, chain_id: str, res_seq: str, out_path: Path) -> bool:
    url = (
        f"https://models.rcsb.org/v1/{pdb_id}/ligand"
        f"?auth_asym_id={chain_id}&auth_seq_id={res_seq}&encoding=sdf"
    )
    try:
        r = requests.get(url, timeout=10)
        if r.ok and r.content.strip():
            out_path.write_bytes(r.content)
            return True
    except Exception as exc:  # noqa: BLE001
        logger.warning(f"Download failed for {url}: {exc}")
    return False


def process_pdb_file(pdb_file_path: Path, smiles_output_dir: Path, sdf_output_dir: Path) -> int:
    try:
        lines = pdb_file_path.read_text().splitlines()
        hetatm_blocks: dict[str, list[str]] = {}
        res_info: dict[str, tuple[str, str]] = {}
        for line in lines:
            if line.startswith("HETATM"):
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = line[22:26].strip()
                key = f"{res_name}_{chain_id}_{res_seq}"
                hetatm_blocks.setdefault(key, []).append(line)
                res_info[key] = (chain_id, res_seq)
            elif line.startswith("END"):
                break
        if not hetatm_blocks:
            return 0

        smiles_file = smiles_output_dir / pdb_file_path.with_suffix(".smi").name
        sdf_file = sdf_output_dir / pdb_file_path.with_suffix(".sdf").name

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir_path = Path(tmp_dir)
            pdb_id = pdb_file_path.stem[:4].lower()

            for key, het_lines in hetatm_blocks.items():
                chain_id, res_seq = res_info[key]
                ligand_sdf = tmp_dir_path / f"{key}.sdf"

                ok = download_ligand_sdf(pdb_id, chain_id, res_seq, ligand_sdf)
                if not ok:
                    ligand_pdb = tmp_dir_path / f"{key}.pdb"
                    ligand_pdb.write_text("\n".join(het_lines) + "\nEND\n")
                    subprocess.run(
                        ["obabel", str(ligand_pdb), "-O", str(ligand_sdf), "-osdf"],
                        capture_output=True,
                        check=False,
                    )

                if ligand_sdf.exists() and ligand_sdf.stat().st_size:
                    subprocess.run(
                        ["obabel", str(ligand_sdf), "-O", str(ligand_sdf.with_suffix('.smi')), "-osmi"],
                        capture_output=True,
                        check=False,
                    )
                    with ligand_sdf.open("rb") as src, sdf_file.open("ab") as dst:
                        dst.write(src.read())
                    smi_tmp = ligand_sdf.with_suffix(".smi")
                    if smi_tmp.exists():
                        smiles = smi_tmp.read_text().strip()
                        with smiles_file.open("a") as smi_out:
                            smi_out.write(f"{smiles}\t{key}\n")

        return len(hetatm_blocks)
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Failed to process {pdb_file_path}: {exc}")
        return 0


def convert_pdb_to_smiles_sdf(cfg) -> None:
    smiles_output_dir = Path(cfg.paths.ligands_dir) / "smiles"
    sdf_output_dir = Path(cfg.paths.ligands_dir) / "sdf"
    smiles_output_dir.mkdir(parents=True, exist_ok=True)
    sdf_output_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = list(Path(cfg.paths.processed_dir).glob("*.pdb"))
    logger.info("\n========== Converting PDB to SMILES and SDF ==========")
    logger.info(f"Processing {len(pdb_files)} PDB files from {cfg.paths.processed_dir}")

    results = Parallel(n_jobs=cfg.n_jobs)(
        delayed(process_pdb_file)(pdb, smiles_output_dir, sdf_output_dir)
        for pdb in tqdm(pdb_files, desc="Processing PDB files")
    )
    logger.info(f"Total proteins processed: {len(pdb_files)}")
    logger.info(f"Total ligands found: {sum(results)}")
