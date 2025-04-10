import os
from pathlib import Path

import joblib
from loguru import logger
from tqdm import tqdm


def count_models_in_file(file_path: Path) -> int:
    models = 0
    try:
        with open(file_path) as f:
            for line in f:
                if line.startswith("MODEL"):
                    models += 1
        return models
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return 0


def keep_only_first_model(file_path: Path):
    tmp_path = file_path.with_suffix(".tmp")
    try:
        with open(file_path) as inp, open(tmp_path, "w") as out:
            in_first = False
            for line in inp:
                if line.startswith("MODEL"):
                    if not in_first:
                        in_first = True
                        out.write(line)
                    else:
                        break
                else:
                    out.write(line)
                if line.startswith("ENDMDL") and in_first:
                    break
        os.replace(tmp_path, file_path)
    except Exception as e:
        logger.error(f"Error rewriting file {file_path}: {e}")
        if tmp_path.exists():
            tmp_path.unlink()


def process_file(file: Path) -> tuple[bool, bool, str]:
    models = count_models_in_file(file)
    multi_model = models > 1
    file_id = file.stem[3:] if file.stem[:3] == "pdb" else file.stem
    if multi_model:
        keep_only_first_model(file)
        modified = True
    else:
        modified = False
    return multi_model, modified, file_id


def remove_multiple_models_from_directory(cfg: object) -> dict:
    input_dir = Path(cfg.paths.processed_dir)
    logger.info("\n========== Checking PDB files for multiple models ==========")
    files = list(input_dir.glob("*.pdb"))
    total_files = len(files)
    logger.info(f"Total PDB files to analyze: {total_files}")

    results = joblib.Parallel(n_jobs=cfg.n_jobs)(
        joblib.delayed(process_file)(file)
        for file in tqdm(files, desc="Processing models")
    )

    multi_model_count = sum(r[0] for r in results)
    modified_files = [r[2] for r in results if r[1]]
    logger.info(f"Found multiple models in: {multi_model_count}")
    logger.info(f"Modified files: {len(modified_files)}")
    logger.info(
        f"Percentage of files retained: {100 - (multi_model_count / total_files * 100):.2f}%"
    )

    logger.info(f"Total files analyzed: {total_files}")
    return {"modified_files": modified_files}
