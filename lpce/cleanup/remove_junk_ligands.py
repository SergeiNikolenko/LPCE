import json
import multiprocessing
from collections import Counter
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from loguru import logger
from tqdm import tqdm


def parse_atom_line(line: str) -> tuple:
    """
    Извлекает координаты и идентификатор лиганда из строки с атомом в PDB файле.

    Args:
        line (str): Строка из PDB файла, описывающая атом.

    Returns:
        tuple: Кортеж с идентификатором лиганда и координатами атома (x, y, z).
    """
    ligand_id = line[17:20].strip()
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())
    return ligand_id, np.array([x, y, z])


def remove_junk_ligands_from_file(input_file_path: Path, junk_ligands: set) -> dict:
    """
    Удаляет мусорные лиганды из PDB файла, если в радиусе 3 Å нет других лигандов.

    Args:
        input_file_path (Path): Путь к PDB файлу.
        junk_ligands (set): Набор идентификаторов мусорных лигандов.

    Returns:
        dict: Словарь с количеством удаленных лигандов.
    """
    try:
        ligand_counts = Counter()
        atoms = []
        ligands_to_remove = set()
        temp_file_path = input_file_path.with_suffix(".tmp")
        changes_made = False
        threshold = 3.0

        # Читаем файл и собираем атомы
        with open(input_file_path) as f_in:
            lines = f_in.readlines()

        for line in lines:
            if line.startswith("HETATM"):
                ligand_id, coordinates = parse_atom_line(line)
                atoms.append((ligand_id, coordinates))

        # Ищем лиганды, которые можно удалить
        for i, (ligand_id, coord) in enumerate(atoms):
            if ligand_id in junk_ligands:
                nearby_ligand_found = False
                for j, (other_ligand_id, other_coord) in enumerate(atoms):
                    if i != j and np.linalg.norm(coord - other_coord) <= threshold:
                        if other_ligand_id not in junk_ligands:
                            nearby_ligand_found = True
                            break
                if not nearby_ligand_found:
                    ligands_to_remove.add(ligand_id)

        # Пишем обновленный файл без удаленных лигандов
        with open(temp_file_path, "w") as f_out:
            for line in lines:
                if line.startswith("HETATM"):
                    ligand_id, _ = parse_atom_line(line)
                    if ligand_id in ligands_to_remove:
                        ligand_counts[ligand_id] += 1
                        changes_made = True
                        continue
                f_out.write(line)

        if changes_made:
            temp_file_path.replace(input_file_path)
        else:
            temp_file_path.unlink()

        return dict(ligand_counts)
    except Exception as e:
        logger.error(f"Failed to process {input_file_path}: {e}")
        return {"error": f"Error processing {input_file_path}: {e}"}


def remove_junk_ligands_from_directory(cfg) -> dict:
    """
    Обрабатывает все PDB файлы в указанной директории и удаляет мусорные лиганды, если они не находятся рядом с другими лигандами.

    Args:
        cfg (DictConfig): Конфигурация с путями к файлам и настройками логирования.

    Returns:
        dict: Результаты обработки, включая количество удаленных лигандов и ошибок.
    """
    input_directory = Path(cfg.paths.processed_dir)
    junk_ligands_file = Path(cfg.output_files.trash_ligands_json)

    logger.info("========== Removing Junk Ligands ==========")

    # Загружаем список мусорных лигандов
    with open(junk_ligands_file) as file:
        junk_ligands = set(json.load(file))

    pdb_files = list(input_directory.rglob("*.pdb"))
    total_files = len(pdb_files)

    logger.info(f"Found {total_files} PDB files in {input_directory}")

    num_cores = multiprocessing.cpu_count() - 1

    # Параллельная обработка файлов
    results = Parallel(n_jobs=num_cores)(
        delayed(remove_junk_ligands_from_file)(pdb_file, junk_ligands)
        for pdb_file in tqdm(
            pdb_files, desc="Removing junk ligands", unit="file", total=total_files
        )
    )

    # Суммируем результаты
    total_ligand_counts = Counter()
    failed_files = []

    for i, result in enumerate(results):
        if isinstance(result, dict) and "error" not in result:
            total_ligand_counts.update(result)
        else:
            failed_files.append(str(pdb_files[i]))

    successful_files = total_files - len(failed_files)

    logger.info(f"Total structures processed: {total_files}")
    logger.info(f"Successfully processed: {successful_files}")
    logger.info(f"Failed to process: {len(failed_files)}")
    logger.info(f"Total ligands removed: {sum(total_ligand_counts.values())}")

    return {
        "total_files": total_files,
        "successful_files": successful_files,
        "failed_files": failed_files,
        "total_ligands_removed": sum(total_ligand_counts.values()),
        "ligands_removed_per_type": dict(total_ligand_counts),
    }
