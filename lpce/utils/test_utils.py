from pathlib import Path
import shutil
import sys
from typing import List, Dict, Any

from loguru import logger

def setup_logging(cfg: Dict[str, Any]) -> None:
    log_file_path = Path(cfg.logging.log_file_tests)
    if log_file_path.exists():
        log_file_path.unlink()

    logger.remove()
    logger.add(sys.stdout, format="{message}", level="INFO")
    logger.add(
        cfg.logging.log_file_tests,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="INFO",
    )

def cleanup_directories(directories: List[Path]) -> None:
    for directory in directories:
        if directory.exists():
            shutil.rmtree(directory)

def setup_test_directories(temp_path: Path) -> Dict[str, Path]:
    directories = {
        "processed": temp_path / "processed",
        "ligands": temp_path / "ligands",
        "logs": temp_path / "logs",
        "final": temp_path / "final",
        "bioml": temp_path / "bioml",
        "separated": temp_path / "separated"
    }
    
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    
    return directories

def update_test_config(cfg: Dict[str, Any], directories: Dict[str, Path]) -> Dict[str, Any]:
    test_cfg = cfg.copy()
    test_cfg.paths.processed_dir = str(directories["processed"])
    test_cfg.paths.bioml_dir = str(directories["bioml"])
    test_cfg.paths.ligands_dir = str(directories["ligands"])
    test_cfg.paths.separated_dir = str(directories["separated"])
    test_cfg.paths.final_dir = str(directories["final"])
    test_cfg.logging.log_file = str(directories["logs"] / "foldseek.log")
    return test_cfg

def copy_results_to_final(
    processed_dir: Path,
    test_cfg: Dict[str, Any],
    final_dirs: Dict[str, Path]
) -> None:
    shutil.copytree(processed_dir, final_dirs["processed"], dirs_exist_ok=True)
    shutil.copytree(test_cfg.paths.bioml_dir, final_dirs["bioml"], dirs_exist_ok=True)
    shutil.copytree(test_cfg.paths.separated_dir, final_dirs["separated"], dirs_exist_ok=True)
    shutil.copytree(test_cfg.paths.final_dir, final_dirs["final"], dirs_exist_ok=True) 