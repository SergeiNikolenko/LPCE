from pathlib import Path


def clean_names(processed_dir):
    processed_dir = Path(processed_dir)
    for file in processed_dir.glob("*.pdb"):
        name_without_ext = file.stem
        if name_without_ext.startswith("pdb"):
            new_name = name_without_ext[3:] + ".pdb"
            new_file_path = file.with_name(new_name)
            file.rename(new_file_path)
            # print(f"Переименован файл: {file.name} -> {new_name}")


def clean_multiple_paths(cfg):
    clean_names(cfg.paths.processed_dir)
    clean_names(cfg.paths.bioml_dir)
    clean_names(cfg.paths.separated_dir)
