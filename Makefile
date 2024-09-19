.PHONY: all run_pipeline extract_complexes remove_dna_rna remove_water remove_junk_ligands convert_pdb filter_ligands

# Запуск полного пайплайна
all: run_pipeline

# Запуск полного пайплайна
run_pipeline:
	python lpce/run_full_pipeline.py

# Извлечение комплексов
extract_complexes:
	python lpce/extraction/extract_complexes.py

# Удаление ДНК/РНК
remove_dna_rna:
	python lpce/cleanup/remove_dna_rna.py

# Удаление молекул воды
remove_water:
	gcc -o lpce/cleanup/remove_water lpce/cleanup/remove_water.c
	python lpce/cleanup/remove_water.py

# Удаление мусорных лигандов
remove_junk_ligands:
	python lpce/cleanup/remove_junk_ligands.py

# Конвертация PDB в SMILES и SDF
convert_pdb:
	python lpce/extraction/convert_pdb_to_smiles_sdf.py

# Фильтрация лигандов
filter_ligands:
	python lpce/cleanup/filter_ligands.py
