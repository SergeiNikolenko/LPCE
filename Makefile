.PHONY: all run_pipeline test tests


all: clean run_pipeline

run_pipeline:
	python lpce/run_full_pipeline.py

test:
	clear
	python lpce/tests/test_pipeline.py

tests:
	clear
	export JUPYTER_PLATFORM_DIRS=1 && pytest lpce/tests/

clean:
	clear
	rm -rf .ruff_cache
	rm -rf /mnt/ligandpro/db/LPCE/processed
	rm -rf /mnt/ligandpro/db/LPCE/ligands
	rm -rf /mnt/ligandpro/db/LPCE/bioml
	rm -rf data/identical_groups.pkl
	rm -rf logs
	rm -rf data/filtered_ligands.json
	rm -rf data/site_info.json
	rm -rf data/grouped_complexes.json
	rm -rf data/removed_ligands_summary.json
