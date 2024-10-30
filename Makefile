.PHONY: all run_pipeline test tests


all: run_pipeline foldseek

run_pipeline:
	clear
	rm -rf /mnt/ligandpro/db/LPCE/processed
	rm -rf /mnt/ligandpro/db/LPCE/ligands
	rm -rf /mnt/ligandpro/db/LPCE/bioml
	rm -rf data/identical_groups.pkl
	rm -rf logs
	python lpce/run_full_pipeline.py

test:
	clear
	python lpce/tests/test_pipeline.py

tests:
	clear
	export JUPYTER_PLATFORM_DIRS=1 && pytest lpce/tests/
