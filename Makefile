.PHONY: all run_pipeline test tests pre-commit clean


all: clean tests_success clean run_pipeline end

run_pipeline:
	clear
	python lpce/run_full_pipeline.py

test:
	clear
	python lpce/tests/test_pipeline.py

tests:
	clear
	export JUPYTER_PLATFORM_DIRS=1 && pytest lpce/tests/

tests_success:
	@echo "Running tests..."
	@if $(MAKE) tests; then \
		echo "Tests passed! Proceeding with the rest of the pipeline..."; \
	else \
		echo "Tests failed! Stopping."; \
		exit 1; \
	fi

clean:
	clear
	rm -rf .pytest_cache
	rm -rf .ruff_cache
	rm -rf /mnt/ligandpro/db/LPCE/processed
	rm -rf /mnt/ligandpro/db/LPCE/ligands
	rm -rf /mnt/ligandpro/db/LPCE/bioml
	rm -rf /mnt/ligandpro/db/LPCE/separated

	rm -rf logs
	mkdir -p logs

	rm -rf data/identical_groups.json
	rm -rf data/filtered_ligands.json
	rm -rf data/site_info.json
	rm -rf data/grouped_complexes.json
	rm -rf data/removed_files.json

end:
	rm -rf /mnt/ligandpro/db/LPCE/raw
	rm -rf /mnt/ligandpro/db/LPCE/processed
	rm -rf /mnt/ligandpro/db/LPCE/ligands
	rm -rf /mnt/ligandpro/db/LPCE/bioml


pre-commit:
	pre-commit run --all-files
