.PHONY: all run_pipeline test tests pre-commit clean tmux

CONFIG_NAME ?= config


all: clean tests_success clean run_pipeline clean end

tmux:
	@tmux new-session -d -s lpce "make all CONFIG_FILE=$(CONFIG_NAME)"
	@tmux attach -t lpce

run_pipeline:
	clear
	python lpce/run_full_pipeline.py $(CONFIG_NAME)

test:
	clear
	python lpce/tests/test_pipeline.py $(CONFIG_NAME)

tests:
	clear
	export JUPYTER_PLATFORM_DIRS=1 && python -m pytest lpce/tests/ --config-name=$(CONFIG_NAME)

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
	rm -rf data/removed_files_tests.json

end:
	rm -rf /mnt/ligandpro/db/LPCE/raw
	rm -rf /mnt/ligandpro/db/LPCE/processed
	rm -rf /mnt/ligandpro/db/LPCE/ligands
	rm -rf /mnt/ligandpro/db/LPCE/bioml

pre-commit: venv
	$(VENV_BIN)/pre-commit run --all-files
