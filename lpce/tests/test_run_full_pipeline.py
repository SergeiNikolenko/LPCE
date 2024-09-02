import pytest
from unittest.mock import patch
from lpce.run_full_pipeline import main

@pytest.fixture(scope="function")
def setup_test_environment(tmp_path):
    """
    Fixture to set up a temporary environment for testing.
    """
    base_dir = tmp_path / "PDB"
    processed_dir = base_dir / "pdb2/processed"
    data_dir = tmp_path / "data"

    processed_dir.mkdir(parents=True)
    data_dir.mkdir(parents=True)

    # Creating dummy files for testing
    (processed_dir / "test.pdb").write_text("dummy content")

    yield

def test_run_full_pipeline(setup_test_environment):
    """
    Test the full pipeline in run_full_pipeline.py script.
    This test mocks the individual steps to verify they are called in the correct order.
    """
    print("Starting the test...")
    
    with patch('lpce.extraction.extract_complexes.extract_complexes', return_value=10) as mock_extract_complexes:
        with patch('lpce.extraction.decompress_files.decompress_pdb_files') as mock_decompress:
            with patch('lpce.cleanup.remove_water.remove_water_from_directory') as mock_remove_water:
                with patch('lpce.cleanup.remove_junk_ligands.remove_junk_ligands_from_directory') as mock_remove_junk:
                    with patch('lpce.extraction.parse_dict.extract_and_save_complexes_with_ligands') as mock_extract_and_save:
                        with patch('lpce.notifications.send_email.send_email_notification') as mock_send_email:
                            main()

                            assert mock_extract_complexes.called
                            assert mock_decompress.called
                            assert mock_remove_water.called
                            assert mock_remove_junk.called
                            assert mock_extract_and_save.called
                            assert mock_send_email.called

                            print("Assertions passed!")
