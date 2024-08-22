from dagster import job
from pipelines.ops.extraction_ops import run_extract_complexes, send_email_notification
from pipelines.ops.decompression_ops import run_decompress_pdb_files
from pipelines.ops.water_removal_ops import run_remove_water
from pipelines.ops.cleanup_ops import run_remove_junk_ligands 

@job
def pdb_extraction_job():
    new_structures = run_extract_complexes()
    decompressed_structures = run_decompress_pdb_files(new_structures)
    no_water = run_remove_water(decompressed_structures)
    run_remove_junk_ligands(no_water)
    send_email_notification(new_structures)