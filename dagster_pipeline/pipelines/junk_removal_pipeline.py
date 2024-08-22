from dagster import job, In, Nothing
from pipelines.ops.water_removal_ops import run_remove_water
from pipelines.ops.cleanup_ops import run_remove_junk_ligands

@job
def junk_removal_job():
    no_water = run_remove_water()
    run_remove_junk_ligands(no_water)
