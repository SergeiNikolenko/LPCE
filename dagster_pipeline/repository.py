from dagster import repository
from pipelines.extraction_pipeline import pdb_extraction_job
from pipelines.junk_removal_pipeline import junk_removal_job

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@repository
def my_pipeline_repository():
    return [pdb_extraction_job, junk_removal_job]
