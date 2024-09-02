from cleanup.remove_water import remove_water_from_directory
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from cleanup.filter_ligands import filter_ligands
from extraction.parse_dict import extract_and_save_complexes_with_ligands


def run():
    #remove_water_from_directory()
    #remove_junk_ligands_from_directory()
    #extract_and_save_complexes_with_ligands()
    filter_ligands()

if __name__ == "__main__":
    run()
