from cleanup.remove_water import remove_water_from_directory
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory

def run():
    remove_water_from_directory()
    remove_junk_ligands_from_directory()

if __name__ == "__main__":
    run()
