from extraction.extract_complexes import extract_complexes
from extraction.decompress_files import decompress_pdb_files
from extraction.convert_pdb_to_smiles_sdf import convert_pdb_to_smiles_sdf
from extraction.parse_dict import extract_and_save_complexes_with_ligands

from cleanup.remove_water import remove_water_from_directory
from cleanup.remove_junk_ligands import remove_junk_ligands_from_directory
from cleanup.remove_dna_rna import remove_dna_rna_from_directory
from cleanup.filter_ligands import filter_ligands

from notifications.send_email import send_email_notification


def main():
    new_structures = extract_complexes()
    #decompress_pdb_files()
    remove_dna_rna_from_directory()
    #remove_water_from_directory()
    #remove_junk_ligands_from_directory()
    #convert_pdb_to_smiles_sdf()
    #extract_and_save_complexes_with_ligands()
    #filter_ligands()
    #send_email_notification(new_structures)


if __name__ == "__main__":
    main()
