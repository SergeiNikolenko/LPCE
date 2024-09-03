from pathlib import Path
import subprocess
from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile
from config.settings import PROCESSED_DIR, OUTPUT_LIG

def process_pdb_file(pdb_file_path, smiles_output_dir, sdf_output_dir):
    """
    Processes a single PDB file to extract HETATM blocks, converts them to SMILES and SDF formats, 
    and saves the results to the specified output directories.

    Args:
        pdb_file_path (Path): Path to the input PDB file.
        smiles_output_dir (Path): Directory where SMILES files will be saved.
        sdf_output_dir (Path): Directory where SDF files will be saved.

    Returns:
        int: The number of ligands found in the PDB file.
    """
    with pdb_file_path.open('r') as file:
        lines = file.readlines()
    
    hetatm_blocks = {}
    ligands_found = 0
    
    for line in lines:
        if line.startswith("HETATM"):
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            key = f"{res_name}_{chain_id}"
            if key not in hetatm_blocks:
                hetatm_blocks[key] = []
            hetatm_blocks[key].append(line)
        elif line.startswith("END"):
            break
    
    if not hetatm_blocks:
        return 0
    
    ligands_found = len(hetatm_blocks)
    
    with tempfile.TemporaryDirectory() as hetatm_dir:
        hetatm_dir = Path(hetatm_dir)
        for key, hetatm_lines in hetatm_blocks.items():
            output_file = hetatm_dir / f"{key}.pdb"
            with output_file.open('w') as out_file:
                out_file.writelines(hetatm_lines)
        
        smiles_file = smiles_output_dir / pdb_file_path.with_suffix('.smi').name
        sdf_file = sdf_output_dir / pdb_file_path.with_suffix('.sdf').name
        
        with smiles_file.open('w') as out_smiles:
            for pdb_file in hetatm_dir.glob("*.pdb"):
                temp_smiles_file = pdb_file.with_suffix('.smi')
                temp_sdf_file = pdb_file.with_suffix('.sdf')
                
                command_smiles = ["obabel", str(pdb_file), "-O", str(temp_smiles_file), "-osmi"]
                subprocess.run(command_smiles, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                command_sdf = ["obabel", str(pdb_file), "-O", str(temp_sdf_file), "-osdf"]
                subprocess.run(command_sdf, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                if temp_smiles_file.exists():
                    with temp_smiles_file.open('r') as temp_file:
                        for line in temp_file:
                            smiles_line = f"{line.strip()}\t{pdb_file.stem}\n"
                            out_smiles.write(smiles_line)
                
                if temp_sdf_file.exists():
                    with sdf_file.open('a') as out_sdf:
                        with temp_sdf_file.open('r') as temp_file:
                            out_sdf.write(temp_file.read())
    
    return ligands_found

def convert_pdb_to_smiles_sdf():
    """
    Processes all PDB files in the input directory (defined in settings) by extracting ligands and converting
    them to SMILES and SDF formats. The results are saved in the output directories defined in settings.

    This function does not require any arguments as it uses paths defined in the configuration file.

    It also prints the total number of processed proteins and the number of ligands found.
    """
    input_dir = Path(PROCESSED_DIR)
    output_dir = Path(OUTPUT_LIG)

    smiles_output_dir = output_dir / "smi_files"
    sdf_output_dir = output_dir / "sdf_files"
    
    smiles_output_dir.mkdir(parents=True, exist_ok=True)
    sdf_output_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_files = list(input_dir.glob("*.pdb"))
    
    results = Parallel(n_jobs=-1)(
        delayed(process_pdb_file)(pdb_file, smiles_output_dir, sdf_output_dir)
        for pdb_file in tqdm(pdb_files, desc="Processing PDB files")
    )
    
    total_proteins = len(pdb_files)
    total_ligands = sum(results)
    
    print(f"Total proteins processed: {total_proteins}")
    print(f"Total ligands found: {total_ligands}")

if __name__ == "__main__":
    convert_pdb_to_smiles_sdf()
