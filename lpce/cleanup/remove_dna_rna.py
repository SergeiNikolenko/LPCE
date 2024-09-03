import os
from pathlib import Path
import joblib
from tqdm import tqdm
from config.settings import PROCESSED_DIR

def remove_dna_rna_from_directory():
    """
    Removes files containing DNA or RNA sequences (A, T, G, C, U) in SEQRES or ATOM lines from the specified directory.
    """
    directory = PROCESSED_DIR
    files = list(Path(directory).glob("*.pdb"))
    total_files = len(files)
    retained_files = 0
    
    def contains_dna_rna_sequence(content):
        nucleotides = {"A", "T", "G", "C", "U"}
        seqres_lines = [line for line in content.splitlines() if line.startswith("SEQRES")]
        atom_lines = [line for line in content.splitlines() if line.startswith("ATOM")]

        # Check SEQRES lines
        for line in seqres_lines:
            sequence = line[19:].split()
            if any(n in nucleotides for n in sequence):
                return True

        # Check ATOM lines for specific nucleotide atoms
        for line in atom_lines:
            if line[17:20].strip() in nucleotides:
                return True

        return False

    def process_file(file):
        with open(file, "r") as f:
            content = f.read()
            if contains_dna_rna_sequence(content):
                os.remove(file)
                return False
        return True

    results = joblib.Parallel(n_jobs=-1)(joblib.delayed(process_file)(file) for file in tqdm(files))
    retained_files = sum(results)

    remaining_percentage = (retained_files / total_files) * 100

    print("=== DNA/RNA Removal Summary ===")
    print(f"Total files analyzed: {total_files:,}")
    print(f"Files retained after removal: {retained_files:,}")
    print(f"Percentage of files retained: {remaining_percentage:.2f}%")
    print("================================")

if __name__ == "__main__":
    remove_dna_rna_from_directory()
