import json
from pathlib import Path


def filter_ligands():
    """
    Filters ligands based on their presence in site
    information and saves updated complexes.
    """
    with open(Path("data/grouped_complexes.json")) as f:
        grouped_complexes = json.load(f)

    with open(Path("data/site_info.json")) as f:
        site_info = json.load(f)

    total_proteins_grouped = len(grouped_complexes)
    total_ligands_grouped = sum(
        len(ligands)
        for chains in grouped_complexes.values()
        for ligands in chains.values()
    )

    total_proteins_site_info = len(site_info)
    ligand_intersections = 0
    ligand_deletions = 0
    remaining_ligands = 0

    for protein, chains in grouped_complexes.items():
        if protein not in site_info:
            remaining_ligands += sum(len(ligands) for ligands in chains.values())
            continue  # Skip protein if it's not in site_info.json
        for chain, ligands in chains.items():
            ligands_to_keep = []
            for ligand_info in ligands:
                ligand_name = ligand_info["ligand"]
                if (
                    ligand_name in site_info[protein]
                    and chain in site_info[protein][ligand_name]
                ):
                    ligands_to_keep.append(ligand_info)
                    ligand_intersections += 1
                else:
                    ligand_deletions += 1
            # Update the list of ligands in the chain
            grouped_complexes[protein][chain] = ligands_to_keep
            remaining_ligands += len(ligands_to_keep)

    # Save updated data if needed
    with open(Path("data/site_removed_complexes.json"), "w") as f:
        json.dump(grouped_complexes, f, indent=4)

    percent_deleted = (ligand_deletions / total_ligands_grouped) * 100

    # Print statistics
    print("=== Protein and Ligand Filtering Summary ===")
    print(f"Total proteins analyzed: {total_proteins_grouped:,}")
    print(f"Total ligands analyzed: {total_ligands_grouped:,}")
    print()
    print(f"Proteins with site info available: {total_proteins_site_info:,}")
    print(f"Relevant ligands found in sites: {ligand_intersections:,}")
    print()
    print(f"Ligands removed during filtering: {ligand_deletions:,}")
    print(f"Percentage of ligands removed: {percent_deleted:.1f}%")
    print()
    print(f"Ligands remaining after filtering: {remaining_ligands:,}")
    print("===========================================")


if __name__ == "__main__":
    filter_ligands()
