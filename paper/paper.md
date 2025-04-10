
**Title:**

**LPCE: An Automated Pipeline for Ligand Processing and Cleaning in Protein Structures**

**Author:**

Sergei A. Nikolenko

**Abstract**

The accurate analysis of protein-ligand interactions is essential in computational biology, drug discovery, and structural bioinformatics. However, existing datasets such as PDBbind often contain inconsistencies, outdated entries, and unwanted components that can hinder computational analyses and model development. We present the Ligand Processing and Cleaning Engine (LPCE), an automated pipeline designed to extract, clean, and process ligands from Protein Data Bank (PDB) structures. LPCE systematically removes unwanted components such as water molecules and non-relevant ligands, filters out redundant and similar structures, and provides a curated dataset suitable for various computational applications. This paper describes the methodology behind LPCE, evaluates its performance, and discusses how it improves upon existing datasets.

**Introduction**

Protein-ligand interactions play a critical role in many biological processes and are a central focus in drug discovery efforts. Structural databases like the Protein Data Bank (PDB) [1] provide a wealth of information on protein-ligand complexes. However, these databases often include redundant entries, inconsistent data annotations, and unwanted components such as water molecules, ions, and non-relevant ligands. Datasets like PDBbind [2], which aggregate protein-ligand complexes from the PDB, inherit these issues and may also be outdated due to the rapid expansion of available structural data.

These inconsistencies can significantly impact computational analyses, including molecular docking, virtual screening, and machine learning model development. Cleaning and standardizing these datasets is essential to ensure reliable and reproducible results. While some efforts have been made to curate these datasets, there remains a need for automated tools that can systematically process and clean protein-ligand structural data.

In this paper, we introduce the Ligand Processing and Cleaning Engine (LPCE), an automated pipeline designed to address these challenges. LPCE provides a systematic approach to extract relevant ligands from PDB structures, remove unwanted components, and filter out redundant or similar entries. By generating a high-quality, up-to-date dataset, LPCE facilitates more accurate computational analyses and model development.
рассказииь про то какеи сейчас есть бд для докинга
почему они проблемыные и тд
почему решили свою делать

вот как подошли к этому

**Methods**

LPCE is composed of several sequential steps, each designed to clean and refine the dataset of protein-ligand complexes. The pipeline is implemented using Python scripts and relies on open-source tools such as Open Babel [3] and Foldseek [4] for specific tasks. Below, we describe each step in detail.

1. **Data Acquisition**

   Protein structures are downloaded from the PDB repository. Structures are stored locally in compressed `.ent.gz` format.

2. **Decompression**

   The compressed PDB files are decompressed to obtain `.pdb` files suitable for processing.

3. **Removal of DNA and RNA Structures**

   Structures containing nucleic acids are identified and removed to focus the dataset on protein-ligand interactions. The pipeline scans the `SEQRES` and `ATOM` records in each PDB file for nucleotide residues (adenine, thymine, guanine, cytosine, uracil) and their derivatives. Files containing these residues are discarded.

4. **Elimination of Multiple Models**

   PDB files containing multiple models (e.g., NMR structures with multiple conformations) are removed. The pipeline counts occurrences of the `MODEL` keyword in each file and discards files with more than one model.

5. **Water Molecule Removal**

   Water molecules are removed to prevent interference with analyses that do not require solvation effects. This step is performed efficiently using a custom C program integrated into the pipeline via a Python wrapper.

6. **Removal of Unwanted Ligands**

   Non-relevant ligands, such as buffer components, ions, and other small molecules not of interest, are removed based on a predefined list (`trash_ligands.json`). This list was curated through exploratory data analysis and manual inspection, including analysis in the Jupyter notebook `ligand_eda.ipynb`.

7. **Conversion to SMILES and SDF Formats**

   Remaining ligand structures are converted into SMILES and SDF formats using Open Babel. This facilitates further computational analyses, descriptor calculations, and compatibility with cheminformatics tools.

8. **Extraction of Protein-Ligand Complexes**

   The pipeline extracts complexes containing both proteins and their associated ligands. Information about binding sites is parsed from PDB files and stored in a JSON format (`site_info.json`) for reference.

9. **Binding Site Filtering**

   Complexes are filtered to retain only those ligands located within known binding sites. If binding site information is unavailable, the complex is retained for further analysis.

10. **Cleanup of Unused Files**

    PDB files that do not contain relevant ligands or have been processed and are no longer needed are deleted to conserve storage space.

11. **Biomolecule Separation**

    Structures are separated into unique biomolecules based on `REMARK 350` information in PDB files. The pipeline extracts information about chains for each biomolecule and removes duplicates or structures where one sequence includes another. Unique structures are retained for further processing.

12. **Protein-Ligand Complex Separation**

    The pipeline identifies and extracts individual protein-ligand complexes:

    - **Identification of Interacting Chains and Ligands**: For each ligand, the pipeline determines which protein chains are within a specified distance threshold (e.g., 4.5 Å). Only chains interacting with the ligand are retained in the complex.

    - **Grouping of Ligands**: Ligands that are close to each other (e.g., within 4 Å) are grouped together to preserve functional ligand groups.

    - **Filtering of Similar Structures**: Complexes are compared based on sequence identity and root-mean-square deviation (RMSD). Structures with high similarity (e.g., sequence identity >98% and RMSD <2 Å) are considered redundant, and only one representative is retained.

    - **Saving of Unique Complexes**: Each unique protein-ligand complex is saved as a separate PDB file, containing only the interacting chains and ligands.

13. **Path Cleanup**

    File names are standardized to a consistent format. For example, files named `pdbXXXX.pdb` are renamed to `XXXX.pdb`, where `XXXX` is the PDB identifier.

14. **Duplicate Removal with Foldseek**

    Foldseek is used to identify and remove redundant structures:

    - **Database Creation**: A Foldseek database is created from the PDB files.

    - **Structure Comparison**: Structures are compared using TM-score and sequence identity thresholds.

    - **Group Formation**: Identical or highly similar structures are grouped together.

    - **Result Storage**: Groups of duplicates are saved in JSON format for reference.

15. **Removal of Similar Structures**

    Within each group identified by Foldseek:

    - Structures are sorted by resolution, with higher-resolution structures preferred.

    - Only the structure with the best resolution is retained.

    - Remaining structures in the group are deleted.

16. **Removal of Surface Ligands**

    Ligands that are not sufficiently buried within the protein structure are removed. For each ligand, the pipeline calculates the minimum distance from each ligand atom to protein atoms. If less than a specified percentage (e.g., 30%) of ligand atoms are within a threshold distance (e.g., 5 Å) of protein atoms, the ligand is considered surface-bound and is removed.

17. **Addition of Hydrogen Atoms**

    Hydrogen atoms are added to ligands using Open Babel. This step prepares ligands for computational analyses that require explicit hydrogens, such as molecular docking and quantum chemical calculations.

18. **Logging of Removed Files**

    A JSON file (`removed_files.json`) is created to log all files removed at each step. This provides transparency and allows for auditing of the cleaning process.

19. **Notification**

    Upon completion of the pipeline, an email notification is sent to the user.

**Results**

Applying LPCE to the PDB dataset resulted in a curated set of protein-ligand complexes suitable for computational analysis. The pipeline processed a total of X structures, with the following outcomes:

- **Structures containing DNA/RNA removed**: Y structures
- **Structures with multiple models removed**: Z structures
- **Water molecules removed**: W instances
- **Unwanted ligands removed**: V ligands
- **Unique protein-ligand complexes obtained**: U complexes

*Note: Replace X, Y, Z, W, V, U with actual numbers obtained after processing.*

**Discussion**

LPCE provides a systematic and automated approach to cleaning and processing protein-ligand structural data. Compared to existing datasets like PDBbind, LPCE offers several advantages:

- **Up-to-date Data**: By processing the latest PDB entries, LPCE ensures that the dataset is current.

- **Automated Cleaning**: The pipeline removes unwanted components and redundant structures, reducing manual effort and potential errors.

- **Customization**: Users can adjust filtering criteria, such as distance thresholds and ligand inclusion lists, to suit specific research needs.

- **Transparency**: Logging of removed files and processing steps allows for reproducibility and auditing.

However, LPCE also has limitations:

- **Dependency on PDB Annotations**: The accuracy of the pipeline relies on the correctness of PDB annotations. Errors in PDB files can propagate through the pipeline.

- **Potential Loss of Relevant Structures**: Strict filtering criteria may inadvertently remove relevant complexes, especially if binding site information is incomplete.

- **Computational Resources**: Processing large datasets requires significant computational resources and storage.

Future work includes:

- **Enhancing Ligand Classification**: Integrating machine learning methods to improve the identification and classification of ligands.

- **Handling of Special Cases**: Extending the pipeline to process membrane proteins, multi-chain complexes, and alternative conformations.

- **User Interface Development**: Creating a user-friendly interface to allow non-expert users to customize and run the pipeline.

**Conclusion**

The Ligand Processing and Cleaning Engine (LPCE) is an effective tool for generating high-quality protein-ligand datasets. By automating the cleaning and processing of PDB structures, LPCE addresses common issues found in public datasets and facilitates more accurate computational analyses. LPCE can serve as a valuable resource for researchers in computational biology, structural bioinformatics, and drug discovery.

**Acknowledgments**

The author thanks the contributors to the LPCE project and acknowledges the support of the open-source community. Special thanks to collaborators who provided insights and feedback during the development of the pipeline.

**References**

[1] Berman, H. M., Westbrook, J., Feng, Z., et al. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235-242.

[2] Liu, Z., Li, Y., Han, L., et al. (2015). PDB-wide collection of binding data: current status of the PDBbind database. *Bioinformatics*, 31(3), 405-412.

[3] O'Boyle, N. M., Banck, M., James, C. A., et al. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3(1), 33.

[4] van Kempen, M., Kim, S. S., Tumescheit, C., Mirdita, M., Gilchrist, C. L. M., & Söding, J. (2022). Foldseek: fast and accurate protein structure search. *bioRxiv*.
