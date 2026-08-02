Automated High-Throughput Docking & MM-GBSA Validation Pipeline
📖 Overview
Relying solely on raw molecular docking scores often leads to false positives due to unnatural ligand strain and the limitations of rigid receptor environments. The "paused scene" of a high docking score can be highly deceiving.

This fully automated, Linux-native computational biology pipeline bridges the gap between high-throughput molecular docking and Molecular Dynamics (MD). It processes raw PDB structures and ligands, performs GPU-accelerated flexible docking, applies strict biophysical filters (H-bond validation & conformational strain), and validates the surviving candidates using MM-GBSA thermodynamic calculations via AmberTools.

✨ Key Features
End-to-End Automation: Provide a list of PDB IDs and ligand names; the pipeline handles download, sanitization, topology generation, docking, and scoring automatically.

Dynamic Flexibility & Pocket Detection: Automatically detects bound ligands (small molecules or peptides), calculates grid box dimensions, and dynamically selects the top interacting residues for flexible docking.

Strict Biophysical Filtering: Rejects "false positive" docking poses by calculating precise H-bond distances (Bio.PDB) and evaluating post-docking conformational strain using Meeko and RDKit (MMFF94).

GPU-Accelerated Docking: Integrates Uni-Dock for blazing-fast, CUDA-accelerated flexible docking.

Seamless AMBER Integration: Features an intelligent Protein Flexibility Updater that safely re-inserts TER gap-breaks into flexible structures, preventing tleap crashes during topology generation.

Smart File-System Validation: A built-in "Smart Checker" verifies the physical output of every single stage. If a stage fails, the pipeline halts to prevent cascading errors.

Data Isolation & Archiving: Automatically extracts the true target name from PDB metadata, creates isolated project directories, moves all relevant files, and purges intermediate junk.
🛠️ Pipeline ArchitectureThe master pipeline (master.py) orchestrates a sequential execution of the following modules:StageScriptDescription0install_tools.shLocal installation script for ADFR Suite and Uni-Dock binaries (does not modify ~/.bashrc).1clean_pro.pyDownloads raw PDBs, identifies native ligands/peptides, calculates Grid Box parameters, isolates chains, and prepares rigid/flex .pdbqt files via ADFR.2ligand.pyFetches 3D structures from PubChem, performs energy minimization (MMFF94), detects rotatable bonds, and prepares .pdbqt files via Meeko.3flex_docking.pyExecutes GPU-accelerated flexible molecular docking using Uni-Dock based on dynamically generated config files.4filter.pyThe Gatekeeper: Parses docked poses, extracts affinity, checks for specific H-bonds, calculates ligand strain, and ranks by SILE (Size-Independent Ligand Efficiency).5rec_pdb.pyRebuilds the flexible receptor structure, updating coordinates for flex-residues and intelligently fixing gap-breaks (TER tags) for AMBER compatibility.6amber.pyPrepares complex topologies (tleap, antechamber, parmchk2) assigning AM1-BCC charges for ligands and FF14SB/GAFF force fields.7gbsa.pyPerforms structure relaxation (Minimization), calculates MM-GBSA binding free energies (MMPBSA.py), executes per-residue decomposition, and generates final CSV reports.
⚙️ Installation & Prerequisites
This pipeline is designed for Linux environments. Ensure you have the CUDA toolkit installed for Uni-Dock acceleration.

1. Setup Conda Environments
To prevent dependency conflicts, it is highly recommended to use two separate environments:

Main Docking Environment:

Bash
conda create -n docking_env python=3.10
conda activate docking_env
conda install -c conda-forge rdkit biopython pandas numpy openbabel
pip install meeko
AmberTools Environment:

Bash
conda create -n amber_env -c conda-forge ambertools
(Ensure the path to amber_env in rec_pdb.py and amber.py matches your local conda installation path).

2. Install Local Binaries
Run the included bash script to download and configure ADFR Suite and Uni-Dock locally:

Bash
chmod +x install_tools.sh
./install_tools.sh
🚀 Usage
Configure Targets: Edit the pdb_ids list in clean_pro.py with your desired targets.

Configure Ligands: Edit the ligands_list in ligand.py.

Run the Pipeline:

Bash
conda activate docking_env
python master.py
The master script will execute all stages sequentially, providing live, color-coded terminal feedback.
📊 Outputs & ReportsOnce the pipeline completes, it dynamically creates a final_results/[TARGET_NAME]/ directory. Intermediate folders (raw_pdb, ligand, clean_proteins) are automatically nuked to save space.Inside the Final_Outputs directory, you will find two primary reports:1. Docking Summary (flex_docking_summary.csv)Provides a rapid overview of the initial docking phase, detailing the detected flexible residues and the raw Vina affinities.Example output for target 8PQ9:Protein_NameLigandFlexible_ResiduesBest_Affinity_kcal/molMAST/STEM CELL GROWTH FACTOR RECEPTOR KITavapritinibCYS_673, LEU_625, MET_624...-11.170MAST/STEM CELL GROWTH FACTOR RECEPTOR KITnilotinibCYS_673, LEU_625, MET_624...-10.940MAST/STEM CELL GROWTH FACTOR RECEPTOR KITcholesterolCYS_673, LEU_625, MET_624...-9.306
2. Thermodynamic Validation (MMGBSA_Final_Ranking.csv)
The ultimate source of truth. Candidates that survived the biophysical strain filter are ranked by their true thermodynamic stability, including Per-Residue Decomposition data to highlight critical binding hotspots.
Example output:

System	Delta_G (kcal/mol)	Top_5_Interacting_Residues
8PQ9_avapritinib	-65.29	CYS673 (-4.2), LEU625 (-3.1)...
8PQ9_nilotinib	-42.71	CYS673 (-3.8), GLU671 (-2.5)...
Candidates performing exceptionally well in this final CSV are highly stable and ready to be exported to GROMACS for 100ns+ production MD simulations.
