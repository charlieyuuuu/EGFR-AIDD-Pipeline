# AI-Assisted Virtual Screening & MD Validation Pipeline for EGFR Inhibitors

## 📌 Project Overview
This repository presents an end-to-end AIDD/CADD pipeline for discovering novel EGFR inhibitors. The workflow integrates **structure-based virtual screening**, **machine learning (QSAR) cross-validation**, and **100 ns all-atom molecular dynamics (MD) simulations**. The primary goal is to eliminate false positives from traditional docking and evaluate the dynamic binding stability of candidate complexes.

## 🛠️ Tech Stack
*   **AI & Cheminformatics:** Scikit-learn (Random Forest), RDKit, Pandas (Morgan Fingerprints)
*   **Structure-Based Design (CADD):** AutoDock Vina (v1.2.5), Meeko, OpenBabel
*   **Molecular Dynamics:** GROMACS 2018 (GPU accelerated), acpype (GAFF2 force field)
*   **Automation & Visualization:** Linux Shell, Python `subprocess`, Matplotlib, PyMOL

## 📂 Pipeline Overview
1.  **Data Preparation**: Filtered 10,734 EGFR-active compounds from the ChEMBL database using Lipinski's Rule of Five via RDKit. Generated 3D conformers and PDBQT formats.
2.  **High-Throughput Virtual Screening (CADD)**: Automated the docking of ~7,500 compounds against the EGFR kinase domain using Python `subprocess` scripts. Selected the Top 50 candidates based on binding affinity.
3.  **AI/QSAR Cross-Validation (AIDD)**: Extracted 2048-bit ECFP4 Morgan fingerprints and trained a Random Forest Regressor (Test R² = 0.682) to predict bioactivity (pIC50). This ML model was used to cross-validate the Top 50 docking hits.
4.  **MD Setup & Simulation (CADD)**: Conducted 100 ns MD simulations on 6 representative complexes using GROMACS. Developed automated scripts for topology merging (handling GAFF2 parameters) and utilized GPU acceleration for high-performance computing.
5.  **Trajectory Analysis**: Evaluated complex stability by calculating ligand RMSD, protein backbone RMSD, and minimum distances using `gmx trjconv`, `gmx rms`, and `gmx mindist`.

## 🔬 Key Findings (AIDD + CADD Synergy)
*Static docking and pure ML predictions have limitations. By combining AI and physics-based methods, we revealed a critical "false positive" case:*
*   The top-ranked docking pose of **CHEMBL6001037** (Vina score: -11.57 kcal/mol) **dissociated** from the binding pocket within 3 ns during the MD simulation.
*   However, a secondary docking pose (pose 2) remained **highly stable** throughout the 100 ns trajectory (ligand RMSD ~ 0.4 nm).
*   *Conclusion:* This highlights the necessity of dynamic validation over static docking scores. The integration of QSAR filtering and MD simulations is essential in modern AIDD workflows to ensure the reliability of hit identification.

## 📊 Results
![Ligand RMSD](results/Fig1_Ligand_RMSD.png)
*Figure 1: Ligand RMSD over 100 ns MD simulations for 6 candidate complexes, demonstrating varying dynamic stabilities.*

## 🚀 Quick Start
1.  Install dependencies: `pip install -r requirements.txt`
2.  Prepare ligands & filter data: `python 01_data_prep.py`
3.  Run high-throughput docking: `python 02_batch_docking.py`
4.  Train QSAR model: `python 03_qsar_model.py`
5.  Analyze MD trajectories & plot: `python 04_md_analysis.py`
