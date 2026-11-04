"""
05_admet_mmpbsa.py
Description:
    This script performs ADMET property prediction for the Top 50 candidates
    using the 'admet-ai' package. It also provides the workflow instructions
    for calculating binding free energy (MM-PBSA) via 'gmx_MMPBSA'.
"""

import os
import pandas as pd

def run_admet_prediction():
    print("=" * 50)
    print("Starting ADMET Prediction Pipeline...")
    print("=" * 50)

    os.environ["CUDA_VISIBLE_DEVICES"] = ""

    try:
        from admet_ai import ADMETModel
    except ImportError:
        print("Error: 'admet-ai' is not installed. Please run: pip install admet-ai")
        return

        print("Loading ADMET Model...")
    model = ADMETModel()

    print("Reading Top 50 docking results and ChEMBL database...")
    top50_path = '../top50_docking.csv' if os.path.exists('../top50_docking.csv') else 'top50_docking.csv'
    raw_db_path = '../CHEMBL203_actives.csv' if os.path.exists('../CHEMBL203_actives.csv') else 'CHEMBL203_actives.csv'

    top50 = pd.read_csv(top50_path)
    raw_db = pd.read_csv(raw_db_path)

    top50 = top50.merge(raw_db[['chembl_id', 'smiles']], on='chembl_id', how='left')
    top50 = top50.dropna(subset=['smiles'])

    print(f"Successfully matched {len(top50)} molecules. Predicting...")
    results = model.predict(smiles=top50['smiles'].tolist())
    
    core_props = results[['hERG', 'Caco2_Wang', 'CYP3A4_Veith', 'Half_Life_Obach']]
    final_df = pd.concat([top50[['chembl_id', 'affinity', 'smiles']], core_props], axis=1)

    out_csv = 'top50_ADMET_evaluated.csv'
    final_df.to_csv(out_csv, index=False)
    print(f"✅ ADMET prediction completed! Saved to {out_csv}")



def print_mmpbsa_instructions():
    print("\n" + "=" * 50)
    print("MM-PBSA Binding Free Energy Calculation Workflow")
    print("=" * 50)
    print("""
To calculate the binding free energy (ΔG) from your 100 ns MD trajectory, 
we use the AMBER-based 'gmx_MMPBSA' tool. Follow these steps:

1. Create a separate conda environment:
   $ conda create -n mmpbsa_env -c conda-forge gmx_mmpbsa python=3.10 -y
   $ conda activate mmpbsa_env

2. Prepare an input file 'mmpbsa.in' (example for GB model):
   -----------------------------------------------------
   &general
     sys_name="EGFR_complex",
     startframe=9900,  # Calculate the last 1 ns of the 100 ns trajectory
     endframe=10000,
     interval=1,
     temperature=298.15,
     PBRadii=4,
   /
   &gb
     igb=5, saltcon=0.150,
   /
   -----------------------------------------------------

3. Run the calculation (ensure files are in the same directory):
   $ gmx_MMPBSA -O -i mmpbsa.in -cs md_pose2.tpr -ci index.ndx \\
                -cg 1 13 -ct md_pose2_center.xtc -cp topol_pose2.top \\
                -o FINAL_RESULTS.dat

4. View the results:
   $ grep "DELTA TOTAL" FINAL_RESULTS.dat

   *Expected Result:* ΔG ≈ -37 kcal/mol (Van der Waals and Electrostatics 
   are the main driving forces for binding).
    """)


if __name__ == "__main__":
    # Run ADMET prediction
    run_admet_prediction()

    # Print MM-PBSA usage guide
    print_mmpbsa_instructions()
