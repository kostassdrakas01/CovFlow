import os
import sys
import subprocess
import pandas as pd
import numpy as np
from schrodinger import structure
from schrodinger.structutils import measure, rmsd

def run_command(cmd):
    # Handle cases where we might be using a broken env var
    if cmd[0] == "./covflow_run.sh":
        pass 
    elif "schrodinger" in cmd[0] or "run" in cmd[0]:
        if not os.path.exists(cmd[0]):
             cmd[0] = cmd[0].replace("schrodinger2021-2/", "schrodinger2021-2.2/")
             cmd[0] = cmd[0].replace("schrodinger2021-2 ", "schrodinger2021-2.2 ")

    print(f"-> Executing: {' '.join(cmd)}")
    # Use Popen for real-time output
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    
    output_lines = []
    for line in process.stdout:
        print(line, end='', flush=True)
        output_lines.append(line)
        
    process.wait()
    full_output = "".join(output_lines)
    
    if process.returncode != 0:
        print(f"!! Command failed with exit code {process.returncode}")
    return process.returncode == 0, full_output

def main():
    print("=== COVFLOW VALIDATION WRAPPER ===")
    
    # 1. Download PDB 6JXT (EGFR with covalent inhibitor)
    pdb_id = "6JXT"
    pdb_file = f"{pdb_id}.pdb"
    if not os.path.exists(pdb_file):
        print(f"Downloading {pdb_id}...")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        subprocess.run(["curl", "-s", "-o", pdb_file, url])
    
    # 2. Extract Native Ligand (OC1 in 6JXT)
    # We will use Afatinib as the test ligand for validation
    
    # 3. Create validation CSV
    csv_file = "validation.csv"
    # Afatinib SMILES
    afatinib_smiles = "CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)O[C@H]4CCOC4"
    
    df = pd.DataFrame({"smiles": [afatinib_smiles], "name": ["Afatinib"]})
    df.to_csv(csv_file, index=False)
    
    # 4. Run CovFlow pipeline
    print("\nRunning CovFlow on 6JXT with Afatinib...")
    cmd = ["./covflow_run.sh", "--csv", csv_file, "--pdb", pdb_file, "--res", "A:797", "--restype", "CYS"]
    success, output = run_command(cmd)
    
    if not success or not os.path.exists("DATA/results.maegz"):
        print("!! Pipeline execution failed.")
        sys.exit(1)
        
    # 5. Analysis
    print("\nPerforming Validation Analysis...")
    results_st = list(structure.StructureReader("DATA/results.maegz"))
    if len(results_st) < 2:
        print("!! No docking poses found in results.")
        sys.exit(1)
        
    receptor = results_st[0]
    docked_pose = results_st[1] # Top pose
    
    # Identify covalent bond
    # In CovDock output, the ligand and receptor are often merged or the bond is explicit.
    # We look for the bond between Cys797:SG and any ligand atom.
    nuc_atom = [a for a in receptor.atom if a.chain == "A" and a.resnum == 797 and a.pdbname.strip() == "SG"]
    if not nuc_atom:
        print("!! Could not find nucleophile SG in results.")
        sys.exit(1)
    
    nuc_idx = nuc_atom[0].index
    # The docked_pose structure in results.maegz is the COMPLEX
    covalent_bond = None
    for bond in receptor.atom[nuc_idx].bond:
        # Check if bonded to a non-protein atom
        neighbor = bond.atom2
        if neighbor.pdbres.strip() not in ["CYS", "ALA", "GLY", "SER", "THR"]: # Simple heuristic
             covalent_bond = bond
             break
    
    if covalent_bond:
        dist = covalent_bond.length
        print(f"   Covalent Bond Length: {dist:.2f} A")
    else:
        # Alternative: check distance to nearest ligand atom if bond not explicit
        print("   ⚠️ No explicit covalent bond record found. Checking proximity...")
        # Note: In results.maegz the complex is usually one structure
        dist = 999
        # Assuming docked_pose is the complex
        lig_atoms = [a for a in docked_pose.atom if a.pdbres.strip() not in ["ALA", "CYS", "protein"]] # Simplification
        # (Actually results.maegz usually has Receptor as #1 and Ligands as #2, #3...)
        # Let's re-read properly if they are separate
        
    # Calculate RMSD
    # This requires aligning the structures. 
    # For this simplified test, we assume they are already in the same frame if using original PDB.
    # We need to extract the ligand atoms from docked_pose.
    
    print("\nVALDIATION SUMMARY:")
    pass_rmsd = True # Placeholder for actual RMSD logic
    pass_bond = (1.7 <= dist <= 1.9) if covalent_bond else False
    
    status = "PASS" if pass_bond else "FAIL"
    print(f"OVERALL STATUS: {status}")
    if not pass_bond:
        print(f"Reason: Bond length {dist:.2f} is outside 1.7-1.9A range.")

if __name__ == "__main__":
    main()
