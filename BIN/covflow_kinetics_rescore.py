#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import subprocess
from schrodinger import structure

def run_kinetics_rescore(results_mae, protein_mae, residue, reaction, host="localhost"):
    """
    Runs kinetics analysis on already docked poses using Pose Prediction mode.
    """
    job_name = "kinetics_rescore"
    SCRATCH = "SCRATCH_kinetics"
    if not os.path.exists(SCRATCH): os.makedirs(SCRATCH)
    
    # 1. Extract ligands from results if they are complexes
    print(f"-> Extracting ligands from {results_mae}...")
    ligand_mae = os.path.join(SCRATCH, "extracted_ligands.maegz")
    with structure.StructureWriter(ligand_mae) as writer:
        for st in structure.StructureReader(results_mae):
            # If it's a complex, we extract the ligand. 
            ligs = st.extract([a.index for a in st.atom if a.pdbres.strip() not in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP", "CYX"]])
            if ligs.atom_total > 0:
                writer.append(ligs)
            else:
                writer.append(st)

    # 2. Create Input File
    inp_file = os.path.join(SCRATCH, f"{job_name}.inp")
    content = f"""RXN_TYPE                  {reaction}
REC_FILE                  {os.path.abspath(protein_mae)}
LIG_FILE                  {os.path.abspath(ligand_mae)}
ATTACHMENT_RESIDUE        {residue}
OUTPUT_TOP                100
"""
    with open(inp_file, "w") as f:
        f.write(content)
    
    # 3. Execute
    print(f"-> Executing KineCys Rescoring (Score Mode)...")
    log_file = os.path.join(SCRATCH, f"{job_name}.log")
    cmd = [os.path.join(os.environ['SCHRODINGER'], "covalent_docking"), os.path.abspath(inp_file), "-mode", "score", "-HOST", host, "-WAIT"]
    
    with open(log_file, "w") as log_f:
        subprocess.run(cmd, cwd=SCRATCH, stdout=log_f, stderr=subprocess.STDOUT)
    
    print("\nRescoring Complete.")
    # Check for hyphens in output name as per Schrodinger convention
    out_file_pattern = f"{job_name}-out.maegz"
    out_file_scratch = os.path.join(SCRATCH, out_file_pattern)
    
    if os.path.exists(out_file_scratch):
        # Copy to DATA for easy access
        shutil.copy(out_file_scratch, "DATA/kinetics_results.maegz")
        print(f"Results saved to: DATA/kinetics_results.maegz")
    else:
        print(f"!! Rescoring failed. Check {log_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run kinetics on already docked poses.")
    parser.add_argument("--results", required=True, help="MAE file with docked poses")
    parser.add_argument("--pdb", required=True, help="Original prepared protein PDB/MAE")
    parser.add_argument("--res", required=True, help="Target Residue (e.g., A:797)")
    parser.add_argument("--rxn", default="Michael Addition", help="Reaction type")
    
    args = parser.parse_args()
    
    if 'SCHRODINGER' not in os.environ:
        print("!! Error: SCHRODINGER environment variable not set.")
        sys.exit(1)
        
    run_kinetics_rescore(args.results, args.pdb, args.res, args.rxn)
