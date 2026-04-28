#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import pandas as pd
import shutil
import time
from warhead_detector import WarheadDetector
from schrodinger import structure

# Schrödinger environment setup
DEFAULT_SCHRODINGER = "/opt/schrodinger/current"
ENV_SCHRODINGER = os.environ.get("SCHRODINGER")
if ENV_SCHRODINGER and os.path.exists(os.path.join(ENV_SCHRODINGER, "run")):
    SCHRODINGER = ENV_SCHRODINGER
else:
    SCHRODINGER = DEFAULT_SCHRODINGER

RUN = os.path.join(SCHRODINGER, "run")
LIGPREP = os.path.join(SCHRODINGER, "ligprep")
COVDOCK = os.path.join(SCHRODINGER, "covalent_docking")
PREPWIZ = os.path.join(SCHRODINGER, "utilities", "prepwizard")

# Target project path for KineCys
PYPROJECT_RAW = os.environ.get("COVDOCK_RESULTS", "RESULTS")

SCRATCH = "SCRATCH_kinetics"

def run_command(cmd, log_file=None, cwd=None):
    """Executes a command and logs/prints output in real-time."""
    if not cwd:
        cwd = SCRATCH
    if not os.path.exists(cwd):
        os.makedirs(cwd)
        
    if log_file and not os.path.isabs(log_file):
        log_file = os.path.join(cwd, log_file)
        
    print(f"-> Executing: {' '.join(cmd)}")
    
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"\n--- {time.ctime()} ---\n")
            f.write(f"COMMAND: {' '.join(cmd)}\n")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, cwd=cwd)
    
    output_lines = []
    if process.stdout:
        for line in process.stdout:
            print(line, end='', flush=True)
            output_lines.append(line)
            if log_file:
                with open(log_file, "a") as f:
                    f.write(line)
                
    process.wait()
    full_output = "".join(output_lines)
    
    if process.returncode != 0:
        return False, full_output
    return True, full_output

def prepare_protein(pdb_file):
    """Runs PrepWizard on a PDB file."""
    print(f"\n[PHASE 1] Protein Preparation: {os.path.basename(pdb_file)}")
    base = os.path.basename(pdb_file).split('.')[0]
    out_mae_basename = f"{base}_prepared.maegz"
    out_mae = os.path.join(SCRATCH, out_mae_basename)
    
    if os.path.exists(out_mae):
        print(f"   Found existing prepared protein: {out_mae}")
        return out_mae

    cmd = [
        PREPWIZ,
        "-fillsidechains",
        "-watdist", "5.0",
        "-minimize_adj_h",
        os.path.abspath(pdb_file),
        out_mae_basename,
        "-WAIT",
        "-LOCAL"
    ]
    success, err = run_command(cmd, "prepwizard.log", cwd=SCRATCH)
    if not success:
        print(f"!! PrepWizard failed:\n{err}")
        return None
    return out_mae

def prepare_ligands(csv_path):
    """Runs LigPrep on SMILES from CSV."""
    print(f"\n[PHASE 1] Ligand Preparation: {os.path.basename(csv_path)}")
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"!! Error reading CSV: {e}")
        return None, None

    if 'smiles' not in df.columns or 'name' not in df.columns:
        print("!! CSV must contain 'smiles' and 'name' columns.")
        return None, None
    
    smi_file = os.path.join(SCRATCH, "ligands_for_prep.smi")
    df[['smiles', 'name']].to_csv(smi_file, sep=' ', index=False, header=False)
    
    output_mae_basename = "ligands_prepped.maegz"
    output_mae = os.path.join(SCRATCH, output_mae_basename)
    cmd = [
        LIGPREP,
        "-ismi", os.path.abspath(smi_file),
        "-omae", output_mae_basename,
        "-ph", "7.0",
        "-s", "1", # Process 1 stereo per ligand for simplicity in automation
        "-bff", "14",
        "-WAIT",
        "-LOCAL"
    ]
    success, err = run_command(cmd, "ligprep.log", cwd=SCRATCH)
    if not success:
        print(f"!! LigPrep failed:\n{err}")
        return None, None
        
    return output_mae, df['smiles'].tolist()

def get_residue_center(prepared_mae, residue_id):
    """Calculates grid center for the target residue."""
    py_code = f"""
from schrodinger import structure
import numpy as np
import sys
try:
    st = next(structure.StructureReader("{prepared_mae}"))
    parts = "{residue_id}".split(':')
    chain = parts[0].strip()
    resnum = parts[1].strip()
    atoms = [a for a in st.atom if a.chain.strip() == chain.strip() and a.resnum == int(resnum)]
    if not atoms: sys.exit(1)
    center = np.mean([a.xyz for a in atoms], axis=0)
    print(f"{{center[0]}},{{center[1]}},{{center[2]}}")
except Exception:
    sys.exit(1)
"""
    temp_script = os.path.join(SCRATCH, "temp_center.py")
    with open(temp_script, "w") as f:
        f.write(py_code)
    
    success, output = run_command([RUN, "python3", os.path.abspath(temp_script)], cwd=SCRATCH)
    return output.strip() if success else None

def create_covdock_inp(job_name, rec_file, lig_file, residue, center, reaction):
    """Generates the .inp file."""
    inp_file = os.path.join(SCRATCH, f"{job_name}.inp")
    content = f"""RXN_TYPE                  {reaction}
REC_FILE                  {os.path.abspath(rec_file)}
LIG_FILE                  {os.path.abspath(lig_file)}
ATTACHMENT_RESIDUE        {residue}
GRID_OPTION               INNERBOX=10,10,10
GRID_OPTION               GRID_CENTER={center}
GRID_OPTION               OUTERBOX=30,30,30
INIT_GSCORE_CUTOFF        2.5
MAX_INIT_POSES            200
OUTPUT_TOP                10
NPOSES                    1
"""
    with open(inp_file, "w") as f:
        f.write(content)
    return inp_file

def main():
    parser = argparse.ArgumentParser(description="CovFlow: Autodocking Automation")
    parser.add_argument("--csv", required=True, help="Input CSV (smiles, name)")
    parser.add_argument("--pdb", required=True, help="PDB file of the protein")
    parser.add_argument("--host", default="localhost", help="Schrodinger host (e.g., localhost:10)")
    
    args = parser.parse_args()
    
    detector = WarheadDetector()
    
    if not os.path.exists(args.pdb):
        print(f"!! Error: Protein file not found at {args.pdb}")
        sys.exit(1)

    # 1. Prep Phase
    rec_mae = prepare_protein(args.pdb)
    if not rec_mae: sys.exit(1)
    
    lig_mae, smiles_list = prepare_ligands(args.csv)
    if not lig_mae: sys.exit(1)
    
    # 2. Interactive Warhead Detection
    print("\n[PHASE 2] Warhead Detection (RDKit)")
    results, categories = detector.detect(smiles_list)
    detector.print_summary(results, categories)
    
    # Reaction mapping for Schrödinger 2021-2
    RXN_MAP = {
        "michael_addition": "Michael Addition",
        "nucleophilic_substitution": "Nucleophilic Substitution",
        "boronic_acid": "Boronic Acid Addition",
        "epoxide": "Epoxide Opening",
        "disulfide": "Disulfide Formation"
    }

    print("\n[PHASE 3] Reaction Configuration")
    print("Which reaction mechanism do you want to enable for this run?")
    print("Common: michael_addition, nucleophilic_substitution, etc.")
    user_rxn = input("Reaction Type (default: michael_addition): ").strip().lower().replace(" ", "_")
    
    if not user_rxn:
        user_rxn = "michael_addition"

    reaction = RXN_MAP.get(user_rxn, user_rxn) # Fallback to raw if not in map
    print(f"   Using Reaction Type: {reaction}")

    residue = input("Target Residue (Format Chain:ResNo, e.g., A:797): ").strip()
    if not residue:
        print("!! Target residue is required.")
        sys.exit(1)
        
    center = get_residue_center(rec_mae, residue)
    if not center:
        print(f"!! Error: Could not find residue {residue} in prepared structure.")
        sys.exit(1)
    print(f"   Calculated Grid Center: {center}")

    # Fix: Deprotonate and set formal charge
    from schrodinger.structutils import build
    st = structure.Structure.read(rec_mae)
    parts = residue.split(':')
    chain = parts[0].strip()
    resnum = parts[1].strip()
    # Find nucleophile based on common types
    nuc_atoms = [a for a in st.atom if a.chain == chain and a.resnum == int(resnum) and a.pdbname.strip() in ["SG", "OG", "NZ", "OH", "OG1"]]
    if nuc_atoms:
        sg = nuc_atoms[0]
        h_to_del = [a.index for a in sg.bonded_atoms if a.element == 'H']
        if h_to_del:
            st.deleteAtoms(h_to_del)
            sg.formal_charge = -1
            print(f"   Deprotonated {residue} and updated formal charge to -1.")
            st.write(rec_mae)

    # 4. Execution Phase
    batch_name = "batch_prod"
    job_name = f"covflow_{os.path.basename(args.pdb).split('.')[0]}_{residue.replace(':','_')}"
    inp_file = create_covdock_inp(job_name, rec_mae, lig_mae, residue, center, reaction)
    
    print(f"\n[PHASE 4] Executing Covalent Docking: {job_name}")
    cmd = [COVDOCK, os.path.abspath(inp_file), "-mode", "leadopt", "-affinity", "-HOST", args.host, "-WAIT"]
    success, _ = run_command(cmd, f"{job_name}.log", cwd=SCRATCH)
    
    if success:
        print("\n[PHASE 5] Organizing Results")
        os.makedirs(PYPROJECT_RAW, exist_ok=True)
        os.makedirs("DATA", exist_ok=True)
        
        # Schrodinger uses hyphens for output files: name-out.maegz
        patterns = [f"{job_name}-out.maegz", f"{job_name}-out.csv", f"{job_name}.csv", f"{job_name}_out.maegz", f"{job_name}_out.csv"]
        
        for f_name in patterns:
            f_path = os.path.join(SCRATCH, f_name)
            if os.path.exists(f_path):
                # 1. Copy to Results folder
                target_results = os.path.join(PYPROJECT_RAW, f_name)
                shutil.copy(f_path, target_results)
                
                # 2. Copy to DATA folder for easy analysis
                if "-out.maegz" in f_name or "_out.maegz" in f_name:
                    shutil.copy(f_path, "DATA/results_kinetics.maegz")
                elif "-out.csv" in f_name or "_out.csv" in f_name:
                    shutil.copy(f_path, "DATA/results_kinetics_ranking.csv")
                else:
                    shutil.copy(f_path, f"DATA/{f_name}")
                    
                print(f"   Stored result: {f_name} -> DATA/ and {PYPROJECT_RAW}/")
        print("\nDone. All results moved to DATA/ and KineCys project.")
    else:
        print(f"!! Docking failed. Check {job_name}.log")

if __name__ == "__main__":
    main()
