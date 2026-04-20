#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import pandas as pd
import shutil
import time
import numpy as np

# Schrodinger imports (assumes running via $SCHRODINGER/run)
try:
    from schrodinger import structure
    from schrodinger.structutils import analyze, build, minimize
except ImportError as e:
    print(f"!! Error: {e}")
    sys.exit(1)

from warhead_detector import WarheadDetector

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
GLIDE = os.path.join(SCHRODINGER, "glide")

# Target project path (Results synchronization)
PYPROJECT_RAW = os.environ.get("COVDOCK_RESULTS", "RESULTS")
SCRATCH = "SCRATCH" # Default

COFACTORS = ["HEM", "MG", "CA", "ZN", "CL", "SO4", "PO4", "HOH", "EDT", "DMS"]
STANDARD_RESIDUES = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP", "CYX"
]
NUCLEOPHILES = {
    "CYS": "SG",
    "SER": "OG",
    "LYS": "NZ",
    "TYR": "OH",
    "THR": "OG1"
}

def run_command(cmd, log_file=None):
    if log_file and not os.path.isabs(log_file):
        log_file = os.path.join(SCRATCH, log_file)
    print(f"-> Executing: {' '.join(cmd)}")
    
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"\n--- {time.ctime()} ---\n")
            f.write(f"COMMAND: {' '.join(cmd)}\n")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    output_lines = []
    for line in process.stdout:
        print(line, end='', flush=True)
        output_lines.append(line)
        if log_file:
            with open(log_file, "a") as f:
                f.write(line)
    process.wait()
    return process.returncode == 0, "".join(output_lines)

def sanitize_protein(pdb_file, target_res_id):
    """
    1. Removes native ligands.
    2. Restores covalent residues.
    3. Calculates grid center.
    """
    print(f"\n[NEW PHASE] Sanitizing Protein & Restoring Bonds")
    # Use StructureReader for better tolerance of problematic PDBs
    with structure.StructureReader(pdb_file) as reader:
        st = next(reader)
    chain, resnum = target_res_id.split(':')
    chain = chain.strip()
    resnum = resnum.strip()
    
    # 1. Find the target residue
    target_atoms = [a for a in st.atom if a.chain == chain and str(a.resnum) == resnum]
    if not target_atoms:
        print(f"!! Error: Could not find residue {target_res_id}")
        return None, None
    
    target_res = target_atoms[0].getResidue()
    
    # 2. Check for existing covalent bonds
    covalent_ligand_atoms = []
    target_atom_objs = [st.atom[i] for i in target_atoms]
    bonds_broken = False
    for atom_obj in target_atom_objs:
        for bond in atom_obj.bond:
            neighbor = bond.atom2
            n_resname = neighbor.pdbres.strip()
            if n_resname not in STANDARD_RESIDUES and n_resname not in COFACTORS:
                print(f"   Detected covalent bond to {n_resname} ({neighbor.resnum})")
                st.deleteBond(atom_obj, neighbor)
                # Restore Hydrogen based on residue
                build.add_hydrogens(st, atom_list=[atom_obj])
                bonds_broken = True

    if bonds_broken:
        print(f"   Broken bonds detected. Locally minimizing target residue to restore native state.")
        # We handle full structure minimization below
    
    # 3. Native Ligand Identification (within 10A)
    # We want to identify the pocket center, but prioritize the nucleophile for the grid.
    target_atom_objs = [st.atom[i] for i in target_atoms]
    coords = [a.xyz for a in target_atom_objs]
    centroid = np.mean(coords, axis=0)
    
    # Always use the Nucleophile (e.g., Cys SG) as the Grid Center for Covalent Docking
    nuc_atom_name = NUCLEOPHILES.get(target_res.pdbres.strip(), "SG")
    nuc_atom_list = [a for a in target_atom_objs if a.pdbname.strip() == nuc_atom_name]
    
    if nuc_atom_list:
        center_coords = nuc_atom_list[0].xyz
        center = f"{center_coords[0]:.3f},{center_coords[1]:.3f},{center_coords[2]:.3f}"
        print(f"   Grid Center (Nucleophile): {center}")
    else:
        # Emergency Fallback to geometric center of residue
        center = f"{centroid[0]:.3f},{centroid[1]:.3f},{centroid[2]:.3f}"
        print(f"   Grid Center (Residue Centroid): {center}")

    # Remove all native ligands/cofactors to avoid clashes
    to_delete = []
    for res in st.residue:
        res_name = res.pdbres.strip()
        if res_name not in COFACTORS and res_name not in STANDARD_RESIDUES:
            to_delete.extend(res.getAtomIndices())
            
    if to_delete:
        st.deleteAtoms(to_delete)
        # Clear existing hydrogens and re-add to fix dangling bonds
        build.delete_hydrogens(st)
        build.add_hydrogens(st)
        print(f"   Deleted {len(to_delete)} atoms from native ligands and restored hydrogens.")
    
    # 4. Local Minimization (Target + 5A) to restore geometry
    try:
        # ASL: Target residue and neighbors within 5.0 A
        minimize_asl = f"(res.n {resnum} and chain {chain}) or (within 5.0 (res.n {resnum} and chain {chain}))"
        minimize.minimize_structure(st, iteration=50, asl=minimize_asl)
    except Exception as e:
        print(f"   ⚠️ Minimization warning (skipping): {e}")
    
    out_file = os.path.join(SCRATCH, "sanitized_protein.maegz")
    st.write(out_file)
    return out_file, center

def prepare_ligands_with_filter(csv_path):
    """LigPrep + Warhead Filtering."""
    print(f"\n[NEW PHASE] Warhead-Aware Ligand Preparation")
    df = pd.read_csv(csv_path)
    detector = WarheadDetector()
    
    valid_smiles = []
    valid_names = []
    
    for _, row in df.iterrows():
        smi = row['smiles']
        name = row['name']
        
        # Check for warhead
        res, cats = detector.detect([smi])
        if not res:
            print(f"   ⚠️ Skipping {name}: No warhead detected.")
            continue
            
        valid_smiles.append(smi)
        valid_names.append(name)

    if not valid_smiles:
        print("!! No valid ligands with warheads found.")
        return None

    # Run LigPrep on valid ligands
    smi_file = os.path.join(SCRATCH, "filtered_ligands.smi")
    with open(smi_file, "w") as f:
        for s, n in zip(valid_smiles, valid_names):
            f.write(f"{s} {n}\n")
    
    output_mae = os.path.join(SCRATCH, "ligands_prepped.maegz")
    output_basename = "ligands_prepped.maegz"
    cmd = [LIGPREP, "-ismi", smi_file, "-omae", output_basename, "-ph", "7.0", "-s", "1", "-WAIT", "-LOCAL"]
    success, _ = run_command(cmd, "ligprep.log")
    
    if os.path.exists(output_basename):
        shutil.move(output_basename, output_mae)
    
    if success and os.path.exists(output_mae):
        # Post-process MAE to label reactive carbon (simple heuristic for acrylamides)
        st_list = list(structure.StructureReader(output_mae))
        for st in st_list:
            # Re-locate warhead and pick the reactive atom
            # For this demo, we label the first atom of the SMARTS match
            for name, pattern_data in detector.patterns.items():
                match = analyze.evaluate_smarts_canvas(st, pattern_data["smarts"])
                if match:
                    # Heuristic: the first atom in Michael Addition smarts is usually the reactive beta carbon
                    reactive_idx = match[0][0]
                    st.property['i_user_reactive_atom'] = reactive_idx
                    break
        with structure.StructureWriter(output_mae) as writer:
            for st in st_list:
                writer.append(st)
        return output_mae
    return None

def run_alanine_scan(rec_mae, lig_mae, center, target_res_id):
    """
    1. Mutate to ALA.
    2. Gridgen.
    3. Glide SP.
    4. Return list of ligands passing -5.0.
    """
    print(f"\n[NEW PHASE] Alanine-Scan Pre-Dock Filter")
    st = structure.Structure.read(rec_mae)
    chain, resnum = target_res_id.split(':')
    
    # 1. Mutate to ALA
    # Mutations are performed on the structure using build.mutate
    target_atoms = [a for a in st.atom if a.chain == chain and a.resnum == int(resnum)]
    if not target_atoms:
        print(f"!! Error: Residue {target_res_id} not found for mutation.")
        return None
    
    # Mutate the residue containing the first atom found
    build.mutate(st, target_atoms[0].index, "ALA")
    out_mutant = os.path.join(SCRATCH, "protein_ala.maegz")
    st.write(out_mutant)
    
    # 2. Gridgen
    grid_inp = os.path.join(SCRATCH, "gridgen.inp")
    grid_file = os.path.join(SCRATCH, "ala_grid.zip")
    with open(grid_inp, "w") as f:
        f.write(f"RECEPTOR_FILE {out_mutant}\n")
        f.write(f"GRID_CENTER {center}\n")
        f.write(f"INNERBOX 10, 10, 10\n")
        f.write(f"OUTERBOX 30, 30, 30\n")
    
    run_command([os.path.join(SCHRODINGER, "glide"), grid_inp, "-WAIT", "-LOCAL"], "gridgen.log")
    
    if os.path.exists("gridgen.zip"):
        shutil.move("gridgen.zip", grid_file)
    elif os.path.exists(os.path.join(SCRATCH, "gridgen.zip")):
        shutil.move(os.path.join(SCRATCH, "gridgen.zip"), grid_file)
    
    # 3. Glide SP
    glide_inp = os.path.join(SCRATCH, "glide_sp.inp")
    glide_out = os.path.join(SCRATCH, "glide_sp_out.maegz")
    with open(glide_inp, "w") as f:
        f.write(f"GRIDFILE {grid_file}\n")
        f.write(f"LIGAND_FILE {lig_mae}\n")
        f.write(f"PRECISION SP\n")
    
    run_command([GLIDE, glide_inp, "-WAIT", "-LOCAL"], "glide_sp.log")
    
    if os.path.exists("glide_sp_pv.maegz"):
        shutil.move("glide_sp_pv.maegz", os.path.join(SCRATCH, "glide_sp_pv.maegz"))
    
    # 4. Filter
    passed_ligands = os.path.join(SCRATCH, "pass_filter.maegz")
    passed_list = []
    
    sp_results = os.path.join(SCRATCH, "glide_sp_pv.maegz")
    if not os.path.exists(sp_results):
        # Fallback if name differs
        sp_results = os.path.join(SCRATCH, "glide_sp_lib.maegz")

    if os.path.exists(sp_results):
        for st in structure.StructureReader(sp_results):
            gscore = st.property.get('r_i_glide_gscore', 0)
            if gscore < -5.0:
                print(f"   Pass: {st.title} (GScore: {gscore:.2f})")
                passed_list.append(st)
            else:
                print(f"   Fail: {st.title} (GScore: {gscore:.2f})")
    
    if not passed_list:
        print("!! No ligands passed the Alanine-Scan filter.")
        return None
        
    with structure.StructureWriter(passed_ligands) as writer:
        for st in passed_list:
            writer.append(st)
    return passed_ligands

def main():
    parser = argparse.ArgumentParser(description="CovFlow Ultra: Advanced Covalent Workflow")
    parser.add_argument("--csv", required=True, help="Input CSV")
    parser.add_argument("--pdb", required=True, help="Protein PDB")
    parser.add_argument("--res", required=True, help="Target Residue ID (e.g. A:797)")
    parser.add_argument("--restype", "--resname", required=True, choices=["CYS", "SER", "LYS", "TYR", "THR"], help="Target residue type")
    parser.add_argument("--rxn", default="michael_addition", help="Reaction type")
    parser.add_argument("--dist", default="2.5", help="Distance constraint in Angstroms")
    parser.add_argument("--host", default="localhost", help="Host")
    
    global SCRATCH
    args = parser.parse_args()
    
    # Dynamic SCRATCH directory
    pdb_base = os.path.basename(args.pdb).split('.')[0]
    SCRATCH = f"SCRATCH_{pdb_base}_{args.res.replace(':','_')}"
    
    if not os.path.exists(SCRATCH): os.makedirs(SCRATCH)

    # Phase 1: Sanitize (Remove native ligands, restore covalent bonds)
    # This must happen BEFORE PrepWizard to remove problematic molecules that cause Lewis errors.
    rec_sanitized, center = sanitize_protein(args.pdb, args.res)
    if not rec_sanitized: sys.exit(1)
    
    # Phase 2: Protein Preparation (Standard PrepWiz)
    rec_prepped_basename = "prepwizard_out.maegz"
    cmd = [PREPWIZ, "-fillsidechains", "-minimize_adj_h", rec_sanitized, rec_prepped_basename, "-WAIT", "-LOCAL"]
    success, _ = run_command(cmd, "prepwizard.log")
    
    rec_prepped = os.path.join(SCRATCH, rec_prepped_basename)
    if os.path.exists(rec_prepped_basename):
        shutil.move(rec_prepped_basename, rec_prepped)
        
    if not success or not os.path.exists(rec_prepped):
        print("!! PrepWizard failed to produce output. Check SCRATCH/prepwizard.log")
        sys.exit(1)
    

    
    # Phase 3: Prepare Ligands + Filter
    lig_filtered = prepare_ligands_with_filter(args.csv)
    if not lig_filtered: sys.exit(1)
    
    # Phase 4: Alanine Scan Filter (REMOVED)
    # The Alanine scan is redundant for known covalent inhibitors and was conflicting with Cys797 mutations.
    # lig_production = run_alanine_scan(rec_sanitized, lig_filtered, center, args.res)
    # if not lig_production: sys.exit(1)
    lig_production = lig_filtered

    # Phase 5: Covalent Docking
    # Parse NUCLEOPHILE_ATOM using ASL based on Target Residue Type
    nuc_atom = NUCLEOPHILES.get(args.restype, "SG")
    chain, resnum = args.res.split(':')
    nuc_asl = f"(res.ptype {args.restype}) AND (res.num {resnum}) AND (atom.ptype {nuc_atom})"

    # Reaction mapping for Schrödinger 2021-2
    RXN_MAP = {
        "michael_addition": "Michael Addition",
        "nucleophilic_substitution": "Nucleophilic Substitution",
        "boronic_acid": "Boronic Acid Addition",
        "epoxide": "Epoxide Opening",
        "disulfide": "Disulfide Formation"
    }
    reaction_str = RXN_MAP.get(args.rxn.lower(), args.rxn)

    job_name = "covdock_final"
    inp_file = os.path.join(SCRATCH, f"{job_name}.inp")
    
    # Custom Robust Michael Addition SMARTS
    # Ligand: Michael acceptor beta carbon
    lig_smarts = "[C:1]=[C:2]-[C:3]=[O,S:4]"
    # Receptor: Nucleophilic C/S/O
    rec_smarts = f"[{nuc_atom}:5]"
    
    with open(inp_file, "w") as f:
        f.write(f"RXN_TYPE {reaction_str}\n")
        f.write(f"REC_FILE {os.path.abspath(rec_prepped)}\n")
        f.write(f"LIG_FILE {os.path.abspath(lig_production)}\n")
        f.write(f"ATTACHMENT_RESIDUE {args.res}\n")
        f.write(f"GRID_OPTION GRID_CENTER={center}\n")
        f.write(f"GRID_OPTION INNERBOX=10,10,10\n")
        f.write(f"GRID_OPTION OUTERBOX=30,30,30\n")
        f.write(f"DIST_CONSTRAINT {args.dist}\n")

    
    print(f"\n[PHASE 5] Executing Final Covalent Docking (Standard Leadopt Mode)")
    
    # Pre-clean stale outputs to prevent false successes if docking fails silently
    stale_patterns = [f"{job_name}-out.maegz", f"{job_name}-out.csv", f"{job_name}.csv", f"{job_name}_out.maegz", f"{job_name}_out.csv"]
    for f in stale_patterns:
        if os.path.exists(f):
            os.remove(f)
            
    # Reverting to 'leadopt' (thorough is not supported in this 2021 version)
    cmd = [COVDOCK, os.path.abspath(inp_file), "-mode", "leadopt", "-HOST", args.host, "-WAIT"]
    success, _ = run_command(cmd, "covdock.log")

    if success:
        print("\n[PHASE 6] Organizing Results")
        os.makedirs(PYPROJECT_RAW, exist_ok=True)
        os.makedirs("DATA", exist_ok=True)
        
        # Schrodinger uses hyphens for output files: name-out.maegz
        patterns = [f"{job_name}-out.maegz", f"{job_name}-out.csv", f"{job_name}.csv", f"{job_name}_out.maegz", f"{job_name}_out.csv"]
        
        for f in patterns:
            if os.path.exists(f):
                # 1. Copy to Dropbox
                target_dropbox = os.path.join(PYPROJECT_RAW, f)
                shutil.copy(f, target_dropbox)
                
                # 2. Copy to DATA folder for easy analysis
                if "-out.maegz" in f or "_out.maegz" in f:
                    shutil.copy(f, "DATA/results.maegz")
                elif "-out.csv" in f or "_out.csv" in f:
                    shutil.copy(f, "DATA/results_ranking.csv")
                elif ".rept" in f:
                    shutil.copy(f, "DATA/results.rept")
                else:
                    shutil.copy(f, f"DATA/{f}")
                    
                print(f"   Stored result: {f} -> DATA/ and Dropbox")

        print("\nWorkflow Complete! Check DATA/ and Dropbox for results.")
    else:
        print(f"\n!! Covalent Docking failed. Check {SCRATCH}/covdock.log for details.")

if __name__ == "__main__":
    main()
