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

def run_command(cmd, log_file=None, cwd=None):
    if not cwd:
        cwd = SCRATCH
    if log_file and not os.path.isabs(log_file):
        log_file = os.path.join(cwd, log_file)
    
    print(f"-> Executing: {' '.join(cmd)}")
    
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        with open(log_file, "a") as f:
            f.write(f"\n--- {time.ctime()} ---\n")
            f.write(f"COMMAND: {' '.join(cmd)}\n")

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, cwd=cwd)
    output_lines = []
    for line in process.stdout:
        print(line, end='', flush=True)
        output_lines.append(line)
        if log_file:
            with open(log_file, "a") as f:
                f.write(line)
    process.wait()
    return process.returncode == 0, "".join(output_lines)

def sanitize_protein(pdb_file, target_res_id, skip_min=False):
    """
    1. Removes native ligands.
    2. Restores covalent residues.
    3. Calculates grid center based on native ligand centroid (if present) or nucleophile.
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
    
    # 3. Native Ligand Identification (for Grid Centering)
    # Identify non-standard residues within 10A of the target residue
    native_ligand_atoms = []
    for res in st.residue:
        res_name = res.pdbres.strip()
        if res_name not in COFACTORS and res_name not in STANDARD_RESIDUES:
            # Check if this residue is near the target
            res_atoms = [st.atom[i] for i in res.getAtomIndices()]
            res_coords = np.mean([a.xyz for a in res_atoms], axis=0)
            target_coords = np.mean([a.xyz for a in target_atom_objs], axis=0)
            if np.linalg.norm(res_coords - target_coords) < 15.0:
                native_ligand_atoms.extend(res.getAtomIndices())
                print(f"   Found potential native ligand: {res_name} ({res.resnum})")
                break # Take the first one found near the pocket

    if native_ligand_atoms:
        coords = [st.atom[i].xyz for i in native_ligand_atoms]
        centroid = np.mean(coords, axis=0)
        center = f"{centroid[0]:.3f},{centroid[1]:.3f},{centroid[2]:.3f}"
        print(f"   Grid Center (Native Ligand Centroid): {center}")
    else:
        # Fallback to the Nucleophile (e.g., Cys SG)
        nuc_atom_name = NUCLEOPHILES.get(target_res.pdbres.strip(), "SG")
        nuc_atom_list = [a for a in target_atom_objs if a.pdbname.strip() == nuc_atom_name]
        
        if nuc_atom_list:
            center_coords = nuc_atom_list[0].xyz
            center = f"{center_coords[0]:.3f},{center_coords[1]:.3f},{center_coords[2]:.3f}"
            print(f"   Grid Center (Nucleophile): {center}")
        else:
            # Emergency Fallback to geometric center of residue
            coords = [a.xyz for a in target_atom_objs]
            centroid = np.mean(coords, axis=0)
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
    if not skip_min:
        try:
            print(f"   Executing local minimization of {target_res_id} + 5.0A...")
            # ASL: Target residue and neighbors within 5.0 A
            minimize_asl = f"(res.n {resnum} and chain {chain}) or (within 5.0 (res.n {resnum} and chain {chain}))"
            minimize.minimize_structure(st, iteration=50, asl=minimize_asl)
        except Exception as e:
            print(f"   ⚠️ Minimization warning (skipping): {e}")
    else:
        print(f"   Skipping local minimization as requested (--no_min)")
    
    out_file = os.path.join(SCRATCH, "sanitized_protein.maegz")
    st.write(out_file)
    return out_file, center

def prepare_ligands_with_filter(csv_path, grid_center=None):
    """LigPrep + Warhead Filtering + Reactive Atom Labeling."""
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
    cmd = [LIGPREP, "-ismi", os.path.abspath(smi_file), "-omae", output_basename, "-ph", "7.0", "-s", "1", "-WAIT", "-LOCAL"]
    success, _ = run_command(cmd, "ligprep.log", cwd=SCRATCH)
    
    # Ligprep with -LOCAL will put file in the cwd (SCRATCH)
    full_output_path = os.path.join(SCRATCH, output_basename)
    if os.path.exists(full_output_path):
        if full_output_path != output_mae:
            shutil.move(full_output_path, output_mae)
    
    if success and os.path.exists(output_mae):
        # Post-process MAE to label reactive carbon
        st_list = list(structure.StructureReader(output_mae))
        
        target_center = None
        if grid_center:
            try:
                target_center = np.array([float(x) for x in grid_center.split(',')])
            except:
                pass

        new_st_list = []
        for st in st_list:
            best_match_idx = None
            min_dist = 999.0
            
            for name, pattern_data in detector.patterns.items():
                matches = analyze.evaluate_smarts_canvas(st, pattern_data["smarts"])
                if matches:
                    # If multiple matches exist, pick the one closest to the grid center
                    if target_center is not None:
                        for match in matches:
                            reactive_idx = match[0] # Heuristic: first atom of match is reactive
                            atom_coords = st.atom[reactive_idx].xyz
                            dist = np.linalg.norm(np.array(atom_coords) - target_center)
                            if dist < min_dist:
                                min_dist = dist
                                best_match_idx = reactive_idx
                    else:
                        # Fallback to first match
                        best_match_idx = matches[0][0]
                    
                    if best_match_idx:
                        st.property['i_user_reactive_atom'] = best_match_idx
                        print(f"   Labeled reactive atom {best_match_idx} for {st.title}")
                        break
            new_st_list.append(st)

        with structure.StructureWriter(output_mae) as writer:
            for st in new_st_list:
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

def generate_affinity_grid(rec_prepped_mae, center, target_res_id):
    """
    Generate a Glide affinity grid for binding affinity scoring.
    
    Args:
        rec_prepped_mae: Path to prepared receptor MAE file
        center: Grid center coordinates (x,y,z)
        target_res_id: Target residue ID (e.g., "A:797")
    
    Returns:
        Path to generated grid file (.zip) or None if failed
    """
    print(f"\n[PHASE 4B] Generating Affinity Grid for Glide Scoring")
    
    # Parse center coordinates
    try:
        cx, cy, cz = map(float, center.split(','))
    except:
        print(f"!! Error: Could not parse center coordinates: {center}")
        return None
    
    grid_name = "covdock_affinity"
    grid_in_file = os.path.join(SCRATCH, f"{grid_name}.in")
    grid_zip_file = os.path.join(SCRATCH, f"{grid_name}_glidegrid.zip")
    
    # Create Glide grid generation input file
    # This configures Glide to generate a grid with affinity potential for binding affinity scoring
    with open(grid_in_file, "w") as f:
        f.write("GRID_CENTER    {:.3f}, {:.3f}, {:.3f}\n".format(cx, cy, cz))
        f.write("INNERBOX       10, 10, 10\n")
        f.write("OUTERBOX       30, 30, 30\n")
        f.write("RECEPTOR       {}\n".format(os.path.abspath(rec_prepped_mae)))
        f.write("GRIDFILE       {}\n".format(os.path.abspath(grid_zip_file)))
        f.write("AFFINITY_GRID  True\n")  # Critical: enable affinity potential generation
    
    # Run Glide to generate grid (use basename since glide runs in SCRATCH directory)
    cmd = [GLIDE, f"{grid_name}.in", "-WAIT"]
    success, _ = run_command(cmd, "grid_generation.log", cwd=SCRATCH)
    
    if success and os.path.exists(grid_zip_file):
        print(f"   ✓ Affinity grid generated: {grid_zip_file}")
        return grid_zip_file
    else:
        print(f"!! Glide grid generation failed. Check {SCRATCH}/grid_generation.log")
        return None

def main():
    parser = argparse.ArgumentParser(description="CovFlow Ultra: Advanced Covalent Workflow")
    parser.add_argument("--csv", required=True, help="Input CSV")
    parser.add_argument("--pdb", required=True, help="Protein PDB")
    parser.add_argument("--res", required=True, help="Target Residue ID (e.g. A:797)")
    parser.add_argument("--restype", "--resname", required=True, choices=["CYS", "SER", "LYS", "TYR", "THR"], help="Target residue type")
    parser.add_argument("--rxn", default="michael_addition", help="Reaction type")
    parser.add_argument("--dist", default="2.5", help="Distance constraint in Angstroms")
    parser.add_argument("--host", default="localhost", help="Host")
    parser.add_argument("--no_min", action="store_true", help="Skip local minimization of target residue (use for native ligands)")
    
    global SCRATCH
    args = parser.parse_args()
    
    # Dynamic SCRATCH directory
    pdb_base = os.path.basename(args.pdb).split('.')[0]
    SCRATCH = f"SCRATCH_{pdb_base}_{args.res.replace(':','_')}"
    
    if not os.path.exists(SCRATCH): os.makedirs(SCRATCH)
    
    # Structure-Specific Awareness (e.g., 7K1I Allosteric Awareness)
    if "7K1I" in args.pdb.upper():
        print("\n[INFO] Structure 7K1I detected: Aware of Inactive/Allosteric-inhibitor-bound state.")
        print("       Grid centering will prioritize native ligand centroid to cover ATP-binding pocket.")

    # Phase 1: Sanitize (Remove native ligands, restore covalent bonds)
    # This must happen BEFORE PrepWizard to remove problematic molecules that cause Lewis errors.
    rec_sanitized, center = sanitize_protein(args.pdb, args.res, skip_min=args.no_min)
    if not rec_sanitized: sys.exit(1)
    
    # Phase 2: Protein Preparation (Standard PrepWiz)
    rec_prepped_basename = "prepwizard_out.maegz"
    # Note: protassign and epik are run by default. Removed -propassign and -epik as they are not valid flags in this version.
    cmd = [PREPWIZ, "-fillsidechains", "-samplewater", "-minimize_adj_h", os.path.abspath(rec_sanitized), rec_prepped_basename, "-WAIT", "-LOCAL"]
    success, _ = run_command(cmd, "prepwizard.log", cwd=SCRATCH)
    
    rec_prepped = os.path.join(SCRATCH, rec_prepped_basename)
        
    if not success or not os.path.exists(rec_prepped):
        print("!! PrepWizard failed to produce output. Check SCRATCH/prepwizard.log")
        sys.exit(1)
    
    # Phase 3: Prepare Ligands + Filter
    lig_filtered = prepare_ligands_with_filter(args.csv, grid_center=center)
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
    
    # Phase 4B: Generate Affinity Grid for Glide Scoring
    affinity_grid = generate_affinity_grid(rec_prepped, center, args.res)
    if not affinity_grid:
        print("!! Warning: Affinity grid generation failed. Proceeding with geometry-only scoring.")
    
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
        # Add affinity grid if successfully generated
        if affinity_grid and os.path.exists(affinity_grid):
            f.write(f"AFFINITY_GRID  {os.path.abspath(affinity_grid)}\n")
        f.write(f"GRID_OPTION GRID_CENTER={center}\n")
        f.write(f"GRID_OPTION INNERBOX=10,10,10\n")
        f.write(f"GRID_OPTION OUTERBOX=30,30,30\n")
        f.write(f"DIST_CONSTRAINT 4.0\n")
        f.write(f"MAX_INIT_POSES 1000\n")

    
    print(f"\n[PHASE 5] Executing Final Covalent Docking (Standard Leadopt Mode + Affinity Scoring)")

    
    # Pre-clean stale outputs to prevent false successes if docking fails silently
    stale_patterns = [f"{job_name}-out.maegz", f"{job_name}-out.csv", f"{job_name}.csv", f"{job_name}_out.maegz", f"{job_name}_out.csv"]
    for f in stale_patterns:
        if os.path.exists(f):
            os.remove(f)
            
    # Reverting to 'leadopt' (thorough is not supported in this 2021 version)
    cmd = [COVDOCK, os.path.abspath(inp_file), "-mode", "leadopt", "-HOST", args.host, "-WAIT"]
    success, _ = run_command(cmd, "covdock.log", cwd=SCRATCH)

    if success:
        print("\n[PHASE 6] Organizing Results")
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
                    shutil.copy(f_path, "DATA/results.maegz")
                elif "-out.csv" in f_name or "_out.csv" in f_name:
                    shutil.copy(f_path, "DATA/results_ranking.csv")
                elif ".rept" in f_name:
                    shutil.copy(f_path, "DATA/results.rept")
                else:
                    shutil.copy(f_path, f"DATA/{f_name}")
                    
                print(f"   Stored result: {f_name} -> DATA/ and {PYPROJECT_RAW}/")

        # [PHASE 7] Red Flag Analysis
        print("\n[PHASE 7] Running Post-Docking Red Flag Analysis")
        results_file = "DATA/results.maegz"
        if os.path.exists(results_file):
            analysis_script = os.path.join(os.path.dirname(__file__), "..", "ANALYSIS", "red_flag_filter.py")
            cmd = [RUN, "python3", os.path.abspath(analysis_script), results_file]
            run_command(cmd, "red_flag_analysis.log", cwd=SCRATCH)

        # [PHASE 8] MedChem & Ligand Efficiency Analysis
        print("\n[PHASE 8] Running MedChem & Ligand Efficiency Analysis")
        if os.path.exists(results_file):
            le_script = os.path.join(os.path.dirname(__file__), "..", "ANALYSIS", "ligand_efficiency.py")
            output_csv = "DATA/medchem_analysis.csv"
            cmd = [RUN, "python3", os.path.abspath(le_script), results_file, "--output", output_csv]
            run_command(cmd, "medchem_analysis.log", cwd=SCRATCH)
            
            # Copy CSV to results folder as well
            if os.path.exists(os.path.join(SCRATCH, output_csv)):
                shutil.copy(os.path.join(SCRATCH, output_csv), os.path.join(PYPROJECT_RAW, "medchem_analysis.csv"))
            elif os.path.exists(output_csv): # If it was written to current dir
                shutil.copy(output_csv, os.path.join(PYPROJECT_RAW, "medchem_analysis.csv"))

        print("\nWorkflow Complete! Check DATA/ and Results for results.")
    else:
        print(f"\n!! Covalent Docking failed. Check {SCRATCH}/covdock.log for details.")

if __name__ == "__main__":
    main()
