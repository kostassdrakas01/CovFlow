#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
from schrodinger import structure
from schrodinger.structutils import analyze

def analyze_red_flags(mae_file, docking_score_threshold=-8.0, prime_energy_threshold=100.0):
    """
    Analyzes docked poses for red flags:
    1. Convergence Failure
    2. Bond Length Outliers (C-S covalent bond)
    3. Energy Mismatch (Good Docking Score vs High Prime Energy)
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    print(f"\n" + "="*60)
    print(f" COVALENT DOCKING RED FLAG ANALYSIS")
    print(f" Source: {os.path.basename(mae_file)}")
    print(f"="*60)
    print(f" {'Ligand Name':25} | {'Score':7} | {'PrimeE':8} | {'Bond':5} | {'Status'}")
    print(f" " + "-"*75)

    flags_found = 0
    total_poses = 0

    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            total_poses += 1
            lig_name = st.title or f"Pose_{total_poses}"
            
            # 1. Check Convergence
            # Properties: b_lp_minimization_converged, b_ff_Minimization_Converged, etc.
            converged = st.property.get('b_lp_minimization_converged', 
                        st.property.get('b_ff_Minimization_Converged', 
                        st.property.get('b_prime_minimization_converged', True)))
            
            # 2. Check Docking Score and Prime Energy
            # Note: Properties are case-sensitive
            docking_score = st.property.get('r_i_docking_score', 
                            st.property.get('r_i_glide_docking_score', 
                            st.property.get('r_i_glide_gscore', 0.0)))
            
            # lp_energy, r_lp_Energy, or prime_energy
            prime_energy = st.property.get('r_lp_Energy', 
                           st.property.get('r_lp_energy', 
                           st.property.get('r_psp_Prime_Energy', 
                           st.property.get('r_prime_energy', 0.0))))
            
            # 3. Bond Length Analysis
            # Find the covalent bond. In CovDock, it's usually between the receptor and ligand.
            
            bond_length = 0.0
            covalent_bond_found = False
            
            # List of standard residues to distinguish protein from ligand
            STANDARD_RES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP", "CYX"]
            
            # General Heuristic: Find any bond between a standard protein residue and a non-standard residue (ligand)
            for atom in st.atom:
                res_name = atom.pdbres.strip()
                if res_name in STANDARD_RES:
                    for bond in atom.bond:
                        neighbor = bond.atom2
                        n_resname = neighbor.pdbres.strip()
                        # If neighbor is NOT part of a standard residue or common cofactor, it's the ligand
                        if n_resname not in STANDARD_RES and n_resname not in ["HOH", "WAT", "MG", "ZN", "CL", "SO4"]:
                            bond_length = bond.length
                            covalent_bond_found = True
                            break
                if covalent_bond_found: break

            # Apply Logic
            red_flags = []
            
            if not converged:
                red_flags.append("CONV_FAIL")
            
            if covalent_bond_found:
                if bond_length > 2.2 or bond_length < 1.6:
                    red_flags.append("BAD_BOND")
            
            if docking_score < docking_score_threshold and prime_energy > prime_energy_threshold:
                red_flags.append("ENERGY_MISMATCH")

            # Formatting
            status = "PASS"
            if red_flags:
                status = "🚩 FLAG: " + ",".join(red_flags)
                flags_found += 1

            bond_str = f"{bond_length:5.2f}" if covalent_bond_found else "N/A"
            print(f" {lig_name:25} | {docking_score:7.2f} | {prime_energy:8.1f} | {bond_str:5} | {status}")

    except Exception as e:
        print(f"!! Error: {e}")

    print(f" " + "-"*75)
    print(f" Summary: Found {flags_found} flagged poses out of {total_poses} total.")
    print(f"="*60 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze docked poses for red flags.")
    parser.add_argument("file", help="Path to the .maegz results file")
    parser.add_argument("--score", type=float, default=-8.0, help="Docking score threshold for energy mismatch")
    parser.add_argument("--energy", type=float, default=100.0, help="Prime energy threshold for energy mismatch")
    args = parser.parse_args()
    
    analyze_red_flags(args.file, args.score, args.energy)
