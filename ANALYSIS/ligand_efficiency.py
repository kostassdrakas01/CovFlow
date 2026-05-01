#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
from schrodinger import structure
from schrodinger.structutils import analyze

def analyze_ligand_efficiency(mae_file, mw_cutoff=560.0, rb_penalty=0.1, output_csv="medchem_analysis.csv"):
    """
    Calculates Ligand Efficiency and Adjusted Scores for docked poses.
    
    Args:
        mae_file: Path to the .maegz results file.
        mw_cutoff: Molecular weight threshold for filtering.
        rb_penalty: Penalty added to the score per rotatable bond.
        output_csv: Name of the output CSV file.
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    print(f"\n" + "="*80)
    print(f" MEDICINAL CHEMISTRY ANALYSIS & LIGAND EFFICIENCY")
    print(f" Source: {os.path.basename(mae_file)}")
    print(f" Filters: MW < {mw_cutoff} Da | RB Penalty: {rb_penalty} per bond")
    print(f"="*80)

    results = []
    total_poses = 0
    passed_filter = 0

    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            total_poses += 1
            lig_name = st.title or f"Pose_{total_poses}"
            
            # Extract only the ligand atoms for accurate MW and HA calculation
            # In covalent docking results, the structure contains both protein and ligand.
            ligands = analyze.find_ligands(st)
            
            if ligands:
                # Sort ligands by weight and pick the largest one
                # (Assuming the largest non-protein/non-water entity is our ligand)
                ligands.sort(key=lambda x: x.st.total_weight, reverse=True)
                lig_st = ligands[0].st
            else:
                # Fallback to ASL if find_ligands fails
                lig_indices = analyze.evaluate_asl(st, "ligand")
                if not lig_indices:
                    lig_indices = analyze.evaluate_asl(st, "not protein and not water and not (res.name HOH,MG,ZN,CL,SO4)")
                
                if lig_indices:
                    lig_st = st.extract(lig_indices)
                else:
                    lig_st = st # Fallback to full structure
            
            # 1. Get Docking Score
            # Properties are case-sensitive and vary by Schrodinger version/application
            docking_score = st.property.get('r_i_docking_score', 
                            st.property.get('r_i_glide_docking_score', 
                            st.property.get('r_i_glide_gscore', 
                            st.property.get('r_lp_Energy', 0.0))))
            
            # 2. Calculate Molecular weight of the LIGAND
            mw = lig_st.total_weight
            
            # 3. Count Heavy Atoms of the LIGAND
            num_heavy_atoms = len(analyze.evaluate_asl(lig_st, "not (atom.ele H)"))
            
            # 4. Count Rotatable Bonds of the LIGAND
            # Check properties first (most accurate if present)
            num_rotatable_bonds = st.property.get('i_i_glide_rotatable_bonds',
                                  st.property.get('i_m_rotatable_bonds',
                                  st.property.get('i_lp_num_rotatable_bonds', None)))
            
            if num_rotatable_bonds is None:
                # Dynamic calculation fallback
                try:
                    num_rotatable_bonds = analyze.count_rotatable_bonds(lig_st)
                except:
                    num_rotatable_bonds = 0
            
            # 5. Calculate Metrics
            le = 0.0
            if num_heavy_atoms > 0:
                le = docking_score / num_heavy_atoms
            
            # Adjusted Score (Penalize flexibility)
            # Since Docking Score is negative (better is more negative), 
            # adding a positive penalty makes it worse.
            adjusted_score = docking_score + (num_rotatable_bonds * rb_penalty)
            
            # 6. Apply MW Filter
            status = "PASS"
            if mw > mw_cutoff:
                status = "FAIL (MW)"
            else:
                passed_filter += 1

            results.append({
                "Ligand": lig_name,
                "Docking_Score": docking_score,
                "MW": mw,
                "Heavy_Atoms": num_heavy_atoms,
                "Rotatable_Bonds": num_rotatable_bonds,
                "LE": le,
                "Adjusted_Score": adjusted_score,
                "Status": status
            })

    except Exception as e:
        print(f"!! Error processing structure: {e}")
        return

    # Create DataFrame and Save
    df = pd.DataFrame(results)
    
    # Sort by Adjusted Score (Lower is better)
    df = df.sort_values(by="Adjusted_Score")
    
    df.to_csv(output_csv, index=False)
    
    # Display Top 10
    print(f"\n {'Ligand Name':25} | {'Score':7} | {'MW':6} | {'HA':3} | {'RB':3} | {'LE':6} | {'AdjScore':8} | {'Status'}")
    print(f" " + "-"*90)
    
    for _, row in df.head(15).iterrows():
        print(f" {row['Ligand']:25} | {row['Docking_Score']:7.2f} | {row['MW']:6.1f} | {int(row['Heavy_Atoms']):3} | {int(row['Rotatable_Bonds']):3} | {row['LE']:6.3f} | {row['Adjusted_Score']:8.2f} | {row['Status']}")

    print(f" " + "-"*90)
    print(f" Summary: {passed_filter}/{total_poses} poses passed MW filter (<{mw_cutoff}).")
    print(f" Results saved to: {output_csv}")
    print(f"="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ligand Efficiency and MedChem filtering for docked poses.")
    parser.add_argument("file", help="Path to the .maegz results file")
    parser.add_argument("--mw", type=float, default=560.0, help="Molecular weight cutoff (default: 560.0)")
    parser.add_argument("--penalty", type=float, default=0.1, help="Penalty per rotatable bond (default: 0.1)")
    parser.add_argument("--output", default="medchem_analysis.csv", help="Output CSV filename")
    
    args = parser.parse_args()
    
    analyze_ligand_efficiency(args.file, args.mw, args.penalty, args.output)
