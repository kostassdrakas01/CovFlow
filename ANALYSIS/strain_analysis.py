import argparse
import sys
import os
import csv
from schrodinger import structure
from schrodinger.structutils import analyze

def analyze_warhead_geometry(mae_file):
    """
    Universal Strain Analysis.
    Finds the warhead (formerly an acrylamide) and measures its torsion.
    Saves results to CSV and determines PASS/FAIL.
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    # SMARTS for the saturated acrylamide (post-reaction): C-C-C=O
    WARHEAD_SMARTS = "[C;H2,H1]-[C;H1,H0]-C=O" 
    
    # Output CSV filename
    base_name = os.path.splitext(mae_file)[0]
    csv_file = f"{base_name}_torsion.csv"

    print(f"\n========================================")
    print(f" UNIVERSAL WARHEAD STRAIN ANALYSIS")
    print(f" Output: {csv_file}")
    print(f"========================================")

    results = []

    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            lig_name = st.title or "Unknown_Ligand"
            
            # Run SMARTS on the structure
            matches = analyze.evaluate_smarts_canvas(st, "[C]-[C]-C=O")
            
            seen_warheads = set()
            # Expanded list including PrepWizard/Epik variants
            STANDARD_AMINO_ACIDS = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", 
                "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                "CYX", "CYM", "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "ARN"
            }
            
            if matches:
                for match in matches:
                    match_atoms = [st.atom[idx] for idx in match]
                    res_name = match_atoms[0].pdbres.strip().upper()
                    
                    # Debug: uncomment to see all residue names being matched
                    print(f"DEBUG: Found match in residue {res_name}")

                    # Filter: Ensure the match is NOT in a standard/modified protein residue
                    if res_name in STANDARD_AMINO_ACIDS:
                        continue
                        
                    # Sort indices to avoid reporting the same warhead twice
                    warhead_id = tuple(sorted(match))
                    if warhead_id in seen_warheads:
                        continue
                    seen_warheads.add(warhead_id)
                        
                    torsion = st.measure(*match)
                    
                    # Criteria: PASS if staggered, FAIL if eclipsed
                    abs_torsion = abs(torsion)
                    status = "FAIL" if 75 < abs_torsion < 135 else "PASS"
                    
                    print(f" {lig_name:25} | {status} | Warhead Torsion: {torsion:6.2f}°")
                    results.append([lig_name, status, f"{torsion:.2f}"])
            
            if not seen_warheads:
                results.append([lig_name, "NO_WARHEAD", "N/A"])

        # Write to CSV
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Ligand", "Strain_Status", "Torsion_Angle"])
            writer.writerows(results)
            
        print(f"\n>> Analysis Complete. {len(results)} entries saved to: {csv_file}")

    except Exception as e:
        print(f"!! Error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze geometric strain of docked ligands.")
    parser.add_argument("file", help="Path to the .maegz results file")
    args = parser.parse_args()
    
    analyze_warhead_geometry(args.file)
