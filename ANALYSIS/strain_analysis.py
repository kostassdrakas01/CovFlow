import argparse
import sys
import os
from schrodinger import structure
from schrodinger.structutils import analyze

def analyze_warhead_geometry(mae_file):
    """
    Universal Strain Analysis.
    Finds the warhead (formerly an acrylamide) and measures its torsion.
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    # SMARTS for the saturated acrylamide (post-reaction): C-C-C=O
    WARHEAD_SMARTS = "[C;H2,H1]-[C;H1,H0]-C=O" 

    print(f"\n========================================")
    print(f" UNIVERSAL WARHEAD STRAIN ANALYSIS")
    print(f"========================================")

    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            lig_name = st.title or "Unknown_Ligand"
            
            # Run SMARTS on the structure (works on both lone ligands and complexes)
            matches = analyze.evaluate_smarts_canvas(st, WARHEAD_SMARTS)
            
            if matches:
                for match in matches:
                    torsion = st.measure(*match)
                    print(f" {lig_name:25} | Warhead Torsion: {torsion:6.2f}°")
            else:
                print(f"❓ {lig_name:25} | No acrylamide-like group found.")

    except Exception as e:
        print(f"!! Error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze geometric strain of docked ligands.")
    parser.add_argument("file", help="Path to the .maegz results file")
    args = parser.parse_args()
    
    analyze_warhead_geometry(args.file)
