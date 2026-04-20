import argparse
import os
import sys
from schrodinger import structure
from schrodinger.structutils import analyze, interactions

def run_custom_analysis(mae_file, sel1_asl, sel2_asl):
    """
    Highly flexible analysis using Schrodinger ASL (Atom Selection Language).
    Calculates the minimum distance between two arbitrary selections.
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    print(f"\n========================================")
    print(f" MASTER ASL ANALYSIS")
    print(f" Sel 1: {sel1_asl}")
    print(f" Sel 2: {sel2_asl}")
    print(f"========================================")

    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            # Evaluate ASL
            try:
                idx1 = analyze.evaluate_asl(st, sel1_asl)
                idx2 = analyze.evaluate_asl(st, sel2_asl)
            except Exception as e:
                print(f"!! ASL Error in {st.title}: {e}")
                continue

            if not idx1:
                print(f"⚠️  {st.title:20} | Selection 1 found 0 atoms.")
                continue
            if not idx2:
                print(f"⚠️  {st.title:20} | Selection 2 found 0 atoms.")
                continue

            # Calculate minimum distance between the two selections
            # get_short_dist returns (index1, index2, distance)
            res = analyze.get_short_dist(st, idx1, st, idx2)
            if res:
                dist = res[2]
                print(f" {st.title:20} | Min Dist: {dist:6.2f} Å")
            else:
                print(f" {st.title:20} | Could not calculate distance.")

    except Exception as e:
        print(f"!! Error processing file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master Analysis: Distance between any two ASL selections.")
    parser.add_argument("file", help="Path to .maegz, .mae, or .pdb file")
    parser.add_argument("--sel1", required=True, help="ASL for Selection 1 (e.g., 'ligand')")
    parser.add_argument("--sel2", required=True, help="ASL for Selection 2 (e.g., 'res.pt A:793')")
    args = parser.parse_args()
    
    run_custom_analysis(args.file, args.sel1, args.sel2)
