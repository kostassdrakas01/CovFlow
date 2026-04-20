import os
import argparse
import sys
from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.structutils.interactions import hbond

def analyze_interactions(mae_file, target_res="A:793"):
    """
    Universal H-bond analysis.
    Works on any .maegz file (auto-detects protein/ligands).
    """
    if not os.path.exists(mae_file):
        print(f"!! Error: File {mae_file} not found.")
        return

    print(f"\n========================================")
    print(f" UNIVERSAL INTERACTION ANALYSIS")
    print(f" Target Residue: {target_res}")
    print(f"========================================")
    
    try:
        st_reader = structure.StructureReader(mae_file)
        
        for st in st_reader:
            lig_name = st.title or "Unknown_Ligand"
            
            # Identify ligand and protein atoms using ASL
            lig_indices = analyze.evaluate_asl(st, "ligand")
            prot_indices = analyze.evaluate_asl(st, "protein")
            
            if not lig_indices:
                # If everything has the same resname, skip the receptor-only entry
                if len(st.atom) > 500:
                    print(f"--- Detected Receptor Registry: {st.title} ---")
                continue
            
            # Calculate H-bonds between the ligand and protein within the structure
            hbonds = hbond.get_hydrogen_bonds(st, atoms1=lig_indices, atoms2=prot_indices)
            
            # Find target residue atoms
            chain, resnum = target_res.split(':')
            relevant = []
            for hb in hbonds:
                a1, a2 = hb[0], hb[1]
                # Check if either atom is our target residue
                if (a1.chain == chain and a1.resnum == int(resnum)) or \
                   (a2.chain == chain and a2.resnum == int(resnum)):
                    relevant.append(hb)

            if relevant:
                print(f" {lig_name:25} | H-Bonds with {target_res}: {len(relevant)}")
            else:
                print(f" {lig_name:25} | No H-bond with {target_res}")
                
    except Exception as e:
        print(f"!! Error processing file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze H-bond interactions from CovDock output.")
    parser.add_argument("file", help="Path to the .maegz results file")
    parser.add_argument("--res", default="A:793", help="Target residue (e.g., A:793)")
    args = parser.parse_args()
    
    analyze_interactions(args.file, args.res)
