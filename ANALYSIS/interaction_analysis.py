import os
import argparse
import sys
from schrodinger import structure
from schrodinger.structutils import analyze, measure
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
            
            # 1. H-bond Analysis
            hbonds = hbond.get_hydrogen_bonds(st, atoms1=lig_indices, atoms2=prot_indices)
            
            # Find target residue atoms
            chain, resnum = target_res.split(':')
            relevant_hbonds = []
            for hb in hbonds:
                a1, a2 = hb[0], hb[1]
                if (a1.chain == chain and a1.resnum == int(resnum)) or \
                   (a2.chain == chain and a2.resnum == int(resnum)):
                    relevant_hbonds.append(hb)

            # 2. Distance-based Analysis (Relaxed Criteria)
            # Find ligand Nitrogens and Target O/N atoms
            lig_n_asl = "ligand and atom.el N"
            target_on_asl = f"(res.n {resnum} and chain {chain}) and (atom.el O,N)"
            
            lig_n_idx = analyze.evaluate_asl(st, lig_n_asl)
            target_on_idx = analyze.evaluate_asl(st, target_on_asl)
            
            min_dist = 99.0
            if lig_n_idx and target_on_idx:
                # Efficient nested loop to find minimum distance
                for i in lig_n_idx:
                    for j in target_on_idx:
                        dist = st.measure(i, j)
                        if dist < min_dist:
                            min_dist = dist

            # Reporting logic: H-bond OR Distance < 3.5A
            hbond_count = len(relevant_hbonds)
            if hbond_count > 0 or min_dist < 3.5:
                status = "PASS"
                detail = f"H-bonds: {hbond_count}, MinDist: {min_dist:.2f}A"
                print(f" {lig_name:25} | {status} | {detail}")
            else:
                print(f" {lig_name:25} | FAIL | No H-bond and MinDist: {min_dist:.2f}A")
                
    except Exception as e:
        print(f"!! Error processing file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze H-bond interactions from CovDock output.")
    parser.add_argument("file", help="Path to the .maegz results file")
    parser.add_argument("--res", default="A:793", help="Target residue (e.g., A:793)")
    args = parser.parse_args()
    
    analyze_interactions(args.file, args.res)
