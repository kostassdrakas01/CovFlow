# CovFlow Reference Manual

Welcome to your organized CovFlow workspace. 

## Folder Structure
- **`/home/konstantinos/Documents/CovFlow/`**: Main workspace.
- **`BIN/`**: Core automated scripts.
- **`ANALYSIS/`**: Universal tools for docking result analysis.
- **`1XKK.pdb`**: Your example protein.
- **`ligands.csv`**: Your example ligand list.

---

## 1. Running New Projects
You can run this workflow on **any** protein and **any** covalent docking reaction.

### Standard Docking (Fast Ranking)
```bash
covflow --csv your_data.csv --pdb your_protein.pdb
```

### Kinetics Docking (K_i and k_inact)
```bash
covflow_kin --csv your_data.csv --pdb your_protein.pdb
```

**Interactive Prompts:**
- **Reaction Type**: Type any Schrodinger reaction (e.g., `michael_addition`, `nucleophilic_substitution`).
- **Target Residue**: Type the chain and number (e.g., `A:797`, `B:102`).

---

## 2. Docking Site & Complex Proteins
CovFlow Ultra includes an automated **Sanitization Phase** to handle complex PDB files.

### HET Groups & Native Ligands
- **Automatic Removal**: The pipeline automatically detects and removes all non-standard HET groups and cofactors from the PDB to prevent clashes.
- **Site Healing**: If your target residue already has a covalent bond to a ligand, CovFlow will break the bond and "heal" the pocket via local minimization before docking.
- **Manual Site Selection**: The docking site is strictly defined by the `--res` argument. It does *not* automatically dock where the native ligand was unless you specify that residue.

### Multiple Chains
- **Full Support**: You can target any residue on any chain (e.g., `--res B:123`).
- **Global Context**: All chains are preserved in the receptor structure to maintain the correct biological environment, but the ligand attachment is restricted to your chosen site.

---

## 3. Analyzing Results
Use these tools on **any** `.maegz` or `.pdb` file (from this workflow or manual Maestro runs).

### Interaction Analysis (H-Bonds)
Check if ligands bind to the hinge (Met793) or any other residue.
```bash
covflow run python3 ANALYSIS/interaction_analysis.py results.maegz --res A:793
```

### Strain Analysis (Torsion)
Check the geometric strain of the acrylamide warhead.
```bash
covflow run python3 ANALYSIS/strain_analysis.py results.maegz
```

### Master Analysis (Custom ASL)
Calculate the distance between **any** selection (using Maestro ASL syntax).
```bash
covflow run python3 ANALYSIS/master_analysis.py results.maegz --sel1 "ligand" --sel2 "res.pt A:793"
```

---

## ASL Selection Cheat Sheet
- **Chain**: `chain A`
- **Residue**: `res.num 797`
- **Specific Atom**: `res.pt A:797 and atom.name SG`
- **Ligand**: `ligand`
- **Water**: `res.nam HOH`
- **Complexes**: `protein and not water`


 ./covflow_run.sh --csv DATA/ligands.csv --pdb DATA/7UKV.pdb --res A:797 --restype CYS
