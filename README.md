# CovFlow Reference Manual

Welcome to your organized CovFlow workspace. 

## Folder Structure
- **`/home/konstantinos/Documents/CovFlow/`**: Main workspace.
- **`BIN/`**: Core automated scripts.
- **`ANALYSIS/`**: Universal tools for docking result analysis.
- **`1XKK.pdb`**: Your example protein.
- **`ligands.csv`**: Your example ligand list.

---

## 🛠️ Setup & Portability
If you are moving this workflow to a new computer, you only need to update **two paths** in the `covflow_run.sh` file:

1.  **`SCHRODINGER_RUN`**: The path to your Schrödinger `run` utility (e.g., `/opt/schrodinger/2021-2/run`).
2.  **`SCRIPT_PATH`**: The absolute path to where you saved `BIN/covflow.py`.

```bash
# Example edit in covflow_run.sh:
SCHRODINGER_RUN="/your/path/to/schrodinger/run"
SCRIPT_PATH="/your/path/to/CovFlow/BIN/covflow.py"
```

---

## 1. Running New Projects
You can run this workflow on **any** protein and **any** covalent docking reaction. No external virtual environment (venv) is required, as the pipeline runs directly via the Schrödinger environment.

> [!NOTE]
> Always run the scripts using the provided wrappers or via `$SCHRODINGER/run python3` to ensure all chemical libraries are correctly loaded.

### Standard Docking (Fast Ranking)
```bash
./covflow_run.sh --csv DATA/ligands.csv --pdb DATA/1XKK.pdb --res A:797 --soften 0.85
```

**Common Options:**
- `--csv`: Path to your ligand list (SMILES and names).
- `--pdb`: Path to the protein structure.
- `--res`: Target residue (e.g., `A:797`).
- `--soften [0.85]`: (Optional) Scales receptor VdW radii to allow water displacement.
- `--no_min`: (Optional) Skips side-chain minimization to preserve native pocket.
- `--center_res [790,793,797]`: (Optional) Comma-separated residue numbers for grid centering. Use this if your protein has different numbering or you want to center on different residues.


**Interactive Prompts:**
- **Reaction Type**: Type any Schrodinger reaction (e.g., `michael_addition`, `nucleophilic_substitution`).
- **Target Residue**: Type the chain and number (e.g., `A:797`, `B:102`).

---

## 2. Docking Site & Complex Proteins
CovFlow Ultra includes an automated **Sanitization Phase** to handle complex PDB files.

### HET Groups & Native Ligands
- **Automatic Removal**: The pipeline automatically detects and removes all non-standard HET groups and cofactors from the PDB to prevent clashes.
- **Site Healing**: If your target residue already has a covalent bond to a ligand, CovFlow will break the bond and "heal" the pocket via local minimization before docking.
- **Hydrated Docking**: Crystal waters within 5.0 Å of the site are **preserved** to maintain critical H-bond networks.
- **Softening (Water Displacement)**: Use `--soften 0.85` to allow large ligands to displace pocket waters by scaling receptor VdW radii.

### Multiple Chains
- **Full Support**: You can target any residue on any chain (e.g., `--res B:123`).
- **Global Context**: All chains are preserved in the receptor structure to maintain the correct biological environment, but the ligand attachment is restricted to your chosen site.

---

### Red Flag Analysis
Automatically check for bad bond geometry, energy strain, and convergence failures.
```bash
./cov_analysis.sh ANALYSIS/red_flag_filter.py DATA/results.maegz
```

### Ligand Efficiency & MedChem Filtering (New)
Calculate Ligand Efficiency (LE) and apply medicinal chemistry filters (MW, flexibility).
```bash
./cov_analysis.sh ANALYSIS/ligand_efficiency.py DATA/results.maegz --mw 560 --penalty 0.1
```

### Interaction Analysis (H-Bonds)
Check if ligands bind to the hinge (Met793) or any other residue.
```bash
./cov_analysis.sh ANALYSIS/interaction_analysis.py DATA/results_7ukw.maegz --res A:793
```

### Strain Analysis (Torsion)
Check the geometric strain of the acrylamide warhead.
```bash
./cov_analysis.sh ANALYSIS/strain_analysis.py results.maegz
```

### Master Analysis (Custom ASL)
Calculate the distance between **any** selection (using Maestro ASL syntax).
```bash
covflow run python3 ANALYSIS/master_analysis.py results.maegz --sel1 "ligand" --sel2 "res.pt A:793"
```

---

## 3. Theoretical Background

Raw **Docking Scores** can often be misleading, especially when screening large libraries. Larger molecules naturally tend to form more interactions and thus achieve "better" (more negative) scores, regardless of their actual binding quality relative to their size.

### The Problem: Size Bias
A very large molecule (e.g., MW > 600) might have a score of -10.0, while a smaller, highly optimized lead (e.g., MW 350) might score -7.0. In raw ranking, the large molecule wins, but in the lab, the smaller one is often a better starting point for drug development.

### The Solution: Ligand Efficiency (LE)
To "ground" your results in reality, CovFlow calculates **Ligand Efficiency (LE)**:
>[!IMPORTANT]
> **Ligand Efficiency (LE) = Docking Score / Number of Heavy Atoms**

By dividing the score by the number of atoms, we identify molecules that provide the "most bang for their buck." A molecule like **Osimertinib** is often more efficient than larger inhibitors, even if the larger ones have a higher absolute score.

### Refining the Ranking
To further prioritize drug-like leads, CovFlow applies:
1.  **MW Filtering ($< 560 \text{ Da}$)**: Ensures the molecule stays within the range of standard oral drugs.
2.  **Flexibility Penalty**: High flexibility (many rotatable bonds) is entropically unfavorable. CovFlow adds a penalty to the score for each rotatable bond to prioritize rigid, pre-organized binders.

---

## 4. ASL Selection Cheat Sheet
- **Chain**: `chain A`
- **Residue**: `res.num 797`
- **Specific Atom**: `res.pt A:797 and atom.name SG`
- **Ligand**: `ligand`
- **Water**: `res.nam HOH`
- **Complexes**: `protein and not water`

---

## 4. Workspace Maintenance
The pipeline generates intermediate logs and temporary structures in `SCRATCH` directories. To keep your workspace clean and organized, use the provided cleanup utility:
```bash
./clean.sh
```
This will safely remove all temporary directories while keeping your core scripts and results intact.
