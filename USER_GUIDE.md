# CovFlow Ultra - User Guide

CovFlow Ultra is an automated pipeline for high-fidelity covalent docking and binding analysis using Schrödinger 2021-2.

## 1. Quick Start Summary

### A. Docking Commands

**Standard Docking (Fast Ranking)**
```bash
./covflow_run.sh --csv DATA/ligands.csv --pdb DATA/1XKK.pdb --res A:797 --restype CYS
```

**Rigid Docking (Prevent Pocket Collapse)**
```bash
./covflow_run.sh --csv DATA/ligands.csv --pdb DATA/1XKK.pdb --res A:797 --restype CYS --no_min
```

**Kinetics Docking (Estimated Affinity/Ki/kinact)**
```bash
./covflow_kinetics_run.sh --csv DATA/ligands.csv --pdb DATA/1XKK.pdb --res A:797 --restype CYS
```

### B. Analysis Commands

**Automated Red Flag Analysis** *(Checks bond geometry and energy strain)*
```bash
./cov_analysis.sh ANALYSIS/red_flag_filter.py DATA/results.maegz
```

**Ligand Efficiency & MedChem Filtering** *(Calculates LE and Adjusted Scores)*
```bash
./cov_analysis.sh ANALYSIS/ligand_efficiency.py DATA/results.maegz --mw 560 --penalty 0.1
```

**H-Bond Analysis (EGFR Hinge: Met793)**
```bash
./cov_analysis.sh ANALYSIS/interaction_analysis.py DATA/results.maegz --res A:793
```

**Warhead Strain/Geometry Analysis**
```bash
./cov_analysis.sh ANALYSIS/strain_analysis.py DATA/results.maegz
```

---

## 2. Precise Pipeline Workflow

Every execution of `covflow_run.sh` follows these precise steps to ensure scientific accuracy:

1.  **Structural Sanitization**:
    - Breaks any existing covalent bonds (e.g., from native inhibitors).
    - Deletes all non-standard residues and cofactors (except those defined in `COFACTORS`).
    - **Intelligent Grid Centering**: Automatically identifies the native ligand in the pocket and calculates its centroid to define the grid center. This ensures full coverage of the binding pocket (e.g., ATP-pocket in kinases) rather than just a residue-centered box.
    - **Hydrated Docking**: All crystal waters within $5.0 \text{\AA}$ of the reaction site are **preserved** and optimized. This maintains critical H-bond networks (e.g., water-mediated hinge binding) that are often lost in "dry" docking workflows.
2.  **Protein Preparation**:
    - Runs `PrepWizard` with `-fillsidechains` and `-samplewater`.
    - Automatically executes **ProtAssign** and **Epik** to optimize the protonation state of the target residue (Cys, Ser, Lys, etc.).
3.  **Ligand Preparation & Site Selection**:
    - Runs `LigPrep` to generate 3D conformers and tautomers.
    - **Multiple Warhead Handling**: Scans for Michael acceptors. If a ligand has multiple identical warheads, the program intelligently selects the one **closest to the receptor's pocket center**.
4.  **Covalent Docking**:
    - Generates a Glide grid and executes `covalent_docking` in `leadopt` mode.
    - Applies a default distance constraint to the nucleophile.
5.  **Automated Red Flag Filtering**:
    - Scans every docked pose for physical anomalies and flags them for manual review.

---

## 3. Post-Docking "Red Flag" Filters

Results are automatically scanned (and can be run manually via `red_flag_filter.py`) to catch common covalent docking artifacts:

- **Convergence Failure (`CONV_FAIL`)**: Flagged if Prime minimization failed, indicating a major steric clash.
- **Bond Length Outliers (`BAD_BOND`)**: Flagged if the covalent bond distance is $> 2.2 \text{\AA}$ or $< 1.6 \text{\AA}$. (Ideal $C-S$ bond $\approx 1.8 \text{\AA}$).
- **Energy Mismatch (`ENERGY_MISMATCH`)**: Flagged if the Docking Score is very good ($< -8.0$) but the Prime Energy is highly positive ($> 100$), indicating the system is under extreme internal strain.

---

## 4. Ligand Efficiency & MedChem Filtering

### The Rationale: Why "Docking Score" Isn't Enough
In drug discovery, the best docking score doesn't always indicate the best drug candidate. Larger molecules often achieve high scores simply by filling space and making more non-specific contacts. This **Size Bias** can lead to prioritizing "greasy," over-sized molecules that have poor pharmacological properties.

### The Solution: Ligand Efficiency (LE)
To find the most "potent for its size" inhibitors, we use **Ligand Efficiency (LE)**:
>[!IMPORTANT]
> **LE = Docking Score / Number of Heavy Atoms**

By normalizing the score to the number of non-hydrogen atoms, we can objectively compare a small fragment to a large lead. Molecules with a high (more negative) LE are considered more "efficient" and are prioritized for optimization.

### Medicinal Chemistry Filters
To ensure the pipeline produces realistic, drug-like results, CovFlow implements:

1.  **MW Filtering ($< 560 \text{ Da}$)**:
    - Automatically flags molecules exceeding the standard molecular weight limit for oral bioavailability.
    - This "grounds" the docking results in the reality of medicinal chemistry.

2.  **Penalizing Flexibility (Adjusted Score)**:
    - **The Problem**: Highly flexible molecules can adopt many conformations. While one conformer might dock well, the entropic cost of "freezing" all those rotatable bonds into the binding pose is very high.
    - **The Solution**: CovFlow calculates an **Adjusted Score** by adding a penalty (default $0.1 \text{ kcal/mol}$) for each **Rotatable Bond**.
    - **Result**: Rigid, pre-organized molecules (like macrocycles or fused systems) are ranked higher than floppy, linear chains even if their raw scores are similar.

**Example Table Output (`medchem_analysis.csv`):**
| Ligand | Score | MW | HA | RB | LE | AdjScore | Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| LN2057 | -9.37 | 476.6 | 34 | 8 | -0.276 | -8.57 | PASS |
| Osimertinib | -6.33 | 500.6 | 37 | 6 | -0.171 | -5.73 | PASS |

---

## 5. Science & Accuracy Details

### A. Sulfur Charge & Protonation Audit
When docking to Cysteine, the script ensures the nucleophile is ready for reaction. Unlike standard automated prep, CovFlow:
1.  **Deletes the Hydrogen** from the nucleophilic atom.
2.  **Audits Formal Charges**: Updates the Sulfur formal charge to `-1` and fixes atom types for OPLS force field consistency.
3.  **Deprotonation State**: Ensures the residue is in the correct state for a Michael Addition, preventing Glide from rejecting valid poses due to valence errors.

### B. Pocket Collapse Prevention (`--no_min`)
For many PDB structures (like EGFR 7UKV), standard side-chain minimization can cause the binding pocket to incorrectly collapse.
- **Use `--no_min`** when you need to preserve the native pocket geometry (e.g., cross-docking or benchmarking).
- **Default (Minimization)**: Good for induced-fit scenarios where the pocket needs to "heal" after a native ligand is removed.

### C. Handling Displaceable Waters (`--soften`)
If you are docking large ligands that you believe should **displace** a pocket water, use the `--soften` flag:
- **`--soften 0.85`**: Reduces the receptor's Van der Waals radii to 85%. This allows ligands to overlap slightly with waters/sidechains during docking.
- **Minimization**: The subsequent Prime minimization step will then resolve the clash by moving the water or the ligand to the most energetically favorable hydrated state.

### D. Grid Sensitivity
If the nucleophile is found to be far from the native pocket center ($> 4.0\text{\AA}$), the script automatically expands the docking **INNERBOX** to ensure the entire pocket is reachable by the ligand.

- **`DATA/`** : Place input CSVs (must have `smiles` and `name`) and protein PDBs here.
- **`BIN/`** : Contains the engine scripts (`covflow.py`, `warhead_detector.py`).
- **`SCRATCH_*/`**: Temporary working directories created per-job to avoid collisions.
- **`RESULTS/`** : Final summarized output.

---

## 6. Maintenance & Cleaning

To clean up your workspace and remove all temporary scratch files, run:
```bash
./clean.sh
```
