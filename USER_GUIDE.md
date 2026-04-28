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

### B. Validation & Test

**Automated Self-Test** *(Runs Afatinib on EGFR-6JXT and verifies bond)*
```bash
$SCHRODINGER/run python3 BIN/validate_pipeline.py
```

### C. Analysis Commands

**H-Bond Analysis (EGFR Hinge: Met793)**
```bash
./cov_analysis.sh ANALYSIS/interaction_analysis.py DATA/results.maegz --res A:793
```

**Warhead Strain/Geometry Analysis**
```bash
./cov_analysis.sh ANALYSIS/strain_analysis.py DATA/results.maegz
```

---

## 2. Setup & Data Folders

- **`DATA/`** : Place input CSVs (must have `smiles` and `name`) and protein PDBs here.
- **`BIN/`** : Contains the engine scripts (`covflow.py`, `warhead_detector.py`).
- **`SCRATCH_*/`**: Temporary working directories created per-job to avoid collisions.
- **`RESULTS/`** : Final summarized output (Auto-synced to Dropbox).

---

## 3. Advanced Features & Fixes

### A. Dynamic Scratch Directories
The pipeline now creates a unique folder for every job:
`Example: SCRATCH_1XKK_A_797/`
This allows you to run multiple docking jobs in parallel without logs overwriting each other.

### B. Sulfur Charge Audit (Critical Fix)
When docking to Cysteine, the script manually deprotonates the nucleophile. Unlike standard PrepWizard, CovFlow Ultra now:
1. Deletes the Hydrogen.
2. Updates the Sulfur formal charge to `-1`.
3. Fixes atom types for OPLS force field consistency.

> [!IMPORTANT]
> This prevents Glide from rejecting poses due to "Lewis structure" or valence errors.

### C. Grid Sensitivity Warnings
If the nucleophile is far from the native pocket center (> 4.0A), the script gives a ⚠️ **WARNING** and automatically expands the docking INNERBOX to `20x20x20` to ensure proper coverage.

### D. Covalent Bond Restoration
During the 'Sanitize' phase, CovFlow detects if your protein PDB already has a covalent bond (e.g., from a co-crystallized inhibitor). It automatically:
1. Breaks the bond.
2. Restores the nucleophile to its native deprotonated state (`-1`).
3. Minimizes the neighboring shell (5A) to "heal" the pocket before new docking.

### E. Pocket Collapse Prevention (`--no_min`)
For certain protein structures, standard side-chain minimization during the preparation cycle can cause the binding pocket to incorrectly collapse or distort. Using the `--no_min` flag instructs the pipeline to bypass this minimization. You **must** use this flag when evaluating rigid cross-benchmark datasets (such as experimental EGFR inhibitors) to preserve the native pocket geometry, ensuring that computational docking scores accurately align with real-world geometric constraints and potency data.

---

## 4. Analysis Significance

1. **`interaction_analysis.py`**:
   - **High significance:** Identifies leads that maintain critical binding contacts with the hinge region. If these are missing, Ki is usually weak.

2. **`strain_analysis.py`**:
   - **High significance:** Flags "strained" warheads (torsion > 30°). Strained acrylamide warheads react slowly (poor kinact).

3. **`master_analysis.py`**:
   - Flexible tool for measuring any distance (e.g., Deep-pocket occupancy or water-bridging) using Maestro ASL.

---

## 5. Maintenance & Cleaning

To clean up your workspace and remove all temporary scratch files, run the automated cleanup utility:
```bash
./clean.sh
```

> [!NOTE]
> Keep the `BIN/` and `DATA/` folders intact during cleanup.
