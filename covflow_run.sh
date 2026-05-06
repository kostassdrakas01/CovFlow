#!/bin/bash
# CovFlow Ultra Wrapper Script
# --- CONFIGURATION (Update these paths for new computers) ---
SCHRODINGER_RUN="/home/konstantinos/Documents/schrodinger2021-2.2/run"
SCRIPT_PATH="/home/konstantinos/Documents/CovFlow/BIN/covflow.py"
# ----------------------------------------------------------

if [ ! -f "$SCHRODINGER_RUN" ]; then
    echo "Error: Schrödinger run utility not found at $SCHRODINGER_RUN"
    exit 1
fi
# Default arguments (can be overridden by environment variables)
CSV_FILE=${CSV_FILE:-"DATA/ligands.csv"}
PDB_FILE=${PDB_FILE:-"DATA/7UKV.pdb"}
TARGET_RES=${TARGET_RES:-"A:797"}
TARGET_RESTYPE=${TARGET_RESTYPE:-"CYS"}
RXN_TYPE=${RXN_TYPE:-"michael_addition"}
DIST_CONSTRAINT=${DIST_CONSTRAINT:-"2.5"}

# Parse arguments manually to allow overriding defaults
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --csv) CSV_FILE="$2"; shift ;;
        --pdb) PDB_FILE="$2"; shift ;;
        --res) TARGET_RES="$2"; shift ;;
        --restype|--resname) TARGET_RESTYPE="$2"; shift ;;
        --rxn) RXN_TYPE="$2"; shift ;;
        --dist) DIST_CONSTRAINT="$2"; shift ;;
        --host) HOST="$2"; shift ;;
        --no_min) NO_MIN="--no_min" ;;
        --soften) SOFTEN="--soften $2"; shift ;;
        --center_res) CENTER_RES="--center_res $2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo "Starting CovFlow..."
"$SCHRODINGER_RUN" python3 "$SCRIPT_PATH" \
    --csv "$CSV_FILE" \
    --pdb "$PDB_FILE" \
    --res "$TARGET_RES" \
    --restype "$TARGET_RESTYPE" \
    --rxn "$RXN_TYPE" \
    --dist "$DIST_CONSTRAINT" \
    $NO_MIN $SOFTEN $CENTER_RES

