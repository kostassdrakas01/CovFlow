#!/bin/bash
# CovFlow Kinetics Wrapper Script
SCHRODINGER_RUN="/home/konstantinos/Documents/schrodinger2021-2.2/run"
SCRIPT_PATH="/home/konstantinos/Documents/CovFlow/BIN/covflow_kinetics.py"

if [ ! -f "$SCHRODINGER_RUN" ]; then
    echo "Error: Schrödinger run utility not found at $SCHRODINGER_RUN"
    exit 1
fi

"$SCHRODINGER_RUN" python3 "$SCRIPT_PATH" "$@"
