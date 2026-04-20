#!/bin/bash
# CovFlow Analysis Wrapper Script
SCHRODINGER_RUN="/home/konstantinos/Documents/schrodinger2021-2.2/run"

if [ ! -f "$SCHRODINGER_RUN" ]; then
    echo "Error: Schrödinger run utility not found at $SCHRODINGER_RUN"
    exit 1
fi

"$SCHRODINGER_RUN" python3 "$@"
