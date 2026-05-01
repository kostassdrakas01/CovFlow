#!/bin/bash
# CovFlow Workspace Cleanup Utility

FORCE=false
if [[ "$1" == "-f" || "$1" == "--force" ]]; then
    FORCE=true
fi

echo "🧹 Cleaning CovFlow Workspace..."

# Find folders starting with 'scratch' (case-insensitive)
SCRATCH_FOLDERS=$(find . -maxdepth 1 -type d -iname "scratch*" | sed 's|./||' | grep -v "^\.$")

if [ -z "$SCRATCH_FOLDERS" ]; then
    echo "No scratch folders found. Checking for temporary files..."
    rm -f *.log *.maegz *.mae *.zip *.rept *.txt
    echo "✅ Workspace cleaned."
else
    if [ "$FORCE" = true ]; then
        rm -rf $SCRATCH_FOLDERS
        rm -f *.log *.maegz *.mae *.zip *.rept *.txt
        echo "✅ Cleanup complete. Folders and temporary files removed."
    else
        echo "The following folders will be removed:"
        echo "$SCRATCH_FOLDERS"
        echo ""
        echo "Note: Temporary Schrödinger files (*.log, *.maegz, etc.) will also be cleared."
        read -p "Are you sure you want to delete these? (y/n): " confirm
        if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
            rm -rf $SCRATCH_FOLDERS
            rm -f *.log *.maegz *.mae *.zip *.rept *.txt
            echo "✅ Cleanup complete. Workspace is now clean."
        else
            echo "Operation cancelled."
        fi
    fi
fi
