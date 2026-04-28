#!/bin/bash
# CovFlow Workspace Cleanup Utility

echo "🧹 Cleaning CovFlow Workspace..."

# List of folders that contain deletable intermediate data
SCRATCH_FOLDERS=$(ls -d SCRATCH* 2>/dev/null)

if [ -z "$SCRATCH_FOLDERS" ]; then
    echo "No scratch folders found. Workspace is already clean."
else
    echo "The following folders will be removed:"
    echo "$SCRATCH_FOLDERS"
    
    # Optional: Confirm before deletion
    read -p "Are you sure you want to delete these folders? (y/n): " confirm
    if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
        rm -rf $SCRATCH_FOLDERS
        echo "✅ Cleanup complete. All intermediate files removed."
    else
        echo "Operation cancelled."
    fi
fi
