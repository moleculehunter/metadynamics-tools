#!/bin/bash

# Find all .gro files that start with "walker_walker" in the current directory
files=$(ls *dn*.gro 2>/dev/null)

# Check if any matching files exist
if [ -z "$files" ]; then
    echo "No walker_walker*.gro files found in this directory!"
    exit 1
fi

# Create a VMD script to load them
echo "# VMD script to load all walker_walker*.gro files" > load_walkers.vmd

# Initialize color index
color_id=0

for file in $files; do
    echo "mol new $file" >> load_walkers.vmd
    echo "mol delrep 0 top" >> load_walkers.vmd
    echo "mol selection \"resname A51 or resname A52 or resname B40 or resname C40\"" >> load_walkers.vmd
    echo "mol representation Lines" >> load_walkers.vmd
    echo "mol color ColorID $color_id" >> load_walkers.vmd 
    echo "mol addrep top" >> load_walkers.vmd

    # Cycle through VMD color IDs (0-32), reset after reaching max
    color_id=$(( (color_id + 1) % 33 ))
done

# Run VMD with the generated script
vmd -e load_walkers.vmd

