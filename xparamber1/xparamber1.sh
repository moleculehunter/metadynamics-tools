#!/bin/bash

# Initialize variables
molfile=""
charge=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -nc) charge="$2"; shift ;;  # Get the value of the charge after the -nc flag
        *.mol2) molfile="$1" ;;     # Get the .mol2 file name
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if both molfile and charge are provided
if [ -z "$molfile" ] || [ -z "$charge" ]; then
    echo "Usage: $0 XXX_syb.mol2 -nc [number]"
    exit 1
fi

# Extract molecule name from the filename (first 3 characters)
name=$(basename "$molfile" | cut -c1-3)

# Run antechamber and parmchk2 with the extracted name and charge
echo "Running antechamber and parmchk2..."
antechamber -i "${molfile}" -fi mol2 -o "${name}.mol2" -fo mol2 -c bcc -s 2 -at gaff -nc "$charge"
parmchk2 -i "${name}.mol2" -f mol2 -o "${name}.frcmod"
antechamber -i "${name}.mol2" -fi mol2 -o "${name}.prepi" -fo prepi -c bcc -s 2 -at gaff -nc "$charge"

# Edit the 5th line of the XXX.prepi file
echo "Editing ${name}.prepi file..."
sed -i "5s/^.../${name}/" "${name}.prepi"

# Copy tleap.in from the library and modify it
echo "Copying and modifying tleap.in..."
cp /path/to/xparamber/library/tleap.in ./tleap.in
sed -i "s/XXX/${name}/g" tleap.in

# Run tleap
echo "Running tleap..."
tleap -f tleap.in

# Run acpype
echo "Running acpype..."
acpype -p "${name}.prmtop" -x "${name}.inpcrd"

echo "Process complete!"

