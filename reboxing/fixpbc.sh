#!/bin/bash

input_gro="cleaned9.gro"
output_gro="fixed9.gro"
tpr="cleaned9.tpr"
mdp="dummy.mdp"
top="A51.top"

echo "[INFO] Tworzę tymczasowy dummy.mdp..."
cat > "$mdp" <<EOF
integrator  = steep
nsteps      = 10
emtol       = 1000.0
EOF

echo "[INFO] Generuję $tpr z $input_gro i $top..."
gmx grompp -f "$mdp" -c "$input_gro" -p "$top" -o "$tpr" -maxwarn 1

# Krok 1: Wyciągnięcie liganda (resid 25), centrowanie, usuwanie PBC
echo "[INFO] Wyciągam i centrowuję ligand (resid 25)..."
echo "25" | gmx trjconv -s "$tpr" -f "$input_gro" -o ligand_raw.gro -pbc mol -center -ur compact <<EOF
0
0
EOF

# Krok 2: Przesunięcie liganda o 3 Å w osi X
echo "[INFO] Przesuwam ligand o 3 Å w osi X i zapisuję do $output_gro..."
vmd -dispdev text -eofexit << VMD_SCRIPT
mol new ligand_raw.gro type gro
set sel [atomselect top "all"]
\$sel moveby {3.0 0.0 0.0}
\$sel writegro $output_gro
quit
VMD_SCRIPT

# Sprzątanie
echo "[INFO] Czyszczę pliki tymczasowe..."
rm -f "$mdp" "$tpr" ligand_raw.gro

echo "[DONE] Ligand został przesunięty i zapisany jako $output_gro"

