#!/bin/bash

topology="A51.top"
ion_conc="0.02"
base_topology="A51.base"
index_file="shared_index.ndx"

# Ustal liczbę plików wejściowych
N=$(ls solvcube*.gro | wc -l)
echo "[INFO] Znaleziono $N plików solvcube*.gro"

# Tworzymy plik dummy.mdp
echo "[INFO] Tworzę dummy.mdp..."
cat > dummy.mdp <<EOF
integrator  = steep
nsteps      = 10
emtol       = 1000.0
EOF

# Tworzymy bazowy plik topologii, który nie będzie modyfikowany
cp "$topology" "$base_topology"
echo "[INFO] Bazowa topologia zapisana jako $base_topology"

# Tworzymy index.ndx z grupą SOL
echo "[INFO] Tworzę index.ndx z grupą SOL..."
echo -e "keep 0\nr SOL\nname 1 SOL\nq" | gmx make_ndx -f solvcube0.gro -o "$index_file" > /dev/null 2>&1

# Pętla iterująca po plikach
for i in $(seq 0 $((N - 1))); do
    infile="solvcube${i}.gro"
    outfile="ionized${i}.gro"
    log="genion_${i}.log"

    echo "──────────────────────────────────────────────"
    echo "[DEBUG] Iteracja $i"
    echo "[DEBUG] Plik wejściowy: $infile"
    echo "[DEBUG] Kopiuję $base_topology jako $topology"

    cp "$base_topology" "$topology"

    if [ ! -f "$infile" ]; then
        echo "[WARN] Plik $infile nie istnieje – pomijam"
        continue
    fi

    echo "[DEBUG] HEAD $infile ↓"
    head -n 10 "$infile"
    echo "[DEBUG] grompp z: $infile → tmp.tpr"

    gmx grompp -f dummy.mdp -c "$infile" -p "$topology" -o tmp.tpr -maxwarn 1 > "grompp_${i}.log" 2>&1

    if [ ! -f tmp.tpr ]; then
        echo "[ERROR] grompp nie wygenerował tmp.tpr – sprawdź grompp_${i}.log"
        continue
    fi

    echo "[DEBUG] genion: tmp.tpr → $outfile"

    echo "1" | gmx genion -s tmp.tpr -o "$outfile" -p "$topology" \
        -pname K -nname CL -neutral -conc "$ion_conc" -n "$index_file" > "$log" 2>&1

    if grep -q "Replacing" "$log"; then
        grep "Replacing" "$log"
        echo "[OK] Jony dodane → $outfile"
    else
        echo "[ERROR] genion nie dodał jonów do $outfile – sprawdź $log"
    fi
done

# Zlicz jony na końcu (debug info)
echo
echo "[INFO] Sprawdzam ilość jonów w plikach wynikowych..."

for file in ionized*.gro; do
    k_count=$(awk '{r=substr($0,6,5); if (r ~ /^ *K *$/) k++} END {print k+0}' "$file")
    cl_count=$(awk '{r=substr($0,6,5); if (r ~ /^ *CL *$/) c++} END {print c+0}' "$file")
    echo "[INFO] $file zawiera: $k_count/$cl_count (K+/Cl-)"
done

rm -f tmp.tpr "$index_file"
echo "[DONE] Jonizacja zakończona. Logi: genion_*.log"

