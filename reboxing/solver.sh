#!/bin/bash

topology="A51.top"
backup_top="A51.top.bak"
water_template="spc216.gro"
solvation_prefix="solvcube"
box_size="8.0"

echo "[INFO] Tworzę backup topologii..."
cp "$topology" "$backup_top"

echo "[INFO] Tworzę cubic box, centrowanie i zalewanie wodą..."

for file in cleaned*.gro; do
    num=$(echo "$file" | sed -E 's/^cleaned([0-9]+)\.gro/\1/')
    echo "[INFO] → Przetwarzam $file → $solvation_prefix${num}.gro"

    gmx editconf -f "$file" -o tmp_centered.gro -box $box_size $box_size $box_size -c > /dev/null 2>&1
    gmx solvate -cp tmp_centered.gro -cs $water_template -o "${solvation_prefix}${num}.gro" -p "$topology" > /dev/null 2>&1
done

# Krok: policz liczbę cząsteczek SOL w każdym pliku
declare -A sol_counts
min_sol=9999999

for file in ${solvation_prefix}*.gro; do
    sol_atoms=$(grep -c "^.\{5\}SOL" "$file")
    sol_mols=$(expr $sol_atoms / 3)
    sol_counts[$file]=$sol_mols
    if [ "$sol_mols" -lt "$min_sol" ]; then
        min_sol=$sol_mols
    fi
done

echo "[INFO] Minimalna liczba cząsteczek SOL: $min_sol"

# Przycinanie .gro do min_sol
for file in ${solvation_prefix}*.gro; do
    num=$(echo "$file" | sed -E 's/^solvcube([0-9]+)\.gro/\1/')
    box_line=$(tail -n 1 "$file")
    total_atoms=$(wc -l < "$file")
    atom_count_line=$((total_atoms - 1))
    atom_lines=$(head -n "$atom_count_line" "$file")

    # Zostaw nagłówek + wszystkie nie-SOL + min_sol * 3 atomów SOL
    echo "[INFO] Przycinam $file ${sol_counts[$file]} → ${min_sol} cząsteczek SOL"

    {
        echo "$atom_lines" | head -n 2  # nagłówek
        echo "$atom_lines" | tail -n +3 | awk -v limit=$((min_sol * 3)) '
            BEGIN { sol=0 }
            {
                if ($0 ~ /^.{5}SOL/) {
                    if (sol >= limit) next
                    sol++
                }
                print
            }
        '
    } > trimmed_tmp.gro

    # Ustaw poprawną liczbę atomów (linia 2)
    natoms=$(wc -l < trimmed_tmp.gro)
    natoms=$((natoms - 2))
    sed -i "2s/.*/$natoms/" trimmed_tmp.gro

    # Dodaj z powrotem box
    echo "$box_line" >> trimmed_tmp.gro
    mv trimmed_tmp.gro "$file"
done

# Ustaw poprawnie A51.top
echo "[INFO] Ustawiam SOL $min_sol w pliku topologii..."

awk -v total=$min_sol '
BEGIN { in_mol = 0 }
{
    if ($0 ~ /^\[ *molecules *\]/) {
        in_mol = 1
        print
        next
    }
    if (in_mol) {
        if ($0 ~ /^\[/) {
            print "SOL", total
            in_mol = 0
        }
        if ($0 !~ /^[[:space:]]*$/ && $1 != "SOL") print
        next
    }
    print
}
END {
    if (in_mol) print "SOL", total
}
' "$backup_top" > "$topology"


echo "[OK] Topologia zaktualizowana."
echo "[DONE] Przygotowanie układów z wodą zakończone."

