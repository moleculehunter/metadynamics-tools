#!/bin/bash

# Zapytaj użytkownika o nazwę układu
read -p "Podaj nazwę układu (np. A52up; B40dn...): " system

# Usuń dwa ostatnie znaki z nazwy systemu
short_system="${system::-2}"

echo "Pełna nazwa układu: $system"
echo "Ligand: $short_system"

# Ustawienie nazw plików
GRO_FILE="${system}1.gro"
NDX_FILE="indexposre.ndx"

# Sprawdzenie, czy plik .gro istnieje
if [ ! -f "$GRO_FILE" ]; then
    echo "Plik $GRO_FILE nie istnieje!"
    exit 1
fi

echo "Tworzenie pliku indeksowego..."

# Tworzenie pliku .ndx
gmx make_ndx -f "$GRO_FILE" -o "$NDX_FILE" << EOF
1 | r 23 24
name 10 DNA-K
2 & ! r 23 24
name 11 K_zew

# Kopiowanie ważnych grup, aby znalazły się na końcu
10
3
7
11
4

# Usunięcie zbędnych grup (1-11)
del 1-11

# Nadanie nazw skopiowanym grupom (które są teraz na dole)
name 1 DNA-K
name 2 ${short_system}
name 3 SOL
name 4 K
name 5 CL

# Stworzenie grupy DNA-K_Ligand, potrzebnej do nałożenia restrainów
1 | 2
name 6 POSRES

# Zakończenie programu
q
EOF

echo "Plik indeksowy $NDX_FILE został utworzony."

echo "Tworzenie pliku posre.itp..."

gmx genrestr -f "$GRO_FILE" -n indexposre.ndx -o posre.itp << EOF
6
q
EOF

echo "Plik posre.itp został utworzony."

#Następnie można odpalić skrypt runrestrained_em.sh

