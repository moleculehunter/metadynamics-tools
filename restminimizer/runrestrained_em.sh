#!/bin/bash

# Zapytaj użytkownika o nazwę układu
read -p "Podaj nazwę układu (np. A52up): " system

# Usuń dwa ostatnie znaki z nazwy systemu
short_system="${system::-2}"

echo "Pełna nazwa układu: $system"
echo "Topologia: $short_system"

# Faza 1: Generowanie plików .tpr
for i in {0..9}; do
    echo "Generowanie ${system}${i}em.tpr..."
    
    gmx grompp -f em.mdp -c ${system}${i}.gro -r ${system}${i}.gro -p ${short_system}.top -n indexposre.ndx -o topol.tpr
    
    mv topol.tpr ${system}${i}em.tpr

    echo "Plik ${system}${i}em.tpr wygenerowany!"
done

# Pauza przed fazą 2
echo "Wszystkie pliki .tpr zostały wygenerowane. Sprawdzamy zawartość folderu..."
ls -l  # Wyświetl listę plików w katalogu 

# Wymuś pauzę – użytkownik musi nacisnąć ENTER, aby kontynuować
read -p "Naciśnij ENTER, aby przejść do minimalizacji wszystkich walkerów..."

# Faza 2: Uruchamianie mdrun
for i in {0..9}; do
    echo "Uruchamianie mdrun dla ${system}${i}em.tpr..."
    
    gmx mdrun -deffnm ${system}${i}em

    echo "Zakończono mdrun dla ${system}${i}em!"
done

echo "Wszystkie symulacje zakończone!"

