xparamber1 to skrypt automatyzujący parametryzację małych ligandów w polu siłowym AMBER
przy użyciu pakietu AmberTools23

Przed użyciem programu należy zainstalować pakiet AmberTools23
conda create --name AmberTools23
conda activate AmberTools23 
conda install -c conda-forge ambertools=23



1. W programie avogadro stwórz cząsteczkę

2. Zapisz cząsteczkę w pliku .mol2 (typ atomów - sybyl)
   Nazwij cząsteczkę używając TYLKO 3 znaków, a następnie dopisz _syb (np. A01_syb.mol2)

3. Stwórz nowy folder i przekopiuj do niego cząsteczkę stworzoną w punkcie 2

4. Do folderu wypakuj archiwum
   (chodzi o pliki tleap.in oraz xparamber1.sh)
   #UWAGA! w wyniku działania programu, plik tleap.in zostaje zmodyfikowany i nie nadaje się do ponownego użycia
   #Można kopiować plik tleap.in zamiast wypakowywać archiwum za każdym razem, jednak trzeba mieć na uwadze,
   aby zawsze kopiować plik który nie jest "zużyty"
   
5. Upewnij się, że pakiet AmberTools23 jest aktywny
   > conda activate AmberTools23

6. Uruchom program wpisując:
   ./xparamber1.sh [nazwa pliku .mol2] -nc [ładunek cząsteczki]
   np. ./xparamber1 A01_syb.mol2 -nc 1
   
Program powinien utworzyć folder A01.amb2gmx, wewnątrz którego znajdują się pliki .gro, .top, .itp, .mdp i rungmx.sh

UWAGA!
Bardzo częstym problemem może być komunikat: "Possible open valence".
Wynika on z (niezbyt umiejętnego) automatycznego przypisywania typów wiązań przez program avogadro
Należy wtedy otworzyć plik A01_syb.mol2 w edytorze tekstowym (np. vim) i zmienić rodzaje wiązań
(Zazwyczaj trzeba zmienić wiązania aromatyczne 'ar' na pojedyncze '1')
(Komunikat 'Possible open valence' wskaże na którym atomie występuje owy problem).
