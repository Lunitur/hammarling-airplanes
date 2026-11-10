# Hammarlingov algoritam

Radna verzija projekta iz numeričke linearne algebre. Cilj je implementirati
postupak iz rada Sorensen–Zhou (2003) za Lyapunovljevu jednadžbu

```text
A*X + X*A' + B*B' = 0,    X = U*U'
```

`U` je gornjetrokutasti faktor, a `A` mora biti stabilna.
Trenutačna implementacija koristi kompleksnu aritmetiku i za realne ulaze.

Napravljeno:

- Sylvesterova jednadžba u realnom Schurovu obliku, shema 2-solve.
- Kompleksni Hammarlingov postupak i povratak faktora u izvornu bazu.
- Mali primjer, računanje reziduala i usporedba s paketom `MatrixEquations`.
- Osnovne provjere na manjim matricama.

Preostaje:

- Realna varijanta Hammarlingova algoritma.
- Sheme 1-solve i E-solve za Sylvesterovu jednadžbu.
- Više numeričkih pokusa, mjerenje vremena i zaključak.
- Primjena na mrežu aerodroma.

## Pokretanje

```sh
nix develop
julia --project=. primjer.jl
julia --project=. tests/runtests.jl
```

U okruženju s već instaliranom Julijom prvo pokrenuti
`julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

Kod je u `sorensen_zhou.jl`, a primjer i bilješke u `primjer.jl`.
Bilježnica `Sorensen-Zhou-Lyap.ipynb` sadrži iste ćelije kao `primjer.jl`.
Može se otvoriti naredbom `jupyter notebook Sorensen-Zhou-Lyap.ipynb`;
odabrati kernel **Julia Hammarling**, koji postavlja `nix develop`.
Nakon izmjene primjera obnoviti bilježnicu:

```sh
nix shell nixpkgs#python3 --command python3 scripts/build_notebook.py
```

Ovo je početna implementacija za male guste matrice, bez zaštitnog skaliranja.
Sylvesterova jednadžba zahtijeva da `A1` i `-A2` nemaju zajedničkih svojstvenih
vrijednosti. Za loše uvjetovane probleme još treba provjeriti ponašanje algoritma.
Izvorni rad je u `Sorensen-Zhou-Lyap.pdf`.

## Kompilacija seminara

Zasebno LaTeX okruženje dostupno je naredbom `nix develop .#latex`.
Za izravnu kompilaciju iz korijena projekta:

```sh
nix develop .#latex --command latexmk -pdf -interaction=nonstopmode -halt-on-error -cd seminar/Hammarlingova_metoda_seminar.tex
```

PDF nastaje u `seminar/Hammarlingova_metoda_seminar.pdf`. Okruženje uključuje
podršku za hrvatski jezik i jezični poslužitelj `texlab`.
