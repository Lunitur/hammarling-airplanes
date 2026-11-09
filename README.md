# Hammarlingov algoritam

Implementacija rada Sorensen–Zhou (2003), priloženog u `Sorensen-Zhou-Lyap.pdf`:

- [Jupyter notebook](Sorensen-Zhou-Lyap.ipynb) s objašnjenjima na hrvatskom.
- [Julia kod po ćelijama](sorensen_zhou.jl), iz kojeg je notebook generiran.
- Sylvesterove sheme `:one`, `:two` i `:eigen` (odjeljak 2 rada).
- Realni i kompleksni Hammarlingov algoritam (algoritmi 3.1 i 3.2).
- Primjer, provjere prema `MatrixEquations` i opcionalna mjerenja vremena.

Glavno sučelje je `U = lyapunov_factor(A, B; arithmetic=:real)` za jednadžbu
`A*X + X*A' + B*B' = 0`, uz **`X = U*U'`**. Za kompleksnu aritmetiku odaberite
`:complex`. `sylvester_sz(A1, A2, D; scheme=:two)` rješava
`A1*X + X*A2 + D = 0`. Ulazne matrice ostaju nepromijenjene.

Za izvršavanje svih primjera i provjera:

```sh
nix develop --command julia --project=. sorensen_zhou.jl
```

Notebook pokrenite naredbama:

```sh
nix develop
jupyter notebook Sorensen-Zhou-Lyap.ipynb
```

Primjena na mrežu aerodroma nalazi se u `Flight-Delays-Lyapunov.ipynb`
(izvor: `flight_delays_lyapunov.jl`, pomoćne funkcije: `flight_network.jl`).
S [Kagglea](https://www.kaggle.com/datasets/usdot/flight-delays) raspakirajte
`flights.csv` i `airports.csv` u `data/flight-delays/`, ili u bilježnici postavite
`download_data = true`. Alternativnu lokaciju određuje varijabla `FLIGHT_DELAYS_DATA`.
Pokrenite ćelije redom u kernelu `Julia Hammarling`. Bilježnica gradi stabilan
difuzijski model za 50 aerodroma, prikazuje kovarijance i uspoređuje izvore
poremećaja, brzine oporavka i solvere. Podaci određuju mrežu; dinamika je
pretpostavljena i ne predstavlja validiranu prognozu stvarnih kašnjenja.

Obnovite novu bilježnicu iz Julia izvora:

```sh
nix shell nixpkgs#python3 --command python3 scripts/build_notebook.py flight_delays_lyapunov.jl Flight-Delays-Lyapunov.ipynb
```

Provjere učitavanja i izgradnje mreže:

```sh
nix develop --command julia --project=. tests/flight_network.jl
```

U notebooku odaberite kernel `Julia Hammarling` i pokrenite ćelije redom.
Pri ulasku u okruženje IJulia se automatski instalira u zasebno Julia okruženje
i registrira kernel koji koristi projektne ovisnosti. Prva instalacija zahtijeva
pristup internetu. Kernel možete odabrati i u VS Codeu.

`run_benchmarks()` pokreće usporedbu vremena. Implementacija je obrazovna,
za male i srednje guste matrice; ne uključuje automatsko skaliranje protiv
overflowa. Ograničenja uvjetovanosti i korekcije pseudokoda objašnjeni su u notebooku.

Nakon promjena Julia ćelija osvježite notebook:

```sh
nix shell nixpkgs#python3 --command python3 scripts/build_notebook.py
```

Razvojno okruženje za implementaciju Hammarlingova algoritma za rješavanje
Lyapunovljeve jednadžbe izravnim računanjem Choleskyjeva faktora rješenja.

```sh
nix develop
julia -i dev.jl
```

Pri ulasku u okruženje automatski se aktivira Julia projekt i instaliraju
ovisnosti iz `Manifest.toml`. Prvo pokretanje zahtijeva pristup internetu.
`flake.lock` zaključava Nix ovisnosti, a `Manifest.toml` Julia pakete.

`dev.jl` učitava:

- `LinearAlgebra`: Schurova i QR dekompozicija, trokutaste matrice, BLAS i LAPACK.
- `Random`: reproducibilni slučajni primjeri preko `Random.seed!`.
- `Test`: provjera reziduala i točnosti implementacije.
- `BenchmarkTools`: mjerenje vremena i alokacija.
- `MatrixEquations`: referentni solveri `MatrixEquations.plyapc` i `MatrixEquations.plyapd`.

Za samu implementaciju osnovne matrične operacije daje `LinearAlgebra`.
Referentni solveri služe za usporedbu rezultata; njihova sučelja opisana su u
[dokumentaciji MatrixEquations](https://andreasvarga.github.io/MatrixEquations.jl/dev/plyapunov.html).

Pokretanje bez interaktivne ljuske:

```sh
nix develop --command julia --project=. -i dev.jl
```
