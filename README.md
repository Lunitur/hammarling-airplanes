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
Koristi se službeni [OPDI v0.0.2 (EUROCONTROL PRC / OpenSky)](https://www.opdi.aero/data),
zadano **siječanj 2025.** Bilježnica automatski preuzima mjesečni Parquet i
[OurAirports referentni CSV](https://raw.githubusercontent.com/davidmegginson/ourairports-data/7dd521a50a421171db1cb69f87657a60b48561c7/airports.csv)
u `data/europe-flights/`. Postojeće datoteke koriste se ponovno; `download_data = false`
isključuje preuzimanje. Alternativnu lokaciju određuje varijabla `FLIGHT_DELAYS_DATA`.
Za ručno preuzimanje spremite
[`flight_list_202501.parquet`](https://www.eurocontrol.int/performance/data/download/OPDI/v002/flight_list/flight_list_202501.parquet)
(oko 34 MB) i navedeni CSV kao `airports.csv`. Parametri `year` i `month` određuju
Parquet datoteku; dostupne mjesece pogledajte na [OPDI portalu](https://www.opdi.aero/flight-list-data).
OurAirports snapshot je zaključan na commit `7dd521a50a421171db1cb69f87657a60b48561c7`.
Lokalne datoteke čuvajte za ponavljanje rezultata jer se OPDI datoteke mogu revidirati.

Pokrenite ćelije redom u kernelu `Julia Hammarling`. Bilježnica gradi stabilan
difuzijski model za 50 aerodroma, prikazuje kovarijance i uspoređuje izvore
poremećaja, brzine oporavka i solvere. Podaci određuju mrežu; dinamika je
pretpostavljena i ne predstavlja validiranu prognozu stvarnih kašnjenja.

OPDI je odabran jer model treba **broj opaženih letova između parova aerodroma**,
a ne mjerena kašnjenja. Učitavaju se samo `dof`, `adep` i `ades`, po grupama redaka.
Oba kraja veze moraju odgovarati OurAirports ICAO oznakama (`icao_code`), imati
valjane koordinate, `continent = EU` i ne smiju biti označeni kao zatvoreni.
Ova geografska definicija uključuje UK i Švicarsku, ali isključuje Cipar i azijski
dio Turske; nije jednaka Europskoj uniji ili području EUROCONTROL-a.
Nepoznate oznake, letovi izvan tog skupa i povratak u isto polazište odbacuju se
i broje u izvještaju. Predviđene oznake `adep_p`/`ades_p` ne koriste se.
Uzorak ovisi o ADS-B pokrivenosti, nije potpuni popis prometa i uključuje i
nekomercijalne letove. Nema podataka o otkazivanju ili kašnjenju prema rasporedu;
`first_seen`/`last_seen` ne tumače se kao vremena na stajanci.

Ostali razmotreni izvori:

- [Hugging Face European Flight Delays](https://huggingface.co/datasets/345rf4gt56t4r3e3/flight-delays-europe-2023-2025):
  `delay_min` je trajanje leta minus medijan rute/mjeseca, ne kašnjenje prema redu
  letenja. Autor izričito upozorava da dataset nije točan istraživački izvor.
- [EUROCONTROL CSV](https://ansperformance.eu/csv/): službene agregacije korisne
  za buduću usporedbu kašnjenja, ali ne daju pojedinačne letove za izgradnju ove mreže.
- [UK CAA punctuality](https://www.caa.co.uk/data-and-analysis/uk-aviation-market/flight-punctuality/uk-flight-punctuality-statistics/):
  podaci o točnosti za UK, užeg geografskog opsega od ovog europskog primjera.

Izvore navodite kao EUROCONTROL PRC / OpenSky Network (OPDI) i OurAirports;
[OPDI opis uporabe](https://www.opdi.aero/about) usmjeren je na nekomercijalnu uporabu,
a [OurAirports](https://ourairports.com/data/) objavljuje podatke u javnoj domeni.
Preuzeti podaci nisu dio ovog Git repozitorija.

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
