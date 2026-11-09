# %% [markdown]
# # Širenje poremećaja u mreži aerodroma
# Primjena Sorensen–Zhou / Hammarlingova algoritma na dataset
# [2015 Flight Delays and Cancellations (Kaggle, US DOT)](https://www.kaggle.com/datasets/usdot/flight-delays).
#
# **Pitanje:** kako mjesto nastanka poremećaja i brzina oporavka utječu na
# varijabilnost cijele mreže? Stvarni letovi određuju veze, a dinamiku zadajemo
# pojednostavljenim difuzijskim modelom. Ne procjenjujemo niti prognoziramo stvarna kašnjenja.
#
# Pokrenite ćelije redom u kernelu **Julia Hammarling** iz ovog projekta.
# Izvor bilježnice je `flight_delays_lyapunov.jl`; pomoćne funkcije su u `flight_network.jl`.

# %%
using LinearAlgebra
using BenchmarkTools
import MatrixEquations

project_root = isfile(joinpath(pwd(), "flight_network.jl")) ? pwd() : @__DIR__
include(joinpath(project_root, "flight_network.jl"))
using .FlightNetwork

module AirportLyapunov
    include(joinpath(@__DIR__, "sorensen_zhou.jl"))
end
using .AirportLyapunov: lyapunov_factor, relative_residual

# %% [markdown]
# ## 1. Podaci i parametri
# Preuzmite i raspakirajte `flights.csv` i `airports.csv` u `data/flight-delays/`.
# Alternativno, uključite `download_data`: preuzimanje cijele Kaggle ZIP arhive
# može potrajati i zauzima nekoliko stotina MB; lokalne CSV datoteke koriste se ponovno.
# Ako javni API traži prijavu ili odbije zahtjev, preuzmite ZIP ručno preko poveznice iznad.
# `airlines.csv` nije potreban. Podaci se ne dodaju u Git.
#
# Početna ćelija učitava izvornu implementaciju te pokreće njezin primjer i provjere.
# Ovdje ne kopiramo algoritam, pa promjene u `sorensen_zhou.jl` vrijede i za ovu analizu.

# %%
data_directory = get(ENV, "FLIGHT_DELAYS_DATA", joinpath(project_root, "data", "flight-delays"))
download_data = false
month = 1
number_of_airports = 50
alpha = 1.0
gamma = 0.2
sigma = 1.0
background_noise = 0.01

ensure_dataset(data_directory; download=download_data)

# %% [markdown]
# ## 2. Čišćenje i agregiranje
# `flights.csv` čitamo redak po redak i zadržavamo samo odabrani mjesec.
# Izbacujemo otkazane letove, nepoznate/nedostajuće oznake i letove sa istim polazištem i odredištem.
# Prihvaćamo samo aerodrome s valjanim koordinatama iz `airports.csv`.
# Neke oznake aerodroma u datasetu su numeričke; ne nagađamo njihovo preslikavanje,
# nego broj odbačenih redaka izričito prikazujemo. To je posebno važno pri promjeni mjeseca.
#
# Rang aerodroma određuje ukupan broj dolazaka i odlazaka u cijelom prihvaćenom
# mjesečnom uzorku. Zatim zadržavamo samo veze između prvih `number_of_airports` čvorova.

# %%
airports = load_airports(joinpath(data_directory, "airports.csv"))
route_data = load_routes(joinpath(data_directory, "flights.csv"), airports; month)
network = build_network(route_data.routes, airports; n=number_of_airports)
display(result_table([
    (mjera="Mjesec", vrijednost=month),
    (mjera="Redaka u mjesecu", vrijednost=route_data.selected),
    (mjera="Otkazani letovi", vrijednost=route_data.cancelled),
    (mjera="Nevaljani ili nepoznati podaci", vrijednost=route_data.invalid),
    (mjera="Prihvaćeni letovi prije izbora čvorova", vrijednost=route_data.accepted),
    (mjera="Aerodroma u modelu", vrijednost=length(network.codes)),
    (mjera="Letova unutar odabrane mreže", vrijednost=Int(sum(network.counts))),
]))

# %% [markdown]
# ## 3. Od mreže do stabilnog sustava
# Neka je $C_{ij}$ broj letova iz aerodroma $i$ u aerodrom $j$. Definiramo
#
# $$\widetilde W=(C+C^\mathsf T)/2,\quad
# s=\max_i\sum_j\widetilde W_{ij},\quad W=\widetilde W/s,$$
# $$D=\operatorname{diag}(W\mathbf1),\quad L=D-W,\quad A=-\alpha L-\gamma I.$$
#
# Simetrizacija zanemaruje smjer letova. Jedinstveni faktor $s$ zadržava relativne
# težine veza i daje usporedivu skalu za $\alpha$ i $\gamma$ unutar ovog eksperimenta.
# Za simetričnu nenegativnu mrežu $L\succeq0$, pa vrijedi
# $\lambda_i(A)=-\alpha\lambda_i(L)-\gamma<0$.
# Član $-\gamma I$ nužan je: sam $-L$ ima nultu svojstvenu vrijednost.
#
# Vrijeme je u proizvoljnim jedinicama modela. Parametri nisu procijenjeni iz vremena
# letova ili kašnjenja. Izdvajanje podmreže zanemaruje sve veze s ostalim aerodromima.

# %%
A = system_matrix(network.L; alpha, gamma)
source_airports = network.codes[1:2]
B = disturbance_matrix(network.codes, source_airports; sigma, background=background_noise)
display((sources=source_airports, size_A=size(A), size_B=size(B),
         largest_eigenvalue=eigmax(Symmetric(A)),
         weight_scale=network.weight_scale))

# %% [markdown]
# ## 4. Stacionarna kovarijanca i Gramian
# Promatramo linearni stohastički sustav
# $$dx=A x\,dt+B\,d\mathcal W_t.$$
# $x_i$ je apstraktno odstupanje od ravnoteže u aerodromu $i$, a $\mathcal W_t$
# sadrži neovisne standardne Wienerove procese. Za stabilan $A$ kovarijanca zadovoljava
# $$AX+XA^\mathsf T+BB^\mathsf T=0.$$
# To je ujedno Gramian upravljivosti determinističkog sustava $\dot x=Ax+Bu$.
# [Referenca za Gramian i Lyapunovljevu jednadžbu](https://www.mathworks.com/help/control/ref/statespacemodel.gram.html).
#
# Svaki odabrani izvor daje stupac $\sigma e_i$ u $B$. Dodajemo i slabu nezavisnu
# pozadinsku pobudu $\eta I$, s `η = background_noise`, tako da
# $$BB^\mathsf T=\sigma^2\sum_{i\in S}e_ie_i^\mathsf T+\eta^2 I.$$
# To izbjegava singularnu kovarijancu kod nedosegljivih modova. Pozadinska pobuda
# jednaka je u svim usporedbama; može se isključiti postavljanjem vrijednosti na nulu.
# Tada singularnost ili loša uvjetovanost mogu otežati faktorizaciju.
#
# $X_{ii}$ mjeri modeliranu varijancu, $X_{ij}$ kovarijancu, a $\operatorname{tr}(X)$
# ukupnu varijabilnost mreže. Te vrijednosti nemaju jedinicu minuta stvarnog kašnjenja.

# %%
U_real = lyapunov_factor(A, B)
U_complex = lyapunov_factor(A, B; arithmetic=:complex)
U_reference = MatrixEquations.plyapc(A, B)
X = Matrix(U_real * U_real')
X_reference = U_reference * U_reference'
solver_results = [
    (solver=name, residual=relative_residual(A, B, U),
     relative_error=norm(U * U' - X_reference) / norm(X_reference))
    for (name, U) in [("Sorensen–Zhou real", U_real),
                      ("Sorensen–Zhou complex", U_complex),
                      ("MatrixEquations.plyapc", U_reference)]
]
display(result_table(solver_results))
@assert maximum(row.residual for row in solver_results) < 1e-10
@assert maximum(row.relative_error for row in solver_results) < 1e-8

# %% [markdown]
# ## 5. Gdje su poremećaji najizraženiji?
# Karta prikazuje geografski raspored bez kartografske podloge; debljina veze ovisi
# o broju letova. Tamniji i veći čvor znači veću varijancu. Prelaskom miša preko
# čvora dobivaju se naziv i vrijednost. Skala boje relativna je tom prikazu.
# Toplinska karta pokazuje apsolutne kovarijance za prvih 20 aerodroma po prometu,
# dok tablica rangira sve čvorove prema varijanci. Grafike su SVG, bez dodatnih servisa.

# %%
display(airport_map(network, diag(X)))
display(covariance_heatmap(network.codes, X))
ranking = sortperm(diag(X); rev=true)
display(result_table([
    (airport=network.codes[i], name=network.airports[i].name,
     monthly_traffic=network.traffic[i], variance=round(X[i, i]; sigdigits=5))
    for i in ranking[1:min(10, length(ranking))]
]))

# %% [markdown]
# ## 6. Jednaka pobuda u prometnom i manje prometnom čvoru
# Uspoređujemo prvi i posljednji aerodrom prema prometu **unutar odabranog skupa**.
# Svaki scenarij ima jedan lokalni izvor iste jačine i jednaku pozadinsku pobudu.
# Posebno izdvajamo varijancu izvan izvornog čvora: ona opisuje koliko se učinak
# proširio na ostatak mreže. Velik lokalni poremećaj ne mora značiti velik utjecaj drugdje.

# %%
source_comparison = map([first(network.codes), last(network.codes)]) do source
    Bs = disturbance_matrix(network.codes, [source]; sigma, background=background_noise)
    Us = lyapunov_factor(A, Bs)
    Xs = Us * Us'
    index = findfirst(==(source), network.codes)
    (source=source, local_variance=Xs[index, index], total_variance=tr(Xs),
     variance_elsewhere=tr(Xs) - Xs[index, index], residual=relative_residual(A, Bs, Us))
end
display(result_table(source_comparison))

# %% [markdown]
# ## 7. Utjecaj brzine oporavka
# Mijenjamo samo $\gamma$, uz istu mrežu, normalizaciju, $\alpha$ i $B$.
# Manji $\gamma$ znači sporiji oporavak. Za ovaj simetrični model
# $$\operatorname{tr}(X)=\tfrac12\operatorname{tr}((\alpha L+\gamma I)^{-1}BB^\mathsf T),$$
# pa ukupna varijanca pada s porastom $\gamma$. To ujedno daje nezavisnu provjeru
# očekivanog ponašanja rezultata, uz provjeru reziduala.

# %%
recovery_results = map([0.02, 0.05, 0.1, 0.2, 0.5, 1.0]) do recovery_rate
    Ag = system_matrix(network.L; alpha, gamma=recovery_rate)
    Ug = lyapunov_factor(Ag, B)
    Xg = Ug * Ug'
    trace_check = tr((-Ag) \ (B * B')) / 2
    @assert isapprox(tr(Xg), trace_check; rtol=1e-8)
    (gamma=recovery_rate, total_variance=tr(Xg),
     max_variance=maximum(real, diag(Xg)), residual=relative_residual(Ag, B, Ug))
end
@assert all(diff(getproperty.(recovery_results, :total_variance)) .< 0)
display(result_table(recovery_results))
display(recovery_plot(recovery_results))

# %% [markdown]
# ## 8. Vrijeme izvođenja
# Pokretanje je opcionalno. BenchmarkTools zagrijava funkcije i mjeri ponovljene
# pozive na istim matricama; prikazujemo medijan. Učitavanje CSV-a i crtanje nisu
# uključeni. Sve tri metode računaju faktor, a Schurova dekompozicija ulazi u mjerenje.
# Simetrični $A$ ima realan spektar, pa ovaj primjer ne ispituje realne Schurove
# blokove veličine 2. Takvi testovi ostaju u izvornoj bilježnici.

# %%
run_timing = false
if run_timing
    timings = [
        (solver="Sorensen–Zhou real", trial=@benchmark(lyapunov_factor($A, $B), samples=10, evals=1)),
        (solver="Sorensen–Zhou complex", trial=@benchmark(lyapunov_factor($A, $B; arithmetic=:complex), samples=10, evals=1)),
        (solver="MatrixEquations.plyapc", trial=@benchmark(MatrixEquations.plyapc($A, $B), samples=10, evals=1)),
    ]
    display(result_table([
        (solver=row.solver, median_ms=median(row.trial).time / 1e6,
         memory_bytes=median(row.trial).memory) for row in timings
    ]))
end

# %% [markdown]
# ## Ograničenja i mogući nastavci
# - Mreža prikazuje broj letova, bez kapaciteta zrakoplova, presjedanja i rasporeda.
# - Simetrizacija i difuzija su pretpostavke modela, a ne empirijski utvrđena dinamika.
# - Pri usporedbi mjeseci normalizacija se ponovno računa; za čistu usporedbu
#   sezonskih promjena trebalo bi zadržati iste čvorove i isti faktor normalizacije.
# - Pozadinski šum osigurava pozitivno definitno rješenje, ali mijenja model.
# - Nakon promjene parametara pokrenite sve ćelije od početka da prikazi budu usklađeni.
#
# Daljnji korak bio bi procijeniti dinamiku iz vremenskih nizova kašnjenja i validirati
# model na izdvojenim podacima. Diskretni model zahtijeva diskretnu Lyapunovljevu jednadžbu,
# koju ovdje implementirani kontinuirani solver ne rješava izravno.
