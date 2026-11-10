# %% [markdown]
# # Hammarlingov algoritam — radna verzija
# Cilj je riješiti Lyapunovljevu jednadžbu
#
# $$AX + XA^* + BB^* = 0$$
#
# računanjem faktora $U$, tako da je $X = UU^*$. Matrica $A$ mora imati sve
# svojstvene vrijednosti u lijevoj poluravnini. Ovdje je napravljena kompleksna
# varijanta algoritma; realna varijanta ostaje za nastavak.
#
# Koristim rad Sorensen–Zhou (2003), priložen u `Sorensen-Zhou-Lyap.pdf`.

# %%
using LinearAlgebra
import MatrixEquations
include(joinpath(@__DIR__, "sorensen_zhou.jl"))

# %% [markdown]
# ## Sylvesterova jednadžba
# Najprije računam $A_1Y + YA_2 + D = 0$ pomoću realnog Schurova oblika.
# Blokovi mogu biti veličine 1 ili 2. Za blokove veličine 2 zasad koristim
# samo shemu 2-solve. Desnu stranu zadajem iz poznatog rješenja radi provjere.

# %%
A1 = [-1.0 4 1; -2 -1 0; 0 0 -3]
A2 = [-2.0 3; -1 -2]
Y_known = [1.0 2; 0 -1; 3 1]
D = -(A1 * Y_known + Y_known * A2)
Y = sylvester_sz(A1, A2, D)
display((greska=norm(Y - Y_known), rezidual=norm(A1 * Y + Y * A2 + D)))

# %% [markdown]
# ## Kompleksni Hammarlingov postupak
# Schurova dekompozicija pretvara $A$ u gornjetrokutastu matricu $R$.
# Faktor računam od zadnjeg stupca prema prvom, uz ažuriranje matrice $B$.
# Nakon povratka u izvornu bazu QR dekompozicija ponovno daje gornji faktor.
#
# Sljedeći primjer ima jedan realni svojstveni korijen i kompleksno konjugirani
# par. Ulazi su realni, ali trenutačna implementacija računa u kompleksnoj aritmetici.

# %%
A = [-2.0 3 1; -3 -2 0.5; 0 0 -1]
B = [1.0 0; 0 1; 1 1]
U = lyapunov_factor(A, B)
X = U * U'
display(X)
display((relativni_rezidual=relative_residual(A, B, U),))

# %% [markdown]
# ## Usporedba s paketom
# Uspoređujem rješenja $X$, jer sam faktor može ovisiti o izboru faza.
# Ovo je provjera na jednom primjeru; za zaključke o točnosti treba više pokusa.

# %%
X_reference = MatrixEquations.lyapc(A, B * B')
display((relativna_greska=norm(X - X_reference) / norm(X_reference),))

# %% [markdown]
# ## Još treba napraviti
# - Razraditi realni Hammarlingov algoritam s blokovima veličine 2.
# - Dodati Sylvesterove sheme 1-solve i E-solve.
# - Isprobati loše uvjetovane probleme i usporediti vremena izvođenja.
# - Dopuniti objašnjenje izvoda i napisati zaključak.
# - Nakon toga povezati algoritam s primjerom mreže aerodroma.
#
# Zasad koristim male guste matrice. Nema skaliranja za vrlo velike ili male
# elemente, a shema s $R^2$ može gubiti točnost. Rezidual provjeravam nakon rješavanja.
#
# Literatura: D. C. Sorensen, Y. Zhou, *Direct methods for matrix Sylvester and
# Lyapunov equations*, Journal of Applied Mathematics 2003:6, 277–303.
