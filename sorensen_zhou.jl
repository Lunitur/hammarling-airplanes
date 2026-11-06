# %% [markdown]
# # Sorensen–Zhou: Sylvesterove i Lyapunovljeve jednadžbe
# Implementacija rada *Direct methods for matrix Sylvester and Lyapunov equations*
# (2003), priloženog kao `Sorensen-Zhou-Lyap.pdf`, algoritmi 2.1, 3.1 i 3.2.
# Ćelije se pokreću redom. Koristimo Float64 i ComplexF64 te guste matrice.
# Konvencija je **X = U Uᴴ**, s gornjetrokutastim faktorom U.

# %%
using LinearAlgebra
using Random
using Test
using BenchmarkTools
import MatrixEquations

# %% [markdown]
# ## 1. Blokovska supstitucija
# Realni Schurov oblik ima dijagonalne blokove veličine 1 ili 2. Njihove granice
# prepoznajemo iz strukturnih nula, bez predznačnog testa iz algoritma 2.1.
# Supstitucija košta O(n²) po desnoj strani. Granice blokova zadržavamo
# i za polinom u R; ne radimo gustu LU faktorizaciju cijelog sustava.

# %%
function schur_blocks(R)
    blocks = UnitRange{Int}[]
    j, n = 1, size(R, 1)
    while j <= n
        k = j < n && !iszero(R[j + 1, j]) ? j + 1 : j
        push!(blocks, j:k)
        j = k + 1
    end
    blocks
end

function quasi_solve(M, C, blocks)
    X = zeros(promote_type(eltype(M), eltype(C)), size(C))
    n = size(M, 1)
    for rows in Iterators.reverse(blocks)
        tail = last(rows)+1:n
        rhs = C[rows, :] - M[rows, tail] * X[tail, :]
        if length(rows) == 1
            pivot = M[first(rows), first(rows)]
            iszero(pivot) && throw(SingularException(first(rows)))
            X[rows, :] = rhs / pivot
        else
            X[rows, :] = M[rows, rows] \ rhs
        end
    end
    X
end

# %% [markdown]
# ## 2. Sylvester: 1-solve, 2-solve i E-solve
# Rješavamo A₁X + XA₂ + D = 0. Schurove transformacije daju
# R₁Y + YR₂ + Q₁ᴴDQ₂ = 0. Stupce računamo slijeva nadesno.
# Za blok S = [a b; c d] shema 2-solve koristi (2.16)–(2.17):
# (R₁² + (a+d)R₁ + (ad-bc)I)[y₁ y₂]
# = [R₁f₁ + df₁ − cf₂, R₁f₂ + af₂ − bf₁].
# R₁² računamo jednom. 1-solve koristi (2.10); E-solve koristi realnu
# bazu kompleksnog svojstvenog vektora iz (2.19)–(2.26).

# %%
function pair_solve(M, Msq, S, F, blocks; scheme=:two)
    a, b, c, d = S[1, 1], S[1, 2], S[2, 1], S[2, 2]
    K = Msq + (a + d) * M + (a * d - b * c) * I
    f1, f2 = F[:, 1], F[:, 2]
    if scheme == :one
        x = quasi_solve(K, reshape(M * f1 + d * f1 - c * f2, :, 1), blocks)
        return hcat(x, (f1 - (M + a * I) * vec(x)) / c)
    elseif scheme == :two
        G = hcat(M * f1 + d * f1 - c * f2, M * f2 + a * f2 - b * f1)
        return quasi_solve(K, G, blocks)
    elseif scheme == :eigen
        E = eigen(S)
        k = findfirst(z -> imag(z) > 0, E.values)
        isnothing(k) && throw(ArgumentError("E-solve zahtijeva kompleksni par"))
        μ, η = real(E.values[k]), imag(E.values[k])
        v = E.vectors[:, k]
        V = hcat(real(v), imag(v))
        G = F * V
        H = (M + μ * I) * G + η * hcat(G[:, 2], -G[:, 1])
        return quasi_solve(K, H, blocks) / V
    end
    throw(ArgumentError("Nepoznata shema: $scheme"))
end

function sylvester_sz(A1::AbstractMatrix{<:Real}, A2::AbstractMatrix{<:Real}, D;
                      scheme=:two)
    scheme in (:one, :two, :eigen) || throw(ArgumentError("Nepoznata shema"))
    m, n = size(A1, 1), size(A2, 1)
    m > 0 && n > 0 && size(A1) == (m, m) && size(A2) == (n, n) && size(D) == (m, n) ||
        throw(DimensionMismatch("Neusklađene ili prazne dimenzije"))
    all(isfinite, A1) && all(isfinite, A2) && all(isfinite, D) ||
        throw(ArgumentError("Ulazi moraju biti konačni"))
    S1, S2 = schur(Matrix{Float64}(A1)), schur(Matrix{Float64}(A2))
    R1, R2 = S1.T, S2.T
    C = S1.Z' * D * S2.Z
    Y = zeros(eltype(C), m, n)
    blocks = schur_blocks(R1)
    Rsq = R1 * R1
    for cols in schur_blocks(R2)
        prev = 1:first(cols)-1
        F = -C[:, cols] - Y[:, prev] * R2[prev, cols]
        Y[:, cols] = length(cols) == 1 ?
            quasi_solve(R1 + R2[first(cols), first(cols)] * I, F, blocks) :
            pair_solve(R1, Rsq, R2[cols, cols], F, blocks; scheme)
    end
    S1.Z * Y * S2.Z'
end

# %% [markdown]
# ## 3. Kompleksni Hammarling — algoritam 3.1
# Nakon kompleksne Schurove dekompozicije rješavamo RP + PRᴴ + BBᴴ = 0.
# Za zadnji red b postavljamo μ = ‖b‖₂, α = √(−2 Re Rⱼⱼ), τ = μ/α.
# Jednadžba (3.26) daje u, a (3.19) ažurira B rangom 1. Normalizacija
# reda b izbjegava izravno kvadriranje njegove norme. B se kopira.

# %%
function hammarling_complex_arithmetic(R, B)
    n = size(R, 1)
    C = Matrix{ComplexF64}(B)
    U = zeros(ComplexF64, n, n)
    for j in n:-1:1
        μ = norm(C[j, :])
        α = sqrt(-2 * real(R[j, j]))
        U[j, j] = μ / α
        if μ > 0 && j > 1
            b = C[j:j, :] / μ
            rhs = α * C[1:j-1, :] * b' + R[1:j-1, j:j] * (μ / α)
            u = UpperTriangular(R[1:j-1, 1:j-1] + conj(R[j, j]) * I) \ (-rhs)
            U[1:j-1, j:j] = u
            C[1:j-1, :] -= α * u * b
        end
    end
    UpperTriangular(U)
end

# %% [markdown]
# ## 4. Realni Hammarling — algoritam 3.2, shema 2-solve
# Zadnji stupac P dobivamo iz (3.31), a zadnja dva iz 2-solve sustava.
# Iz P₂₂ = T Tᵀ izvlačimo gornji faktor odozdo prema gore:
# τ₂₂ = √p₂₂, τ₁₂ = p₁₂/τ₂₂, τ₁₁ = √(p₁₁−τ₁₂²).
# Nakon izdvajanja zadnjeg stupca oduzimamo njegov doprinos od prethodnog
# stupca P. Taj korak nije eksplicitno naveden u pseudokodu na str. 297.
# Slučaj j=1 pripada skalarnoj grani. Nulti pivot dopuštamo samo uz numerički
# nulti ostatak stupca; znatno negativan pivot prijavljujemo.

# %%
function extract_column!(U, C, j)
    p = U[j, j]
    tol = 100 * eps(Float64) * maximum(abs, U[1:j, j]; init=0.0)
    p < -tol && throw(DomainError(p, "Negativan pivot: gubitak numeričke točnosti"))
    if p <= 0
        maximum(abs, U[1:j-1, j]; init=0.0) <= tol ||
            throw(DomainError(p, "Nulti pivot uz nenulti ostatak stupca"))
        U[1:j, j] .= 0
        return
    end
    U[j, j] = sqrt(p)
    U[1:j-1, j] /= U[j, j]
    C[j, :] /= U[j, j]
    C[1:j-1, :] -= U[1:j-1, j:j] * C[j:j, :]
end

function hammarling_real_arithmetic(R, B)
    n = size(R, 1)
    C = Matrix{Float64}(B)
    U = zeros(n, n)
    Rsq = R * R
    blocks = schur_blocks(R)
    for k in length(blocks):-1:1
        cols = blocks[k]
        j = last(cols)
        M = R[1:j, 1:j]
        F = -C[1:j, :] * C[cols, :]'
        if length(cols) == 1
            U[1:j, cols] = quasi_solve(M + R[j, j] * I, F, blocks[1:k])
            extract_column!(U, C, j)
        else
            U[1:j, cols] = pair_solve(M, Rsq[1:j, 1:j], R[cols, cols]', F, blocks[1:k])
            U[j-1, j] = (U[j-1, j] + U[j, j-1]) / 2
            U[j, j-1] = 0
            extract_column!(U, C, j)
            U[1:j-1, j-1] -= U[1:j-1, j] * U[j-1, j]
            extract_column!(U, C, j-1)
        end
    end
    UpperTriangular(U)
end

# %% [markdown]
# ## 5. Povratak u izvornu bazu
# X = (QU)(QU)ᴴ, ali QU općenito nije trokutast. RQ dekompoziciju provodimo
# QR dekompozicijom (QU)ᴴ s obrnutim redoslijedom stupaca. Time vraćamo gornji
# faktor bez formiranja X. Fazama stupaca dobivamo nenegativnu dijagonalu.
# Javno sučelje provjerava stabilnost (Re λ < 0), dimenzije i konačnost ulaza.
# Funkcije `hammarling_*_arithmetic` su jezgre za stabilni Schurov oblik.

# %%
function upper_factor(F)
    R = Matrix(qr(Matrix(F')[:, end:-1:1]).R)
    U = Matrix(reverse(R'; dims=(1, 2)))
    for j in axes(U, 2)
        d = U[j, j]
        iszero(d) || (U[:, j] *= conj(d) / abs(d))
    end
    UpperTriangular(U)
end

function lyapunov_factor(A::AbstractMatrix, B::AbstractMatrix; arithmetic=:real)
    n = size(A, 1)
    n > 0 && size(A) == (n, n) && size(B, 1) == n ||
        throw(DimensionMismatch("A mora biti neprazna kvadratna, a B imati n redaka"))
    arithmetic in (:real, :complex) || throw(ArgumentError("Aritmetika: :real ili :complex"))
    arithmetic == :real && !(eltype(A) <: Real && eltype(B) <: Real) &&
        throw(ArgumentError("Za kompleksne ulaze odaberite arithmetic=:complex"))
    all(isfinite, A) && all(isfinite, B) || throw(ArgumentError("Ulazi moraju biti konačni"))
    T = arithmetic == :real ? Float64 : ComplexF64
    S = schur(Matrix{T}(A))
    all(z -> real(z) < 0, S.values) || throw(ArgumentError("A nije Hurwitz stabilna"))
    C = S.Z' * Matrix{T}(B)
    U = arithmetic == :real ? hammarling_real_arithmetic(S.T, C) :
                             hammarling_complex_arithmetic(S.T, C)
    upper_factor(S.Z * U)
end

function relative_residual(A, B, U)
    X = U * U'
    D = B * B'
    denominator = 2 * norm(A) * norm(X) + norm(D)
    r = norm(A * X + X * A' + D)
    iszero(denominator) ? r : r / denominator
end

# %% [markdown]
# ## 6. Reproducibilan primjer
# A ima realni svojstveni korijen i kompleksno konjugirani par. Referentni
# `MatrixEquations.plyapc` koristi istu konvenciju X = U Uᴴ.

# %%
A = [-2.0 3 1; -3 -2 0.5; 0 0 -1]
B = [1.0 0; 0 1; 1 1]
Ur = lyapunov_factor(A, B)
Uc = lyapunov_factor(A, B; arithmetic=:complex)
Uref = MatrixEquations.plyapc(A, B)
example_results = (
    real_residual=relative_residual(A, B, Ur),
    complex_residual=relative_residual(A, B, Uc),
    relative_solution_error=norm(Ur * Ur' - Uref * Uref') / norm(Uref * Uref'),
)
display(example_results)

# %% [markdown]
# ## 7. Provjere
# Provjeravamo oba tipa Schurovih blokova, skalarni završetak, p<n i p>n,
# kompleksne koeficijente, singularno rješenje i neispravne ulaze.
# Poznato rješenje Sylvesterove jednadžbe neovisno provjerava sve tri sheme.

# %%
@testset "Sorensen–Zhou" begin
    rng = MersenneTwister(2026)
    @testset "Lyapunov real/complex, tall/wide B" begin
        for n in (1, 2, 5, 8), p in (1, n, n + 3)
            H = randn(rng, n, n)
            Atest = H - (maximum(real, eigvals(H)) + 1) * I
            Btest = randn(rng, n, p)
            Acopy, Bcopy = copy(Atest), copy(Btest)
            Xref = MatrixEquations.lyapc(Atest, Btest * Btest')
            for arithmetic in (:real, :complex)
                U = lyapunov_factor(Atest, Btest; arithmetic)
                @test istriu(U)
                @test minimum(real, diag(U)) >= 0
                @test relative_residual(Atest, Btest, U) < 1e-11
                @test U * U' ≈ Xref rtol=1e-9
            end
            @test Atest == Acopy && Btest == Bcopy
        end
    end
    @testset "Semidefinitni i kompleksni slučajevi" begin
        for Atest in (Matrix(Diagonal([-1., -2, -3])), [-1. 2 0; -2 -1 0; 0 0 -3])
            for Btest in (zeros(3, 0), zeros(3, 2), reshape([0., 0, 1], 3, 1))
                for arithmetic in (:real, :complex)
                    U = lyapunov_factor(Atest, Btest; arithmetic)
                    @test relative_residual(Atest, Btest, U) < 1e-12
                end
            end
        end
        H = randn(rng, ComplexF64, 5, 5)
        Ac = H - (maximum(real, eigvals(H)) + 1) * I
        Bc = randn(rng, ComplexF64, 5, 2)
        U = lyapunov_factor(Ac, Bc; arithmetic=:complex)
        @test relative_residual(Ac, Bc, U) < 1e-12
        @test U * U' ≈ MatrixEquations.lyapc(Ac, Bc * Bc') rtol=1e-10
    end
    @testset "Sylvesterove sheme" begin
        A1 = [-1. 4 1; -2 -1 0; 0 0 -3]
        A2 = [-2. 3; -1 -2]
        Xtrue = randn(rng, 3, 2)
        D = -(A1 * Xtrue + Xtrue * A2)
        for scheme in (:one, :two, :eigen)
            @test sylvester_sz(A1, A2, D; scheme) ≈ Xtrue rtol=1e-12
            A3 = Matrix(Diagonal([-2., -4]))
            @test sylvester_sz(A1, A3, -(A1 * Xtrue + Xtrue * A3); scheme) ≈ Xtrue
        end
    end
    @testset "Ulazne provjere" begin
        @test_throws ArgumentError lyapunov_factor([1.;;], [1.;;])
        @test_throws ArgumentError lyapunov_factor([0.;;], [1.;;])
        @test_throws ArgumentError lyapunov_factor([NaN;;], [1.;;])
        @test_throws ArgumentError lyapunov_factor([-1.;;], [1.;;]; arithmetic=:other)
        @test_throws DimensionMismatch lyapunov_factor(zeros(2, 3), zeros(2, 1))
        @test_throws ArgumentError lyapunov_factor([-1. + im;;], [1.;;])
    end
end

# %% [markdown]
# ## 8. Mjerenje vremena
# Mjerenje uključuje Schurovu dekompoziciju i povratak u izvornu bazu.
# Uspoređujemo izvedbu na istom problemu; ne pretpostavljamo prednost u brzini
# iz rezultata za MATLAB iz 2003. Paketni solver koristi optimizirane rutine.
# `run_benchmarks()` se poziva ručno kako pokretanje svih ćelija ostane kratko.

# %%
function run_benchmarks(; n=100, p=5)
    rng = MersenneTwister(42)
    H = randn(rng, n, n)
    A = H - (maximum(real, eigvals(H)) + 1) * I
    B = randn(rng, n, p)
    (
        real=@benchmark(lyapunov_factor($A, $B), samples=10, evals=1),
        complex=@benchmark(lyapunov_factor($A, $B; arithmetic=:complex), samples=10, evals=1),
        reference=@benchmark(MatrixEquations.plyapc($A, $B), samples=10, evals=1),
    )
end

# %% [markdown]
# ## Literatura i ograničenja
# D. C. Sorensen, Y. Zhou (2003), *Direct methods for matrix Sylvester and
# Lyapunov equations*, Journal of Applied Mathematics 2003:6, 277–303.
# Izvor: `Sorensen-Zhou-Lyap.pdf`, algoritmi 2.1, 3.1 i 3.2.
# Sylvesterov problem zahtijeva da A₁ i −A₂ nemaju zajedničkih svojstvenih
# vrijednosti. Bliski spektri i nenormalne matrice mogu dati loše uvjetovan problem.
# Ovo je obrazovna implementacija za male i srednje guste probleme. Nema
# automatskog skaliranja protiv overflowa; polinom R² može gubiti točnost.
# Rezidual uvijek provjerite. Za ekstremne skale koristite paketni solver.
