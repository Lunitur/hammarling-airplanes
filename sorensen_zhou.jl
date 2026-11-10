using LinearAlgebra

function schur_blocks(R)
    blocks = UnitRange{Int}[]
    j = 1
    while j <= size(R, 1)
        k = j < size(R, 1) && !iszero(R[j + 1, j]) ? j + 1 : j
        push!(blocks, j:k)
        j = k + 1
    end
    blocks
end

function quasi_solve(M, C, blocks)
    X = zeros(promote_type(eltype(M), eltype(C)), size(C))
    for rows in Iterators.reverse(blocks)
        tail = last(rows)+1:size(M, 1)
        rhs = C[rows, :] - M[rows, tail] * X[tail, :]
        X[rows, :] = M[rows, rows] \ rhs
    end
    X
end

function pair_solve(R, Rsq, S, F, blocks)
    a, b, c, d = S[1, 1], S[1, 2], S[2, 1], S[2, 2]
    M = Rsq + (a + d) * R + (a * d - b * c) * I
    f1, f2 = F[:, 1], F[:, 2]
    C = hcat(R * f1 + d * f1 - c * f2, R * f2 + a * f2 - b * f1)
    quasi_solve(M, C, blocks)
end

# Realni Schurovi blokovi rješavaju se shemom 2-solve iz rada.
function sylvester_sz(A1::AbstractMatrix{<:Real}, A2::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real})
    m, n = size(A1, 1), size(A2, 1)
    m > 0 && n > 0 && size(A1) == (m, m) && size(A2) == (n, n) && size(D) == (m, n) ||
        throw(DimensionMismatch("Neusklađene dimenzije matrica"))
    all(isfinite, A1) && all(isfinite, A2) && all(isfinite, D) ||
        throw(ArgumentError("Elementi matrica moraju biti konačni"))

    S1, S2 = schur(Matrix{Float64}(A1)), schur(Matrix{Float64}(A2))
    R1, R2 = S1.T, S2.T
    C = S1.Z' * D * S2.Z
    Y = zeros(m, n)
    blocks = schur_blocks(R1)
    Rsq = R1 * R1
    for cols in schur_blocks(R2)
        prev = 1:first(cols)-1
        F = -C[:, cols] - Y[:, prev] * R2[prev, cols]
        Y[:, cols] = length(cols) == 1 ?
            quasi_solve(R1 + R2[first(cols), first(cols)] * I, F, blocks) :
            pair_solve(R1, Rsq, R2[cols, cols], F, blocks)
    end
    S1.Z * Y * S2.Z'
end

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

# QR s obrnutim stupcima daje gornji faktor za X = F*F', bez formiranja X.
function upper_factor(F)
    R = Matrix(qr(Matrix(F')[:, end:-1:1]).R)
    U = Matrix(reverse(R'; dims=(1, 2)))
    for j in axes(U, 2)
        d = U[j, j]
        iszero(d) || (U[:, j] *= conj(d) / abs(d))
    end
    UpperTriangular(U)
end

"""Računa gornji faktor U, X = U*U', za A*X + X*A' + B*B' = 0."""
function lyapunov_factor(A::AbstractMatrix, B::AbstractMatrix)
    n = size(A, 1)
    n > 0 && size(A) == (n, n) && size(B, 1) == n ||
        throw(DimensionMismatch("A mora biti kvadratna, a B imati isti broj redaka"))
    all(isfinite, A) && all(isfinite, B) ||
        throw(ArgumentError("Elementi matrica moraju biti konačni"))

    S = schur(Matrix{ComplexF64}(A))
    all(z -> real(z) < 0, S.values) || throw(ArgumentError("A mora biti stabilna"))
    C = S.Z' * Matrix{ComplexF64}(B)
    U = hammarling_complex_arithmetic(S.T, C)
    upper_factor(S.Z * U)
end

function relative_residual(A, B, U)
    X = U * U'
    D = B * B'
    denominator = 2 * norm(A) * norm(X) + norm(D)
    residual = norm(A * X + X * A' + D)
    iszero(denominator) ? residual : residual / denominator
end
