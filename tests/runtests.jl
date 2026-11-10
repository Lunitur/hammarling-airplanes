using Test
using LinearAlgebra
using Random
import MatrixEquations
include(joinpath(@__DIR__, "..", "sorensen_zhou.jl"))

@testset "Sylvesterova jednadžba" begin
    A1 = [-1.0 4 1; -2 -1 0; 0 0 -3]
    Y = [1.0 2; 0 -1; 3 1]
    for A2 in ([-2.0 3; -1 -2], Matrix(Diagonal([-2.0, -4])))
        D = -(A1 * Y + Y * A2)
        @test sylvester_sz(A1, A2, D) ≈ Y rtol=1e-12
    end
end

@testset "Kompleksni Hammarlingov postupak" begin
    @test lyapunov_factor(reshape([-2.0], 1, 1), reshape([4.0], 1, 1))[1, 1] ≈ 2
    rng = MersenneTwister(42)
    for T in (Float64, ComplexF64), n in (2, 5), p in (1, n + 1)
        H = randn(rng, T, n, n)
        A = H - (maximum(real, eigvals(H)) + 1) * I
        B = randn(rng, T, n, p)
        A_before, B_before = copy(A), copy(B)
        U = lyapunov_factor(A, B)
        X_reference = MatrixEquations.lyapc(A, B * B')
        @test istriu(U)
        @test minimum(real, diag(U)) >= 0
        @test relative_residual(A, B, U) < 1e-12
        @test U * U' ≈ X_reference rtol=1e-10
        @test A == A_before && B == B_before
    end

    A = Matrix(Diagonal([-1.0, -2, -3]))
    for B in (zeros(3, 1), reshape([0.0, 0, 1], 3, 1))
        U = lyapunov_factor(A, B)
        @test U * U' ≈ MatrixEquations.lyapc(A, B * B')
        @test relative_residual(A, B, U) < 1e-12
    end
    @test_throws ArgumentError lyapunov_factor(reshape([1.0], 1, 1), ones(1, 1))
    @test_throws DimensionMismatch lyapunov_factor(zeros(2, 3), ones(2, 1))
end
