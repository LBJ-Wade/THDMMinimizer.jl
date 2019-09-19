using THDMMinimizer
using LinearAlgebra
using ForwardDiff: Dual, partials, value
using SpecialFunctions
using Test

include("test_utils.jl")

@testset "THDMMinimizer.jl" begin
    @testset "test jacobi" begin
        M1 = test_matrix1(rand(), rand())
        M2 = test_matrix2(rand(), rand())
        M3 = test_matrix3(rand(), rand())
        Ms = [M1, M2, M3]

        for M in Ms
            evals = sort(eigvals(M))
            evals_jacobi = sort(jacobi(M))
            for (eval, evalj) in zip(evals, evals_jacobi)
                @test abs(eval - evalj) < 1e-10
            end
        end
    end

    @testset "test jacobi derivs" begin
        # Test matrix 1
        a1, b1 = rand(), rand()
        devals1aj = sort(jacobi(test_matrix1(
            Dual{Float64}(a1, 1),
            Dual{Float64}(b1, 0)
        )))
        devals1bj = sort(jacobi(test_matrix1(
            Dual{Float64}(a1, 0),
            Dual{Float64}(b1, 1)
        )))
        devals1a = sort(test_matrix1_eval_derivs(a1, b1, :a))
        devals1b = sort(test_matrix1_eval_derivs(a1, b1, :b))
        @test abs(partials(devals1aj[1], 1) - partials(devals1a[1], 1)) < 1e-10
        @test abs(partials(devals1aj[2], 1) - partials(devals1a[2], 1)) < 1e-10
        @test abs(partials(devals1bj[1], 1) - partials(devals1b[1], 1)) < 1e-10
        @test abs(partials(devals1bj[2], 1) - partials(devals1b[2], 1)) < 1e-10
        # Test matrix 2
        a2, b2 = rand(), rand()
        devals2aj = sort(jacobi(test_matrix2(
            Dual{Float64}(a2, 1),
            Dual{Float64}(b2, 0)
        )))
        devals2bj = sort(jacobi(test_matrix2(
            Dual{Float64}(a2, 0),
            Dual{Float64}(b2, 1)
        )))
        devals2a = sort(test_matrix2_eval_derivs(a2, b2, :a))
        devals2b = sort(test_matrix2_eval_derivs(a2, b2, :b))
        @test abs(partials(devals2aj[1], 1) - partials(devals2a[1], 1)) < 1e-10
        @test abs(partials(devals2aj[2], 1) - partials(devals2a[2], 1)) < 1e-10
        @test abs(partials(devals2bj[1], 1) - partials(devals2b[1], 1)) < 1e-10
        @test abs(partials(devals2bj[2], 1) - partials(devals2b[2], 1)) < 1e-10
        # Test matrix 3
        a3, b3 = rand(), rand()
        devals3aj = sort(jacobi(test_matrix3(
            Dual{Float64}(a3, 1),
            Dual{Float64}(b3, 0)
        )))
        devals3bj = sort(jacobi(test_matrix3(
            Dual{Float64}(a3, 0),
            Dual{Float64}(b3, 1)
        )))
        devals3a = sort(test_matrix3_eval_derivs(a3, b3, :a))
        devals3b = sort(test_matrix3_eval_derivs(a3, b3, :b))
        @test abs(partials(devals3aj[1], 1) - partials(devals3a[1], 1)) < 1e-10
        @test abs(partials(devals3aj[2], 1) - partials(devals3a[2], 1)) < 1e-10
        @test abs(partials(devals3bj[1], 1) - partials(devals3b[1], 1)) < 1e-10
        @test abs(partials(devals3bj[2], 1) - partials(devals3b[2], 1)) < 1e-10
    end
end
