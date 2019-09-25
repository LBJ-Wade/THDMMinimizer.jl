function rge_system(du, u, p, μ)
    m112, m122, m222, λ1, λ2, λ3, λ4, λ5, yt, gp, g = u

    γ1 = 9 / 4 * g^2 + 3 / 4 * gp^2
    γ2 = 9 / 4 * g^2 + 3 / 4 * gp^2 - 3yt^2
    gs = 8.88577 / sqrt(21.612 + 7.0 * log(μ))
    au = -8gs^2 - 9 / 4 * g^2 - 17 / 12 * gp^2

    du[1] = 6λ1 * m112 + (4λ3 + 2λ4) * m222 - 2γ1 * m112
    du[2] = m122 * (2λ3 + 4λ4 + 6λ5 - γ1 - γ2)
    du[3] = 6λ2 * m222 + (4λ3 + 2λ4) * m112 - 2γ2 * m222

    du[4] = 12λ1^2 + 4λ3^2 + 4λ3 * λ4 + 2λ4^2 + 2λ5^2 + 9 / 4 * g^4 +
            3 / 2 * g^2 * gp^2 + 3 / 4 * gp^4 - 4γ1 * λ1
    du[5] = 12λ2^2 + 4λ3^2 + 4λ3 * λ4 + 2λ4^2 + 2λ5^2 + 9 / 4 * g^4 +
            3 / 2 * g^2 * gp^2 + 3 / 4 * gp^4 - 4γ2 * λ2 - 12yt^4
    du[6] = (λ1 + λ2) * (6λ3 + 2λ4) + 4λ3^2 + 2λ4^2 + 2λ5^2 + 9 / 4 * g^4 -
            3 / 2 * g^2 * gp^2 + 3 / 4 * gp^2 - 2 * (γ1 + γ2) * λ3
    du[7] = 2 * (λ1 + λ2) * λ4 + 8λ3 * λ4 + 4λ4^2 + 8λ5^2 - 2 * (γ1 + γ2) * λ4 +
            3g^2 * gp^2
    du[8] = 2λ5 * (λ1 + λ2 + 4λ3 + 6λ4 - γ1 - γ2)

    du[9] = au * yt + 9.0 / 2.0 * yt^3
    du[10] = 7gp^3
    du[11] = -3g^3

    du[:] = du[:] ./ (16π^2 * μ)
end

function run_parameters(params::Params{Float64}, μ1::Float64, μ2::Float64)
    u0 = [
        params.m112,
        params.m122,
        params.m222,
        params.λ1,
        params.λ2,
        params.λ3,
        params.λ4,
        params.λ5,
        params.yt,
        params.gp,
        params.g
    ]
    prob = ODEProblem(rge_system, u0, (μ1, μ2))
    sol = DifferentialEquations.solve(prob, Tsit5(), retol=1e-8, abstol=1e-8)
end
