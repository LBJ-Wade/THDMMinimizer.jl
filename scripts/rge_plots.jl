using THDMMinimizer
using DelimitedFiles
using LaTeXStrings
import PyPlot; const plt = PyPlot;
using Printf
using DifferentialEquations

data1 = readdlm(string(@__DIR__) * "/../data/verified_a1.csv", ',', skipstart=1);
data2 = readdlm(string(@__DIR__) * "/../data/verified_a2.csv", ',', skipstart=1);

function row_to_vacs_and_params(row; tol=1e-5)
    if length(row) == 18
        pars = Params([[row[i] for i in 1:12]; 0.0])
        nvac = Vacuum(row[13], row[14], row[15], 0.0, NotSet)
        cbvac = Vacuum(row[16], row[17], row[18], 0.0, NotSet)
        categorize_vacuum!(nvac, pars; tol=tol)
        categorize_vacuum!(cbvac, pars; tol=tol)
        nvac.potential = potential_eff(Fields(nvac), pars)
        cbvac.potential = potential_eff(Fields(cbvac), pars)
        return nvac, cbvac, pars
    elseif length(row) == 15
        pars = Params([row[i] for i in 1:12])
        vac = Vacuum(row[13], row[14], row[15], 0.0, NotSet)
        categorize_vacuum!(vac, pars; tol=tol)
        vac.potential = potential_eff(Fields(vac), pars)
        return vac, pars
    else
        throw("row doesn't contain correct number of columns")
    end
end

function run_vacuua(row; μf::Float64=400.0, tol=1e-5)
    nvac, cbvac, params = row_to_vacs_and_params(row)

    μs = range(params.μ, stop=μf, length=100);
    sols = run_parameters(params, μs[1], μs[end]);
    nvacs = [nvac]
    cbvacs = [cbvac]

    for μ in μs[2:end]
        ps = sols(μ);
        ps = [[ps[1:8]; μ; ps[9:end]]; 0.0]
        pars = Params(ps)
        nvac_new = deepcopy(nvacs[end])
        cbvac_new = deepcopy(cbvacs[end])
        minimize!(nvac_new, pars; method=Optim.NewtonTrustRegion())
        minimize!(cbvac_new, pars; method=Optim.NewtonTrustRegion())
        categorize_vacuum!(nvac_new, pars; tol=tol)
        categorize_vacuum!(cbvac_new, pars; tol=tol)
        nvacs = [nvacs; nvac_new]
        cbvacs = [cbvacs; cbvac_new]
    end
    μs, nvacs, cbvacs
end

function num_to_sn(num)
    strs = split(@sprintf("%.2E", num), "E")
    if strs[2] == "+00"
        return strs[1]
    elseif '+' in strs[2]
        strexp = strip(split(strs[2], "+")[2], '0')
        return "$(strs[1])\\times10^{$strexp}"
    elseif '-' in strs[2]
        strexp = strip(split(strs[2], "-")[2], '0')
        return "$(strs[1])\\times10^{-$strexp}"
    end
end

function create_title_string(pars)
    m112 = num_to_sn(pars.m112)
    m122 = num_to_sn(pars.m122)
    m222 = num_to_sn(pars.m222)
    λ1 = num_to_sn(pars.λ1)
    λ2 = num_to_sn(pars.λ2)
    λ3 = num_to_sn(pars.λ3)
    λ4 = num_to_sn(pars.λ4)
    λ5 = num_to_sn(pars.λ5)
    strms = LaTeXString("m_{11}^2 =" * m112 * L", m_{12}^2 = " * m122 *
             ", m_{12}^2 = " * m122)
    strλs = (L"\lambda_{1} = " * λ1 * L", \lambda_{2} = " * λ2 *
             L", \lambda_{3} = " * λ3 * L", \lambda_{4} = " * λ4 *
             L", \lambda_{5} = " * λ5)
    strms * strλs
end

# good points:
# 16, 22, 26, 31,34, 51, 60, 71, 72, 83, 86, 91, 927, 924, 919, 917, 909, 900, 897
point = 1890#good2pts2[15]
plt.close_figs()
begin
    μf=200.0
    _, _, pars = row_to_vacs_and_params(data1[point,:])
    μs, nvacs, cbvacs = run_vacuua(data1[point,:]; μf=μf)
    veffsn = [nvac.potential for nvac in nvacs]
    veffscb = [cbvac.potential for cbvac in cbvacs]
    plt.figure(dpi=100)
    plt.plot(μs, veffsn, lw=2, label=L"$V_{\mathrm{EW}}$")
    plt.plot(μs, veffscb, lw=2, label=L"$V_{\mathrm{CB}}$")
    plt.xlabel(L"$\mu \ (\mathrm{GeV})$", fontsize=16)
    plt.legend(fontsize=14)
    plt.xlim([pars.μ, μf])
    #plt.text(250, -8.7e8, latexstring("m_{11}^2 = " * num_to_sn(pars.m112)), fontsize=12)
    #plt.text(250, -8.76e8, latexstring("m_{12}^2 = " * num_to_sn(pars.m122)), fontsize=12)
    #plt.text(250, -8.82e8, latexstring("m_{22}^2 = " * num_to_sn(pars.m222)), fontsize=12)
    #plt.text(253.5, -8.88e8, latexstring("\\lambda_{1} = " * num_to_sn(pars.λ1)), fontsize=12)
    #plt.text(253.5, -8.94e8, latexstring("\\lambda_{2} = " * num_to_sn(pars.λ2)), fontsize=12)
    #plt.text(253.5, -9.0e8, latexstring("\\lambda_{3} = " * num_to_sn(pars.λ3)), fontsize=12)
    #plt.text(253.5, -9.06e8, latexstring("\\lambda_{4} = " * num_to_sn(pars.λ4)), fontsize=12)
    #plt.text(253.5, -9.12e8, latexstring("\\lambda_{5} = " * num_to_sn(pars.λ5)), fontsize=12)
    plt.gcf()
    #plt.savefig("cb_rge_" * string(point) * ".pdf")
end

# good points: 456 = good2pts[46], 533 = good2pts[52]
# 537 = good2pts[53], 767 = good2pts[76], 113 = good2pts[1042]
point = good2pts[113]
plt.close_figs()
begin
    μf=400.0
    _, _, pars = row_to_vacs_and_params(data2[point,:])
    μs, nvacs, cbvacs = run_vacuua(data2[point,:]; μf=μf)
    veffsn = [nvac.potential for nvac in nvacs]
    veffscb = [cbvac.potential for cbvac in cbvacs]
    plt.figure(dpi=100)
    plt.plot(μs, veffsn, lw=2, label=L"$V_{\mathrm{EW}}$")
    plt.plot(μs, veffscb, lw=2, label=L"$V_{\mathrm{CB}}$")
    plt.xlabel(L"$\mu \ (\mathrm{GeV})$", fontsize=16)
    plt.legend(fontsize=14)
    plt.xlim([pars.μ, μf])
    #plt.text(250, -6.34e8, latexstring("m_{11}^2 = " * num_to_sn(pars.m112)), fontsize=12)
    #plt.text(250, -6.4e8, latexstring("m_{12}^2 = " * num_to_sn(pars.m122)), fontsize=12)
    #plt.text(250, -6.46e8, latexstring("m_{22}^2 = " * num_to_sn(pars.m222)), fontsize=12)
    #plt.text(253.5, -6.52e8, latexstring("\\lambda_{1} = " * num_to_sn(pars.λ1)), fontsize=12)
    #plt.text(253.5, -6.58e8, latexstring("\\lambda_{2} = " * num_to_sn(pars.λ2)), fontsize=12)
    #plt.text(253.5, -6.64e8, latexstring("\\lambda_{3} = " * num_to_sn(pars.λ3)), fontsize=12)
    #plt.text(253.5, -6.70e8, latexstring("\\lambda_{4} = " * num_to_sn(pars.λ4)), fontsize=12)
    #plt.text(253.5, -6.76e8, latexstring("\\lambda_{5} = " * num_to_sn(pars.λ5)), fontsize=12)
    plt.gcf()
    #plt.savefig("normal_rge_" * string(point) * ".pdf")
end

function find_pts_with_sm_vevs(data)
    sm_vevs = Array{Int,1}(undef,0)
    for i in 1:length(data[:, 1])
        nvac,_,pars=row_to_vacs_and_params(data[i, :])
        val = 246.0 - sqrt(nvac.v1^2 + nvac.v2^2)
        if abs(val) < 1e-2
            sm_vevs = [sm_vevs; Int(i)]
        end
    end
    sm_vevs
end

good2pts = find_pts_with_sm_vevs(data2)

good2pts2 = find_pts_with_sm_vevs(data1)


struct Dual{T<:Real} <:Real
    val::T # real component of the dual number.
    eps::T # infinitesimal component of the dual number. eps^2 = 0
end

function Base.:*(z::Dual{T}, w::Dual{T}) where T<:Real
    Dual{T}(z.val * w.val, z.val * w.eps + z.eps * w.val)
end


function Base.sin(z::Dual{T}) where T<:Real
    Dual{T}(sin(z.val), cos(z.eps))
end

function Base.cos(z::Dual{T}) where T<:Real
    Dual{T}(cos(z.val), -sin(z.eps))
end

d1 = Dual{Float64}(1.0, 1.0)
d2 = Dual{Float64}(2.0, 0.0)

d1 * d2

@show sin(d1)

sin(1)
cos(1)

d3 = Dual{Dual{Float64}}(Dual{Float64}(1.0, 0.0), Dual{Float64}(0.0, 1.0))

@show sin(d3)
