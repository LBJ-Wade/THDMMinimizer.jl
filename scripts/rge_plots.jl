using THDMMinimizer
using DelimitedFiles
using LaTeXStrings
import PyPlot; const plt = PyPlot;
using Printf
using DifferentialEquations

data1 = readdlm(string(@__DIR__) * "/../data/verified_pos_mass_a1.csv", ',', skipstart=1);
data2 = readdlm(string(@__DIR__) * "/../data/verified_pos_mass_a2.csv", ',', skipstart=1);

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
smidxs = filter(i->abs(sqrt(data1[i,13]^2 + data1[i,14]^2) -246.0) < 1.0, i:length(data1[:,1]));
i = 10
point = data1[202,:]
begin
    μf=400.0
    point = data1[175,:]
    _, _, pars = row_to_vacs_and_params(point)
    μs, nvacs, cbvacs = run_vacuua(point; μf=μf)
    veffsn = [nvac.potential for nvac in nvacs]
    veffscb = [cbvac.potential for cbvac in cbvacs]
    plt.figure(dpi=100)
    plt.plot(μs, veffsn, "--", lw=2, label=L"$V_{\mathrm{EW}}$")
    plt.plot(μs, veffscb, lw=2, label=L"$V_{\mathrm{CB}}$")
    plt.xlabel(L"$\mu \ (\mathrm{GeV})$", fontsize=16)
    plt.legend(fontsize=14)
    plt.xlim([pars.μ, μf])
    plt.gcf()
    plt.savefig("figures/rgecb/cb_rge_" * string(175) * ".pdf")
    plt.close_figs()
end

nvac, cbvac, pars = row_to_vacs_and_params(data1[77,:]);
306, 214, 209, 208, 175

smidxsn = filter(i->abs(sqrt(data2[i,13]^2 + data2[i,14]^2) -246.0) < 1.0, i:length(data2[:,1]));
point = data2[smidxsn[1],:]

smidxsn[1]

plt.close_figs()
begin
    μf=400.0
    _, _, pars = row_to_vacs_and_params(point)
    μs, nvacs, cbvacs = run_vacuua(point; μf=μf)
    veffsn = [nvac.potential for nvac in nvacs]
    veffscb = [cbvac.potential for cbvac in cbvacs]
    plt.figure(dpi=100)
    plt.plot(μs, veffsn, lw=2, "--", label=L"$V_{\mathrm{EW}}$")
    plt.plot(μs, veffscb, lw=2, label=L"$V_{\mathrm{CB}}$")
    plt.xlabel(L"$\mu \ (\mathrm{GeV})$", fontsize=16)
    plt.legend(fontsize=14)
    plt.xlim([pars.μ, μf])
    plt.gcf()
    plt.savefig("normal_rge_" * string(smidxsn[1]) * ".pdf")
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
