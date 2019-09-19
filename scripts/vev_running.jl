using THDMMinimizer
using CSV
using DataFrames
using Plots

data1 = CSV.read(string(@__DIR__) * "/../data/typea1.csv")
data2 = CSV.read(string(@__DIR__) * "/../data/typea2.csv")

function run_vacuua(datanum, type; μf::Float64=400.0)
    global data2
    global data1
    datapt = (type == :normal ? data2[datanum,:] : data1[datanum,:])

    params = Params([float(datapt[i]) for i in 1:12]);
    nvac = Vacuum(datapt[13], datapt[14], datapt[15], 0.0, NotSet);
    cbvac = Vacuum(datapt[16], datapt[17], datapt[18], 0.0, NotSet);
    nvac.potential = potential_eff(Fields(nvac), params);
    cbvac.potential = potential_eff(Fields(cbvac), params);
    categorize_vacuum!(nvac, one_loop_masses(nvac, params))
    categorize_vacuum!(cbvac, one_loop_masses(cbvac, params))

    μs = range(params.μ, stop=μf, length=100);
    sols = run_parameters(params, μs[1], μs[end]);
    convergedn = minimize!(nvac, params)
    convergedcb = minimize!(cbvac, params)
    nvacs = [nvac]
    cbvacs = [cbvac]
    for μ in μs[2:end]
        ps = sols(μ);
        ps = [ps[1:8]; μ; ps[9:end]]
        pars = Params(ps)
        _nvac = deepcopy(nvacs[end])
        _cbvac = deepcopy(cbvacs[end])
        convergedn = minimize!(_nvac, pars)
        convergedcb = minimize!(_cbvac, pars)

        categorize_vacuum!(_nvac, one_loop_masses(_nvac, pars))
        categorize_vacuum!(_cbvac,  one_loop_masses(_cbvac, pars))
        nvacs = [nvacs; _nvac]
        cbvacs = [cbvacs; _cbvac]
    end
    μs, nvacs, cbvacs
end

μs, nvacs, cbvacs = run_vacuua(3, :cb; μf=100.0)

veffsn = [_nvac.potential for _nvac in nvacs]
veffscb = [_cbvac.potential for _cbvac in cbvacs]

plot(μs, veffsn)
plot!(μs, veffscb)
