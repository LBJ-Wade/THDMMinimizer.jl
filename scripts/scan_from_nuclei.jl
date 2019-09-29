"""
Perturb and Reminimize:
-----------------------

Strategy 1: preserving v₁² + v₂² = v²
=====================================
1. Tweak the vacuaa by a small amount:
    - for normal, rotate using β → β + δβ
    - for CB, v₁ → v₁ + δv₁, v₂ → v₂ + δv₂, α → α + δα
2. Try to solve root equations.
3. If success, then find all other minima at new vevs
4. Determine if still a counter example and repeat

Strategy 2:
===========
1. Tweak parameters by a small amount:
2. Try to solve root equations.
3. If success, then find all other minima at new vevs
4. Determine if still a counter example and repeat
"""

using THDMMinimizer
using CSV
using ProgressMeter
using StatsBase

const fileidx = rand(1:10000); # choose random number for file extension.

function row_to_vacs_and_params(row; tol=1e-5)
    if length(row) == 18
        pars = Params([[float(row[i]) for i in 1:12]; 0.0])
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

struct MaxIter <: Exception end

function strategy_one(nvacuum::Vacuum, cbvacuum::Vacuum, params::Params; maxiter=100, tol=1e-5)
    shift() = (1 + 1e-2 * (2rand()-1))
    # generate new normal vacuum
    β = atan(nvacuum.v2 / nvacuum.v1) * shift() # shift by 5 -> -5%
    v = sqrt(nvacuum.v1^2 + nvacuum.v2^2)
    nvac = Vacuum(v * cos(β), v * sin(β), 0.0, 0.0, NotSet)
    # generate new CB vacuum
    cbvac = Vacuum(
        cbvacuum.v1 * shift(),
        cbvacuum.v2 * shift(),
        cbvacuum.α * shift(),
        0.0,
        NotSet
    )
    pars = deepcopy(params)
    iter = 0
    while iter < maxiter
        try
            try_find_root_eff!(nvac, cbvac, pars)
            good_ngrads = check_gradients(nvac, pars;tol=tol)
            good_cbgrads = check_gradients(cbvac, pars;tol=tol)
            perturbative = check_are_perturbative(pars)
            bounded = is_bounded(pars)
            if good_ngrads && good_cbgrads && perturbative && bounded
                good_nmasses = check_one_loop_masses(nvac, pars, :normal; tol=tol)
                good_cbmasses = check_one_loop_masses(cbvac, pars, :cb;tol=tol)
                if good_nmasses && good_cbmasses
                    categorize_vacuum!(nvac, pars;tol=tol)
                    categorize_vacuum!(cbvac, pars;tol=tol)
                    nvac.potential = potential_eff(Fields(nvac), pars)
                    cbvac.potential = potential_eff(Fields(cbvac), pars)
                    return nvac, cbvac, pars
                end
            end
        catch e
            if e isa JacobiMaxIterError
                iter += 1
                continue
            else
                throw(e)
            end
        end
        iter += 1
    end
    throw(MaxIter())
end

function find_new_ces(data_seed)
    fname_a1 = "typea1_$fileidx.csv"
    fname_a2 = "typea2_$fileidx.csv"

    seed = row_to_vacs_and_params(data_seed)
    for i in 1:20
        found = false
        while !found
            try
                nvac, cbvac, pars = strategy_one(seed[1], seed[2], seed[3]; maxiter=100)
                vacs = [nvac, cbvac]
                vacs = [vacs; find_new_minima(pars; tol=1e-5,method=Optim.BFGS())]
                type = catagorize_results(vacs)
                if type == :typea1
                    deepestn = find_deepest_normal_min(vacs)
                    deepestcb = find_deepest_cb_min(vacs)
                    write_results_to_file(fname_a1, deepestn, deepestcb, pars)
                    found = true
                    seed = (deepestn, deepestcb, pars)
                elseif type == :typea2
                    deepestn = find_deepest_normal_min(vacs)
                    deepestcb = find_deepest_cb_min(vacs)
                    write_results_to_file(fname_a2, deepestn, deepestcb, pars)
                    found = true
                    seed = (deepestn, deepestcb, pars)
                end
            catch e
                if e isa JacobiMaxIterError
                    continue
                elseif e isa MaxIter
                    continue
                else
                    throw(e)
                end
            end
        end
    end
end

function write_headers()
    fname_a1 = "typea1_$fileidx.csv"
    fname_a2 = "typea2_$fileidx.csv"

    header = "m112,m122,m222,"
    header *= "lam1,lam2,lam3,lam4,lam5,"
    header *= "mu,yt,gp,g,"
    header_a = header * "nvev1,nvev2,nvev3,cbvev1,cbvev2,cbvev3\n"

    open(fname_a1, "w") do io
        write(io, header_a)
    end
    open(fname_a2, "w") do io
        write(io, header_a)
    end
end

function scan(type::Symbol)
    if type == :a1
        data = CSV.read(string(@__DIR__) * "/../data/seeds_a1.csv");
    else
        data = CSV.read(string(@__DIR__) * "/../data/seeds_a2.csv");
    end

    write_headers()
    @showprogress 1 "Finding new CEs" for i in 1:length(data[:,1])
        for j in 1:5
            find_new_ces(data[i,:])
        end
    end
end
