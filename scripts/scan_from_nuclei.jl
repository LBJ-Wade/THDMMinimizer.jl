using THDMMinimizer
using CSV

data1 = CSV.read(string(@__DIR__) * "/../data/verified_a1.csv", ',', skipstart=1);
data2 = CSV.read(string(@__DIR__) * "/../data/verified_a2.csv", ',', skipstart=1);

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

"""
Perturb and Reminimize:
-----------------------

1. Tweak the vacuaa by a small amount:
    - for normal, rotate using β → β + δβ
    - for CB, v₁ → v₁ + δv₁, v₂ → v₂ + δv₂, α → α + δα
2. Try to solve root equations.
3. If success, then find all other minima

"""
