using THDMMinimizer
using CSV
using DataFrames

df_a1 = DataFrame(CSV.read(string(@__DIR__) * "/../data/candidates_a1.csv"))
df_a2 = DataFrame(CSV.read(string(@__DIR__) * "/../data/candidates_a2.csv"))
df_b = DataFrame(CSV.read(string(@__DIR__) * "/../data/candidates_b.csv"))
df_c = DataFrame(CSV.read(string(@__DIR__) * "/../data/candidates_c.csv"))


function row_to_vacs_and_params(row; tol=1e-5)
    if length(row) == 18
        pars = Params([row[i] for i in 1:12])
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

function verify_data(row; tol=1e-5)
    if length(row) == 18
        nvac, cbvac, pars = row_to_vacs_and_params(row; tol=tol)
        goodnvac = check_gradients(nvac, pars; tol=tol)
        goodnvac = goodnvac && check_one_loop_masses(nvac, pars, :normal; tol=tol)
        goodcbvac = check_gradients(cbvac, pars; tol=tol)
        goodcbvac = goodcbvac && check_one_loop_masses(cbvac, pars, :cb; tol=tol)
        perturbative = check_are_perturbative(pars)
        bounded = is_bounded(pars)
        return goodnvac && goodcbvac && perturbative && bounded
    elseif length(row) == 15
        vac, pars = row_to_vacs_and_params(row; tol=tol)
        goodvac = check_gradients(vac, pars; tol=tol)
        goodvac = goodvac && check_one_loop_masses(vac, pars, :normal; tol=tol)
        goodpars = check_are_perturbative(pars)
        goodpars = goodpars && is_bounded(pars)
        return goodvac && goodpars
    else
        throw("row doesn't contain correct number of columns")
    end
end


df_a1_verified = filter(verify_data, df_a1)
df_a2_verified = filter(verify_data, df_a2)
df_b_verified = filter(verify_data, df_b)
df_c_verified = filter(verify_data, df_c)

CSV.write(string(@__DIR__) * "/../data/verified_a1.csv", df_a1_verified)
CSV.write(string(@__DIR__) * "/../data/verified_a2.csv", df_a2_verified)
CSV.write(string(@__DIR__) * "/../data/verified_b.csv", df_b_verified)
CSV.write(string(@__DIR__) * "/../data/verified_c.csv", df_c_verified)
