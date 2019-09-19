"""
    generate_params_and_vacs(μ::Float64=HIGGS_VEV)

Generate parameters and vacuua (normal and CB) such that the parameters are
bounded and perturbative and the tree-level potential is approximately
extremized at both the normal and CB configurations.
"""
function generate_params_and_vacs(μ::Float64=HIGGS_VEV)
    done = false
    while !done
        nvac = generate_normal_vac(HIGGS_VEV)
        cbvac = generate_cb_vac(HIGGS_VEV)
        pars = Params(μ)
        try
            solve_tree_eqns!(nvac, cbvac, pars)
            are_perturbative = check_are_perturbative(pars)
            bounded = is_bounded(pars)
            if are_perturbative && bounded
                set_top_yukawa!(nvac, pars)
                return nvac, cbvac, pars
            end
        catch e
            if e isa JacobiMaxIterError
                continue
            else
                throw(e)
            end
        end
    end
end
