"""
    generate_params_and_vacs(μ::Float64=HIGGS_VEV)

Generate parameters and vacuua (normal and CB) such that the parameters are
bounded and perturbative and the tree-level potential is
extremized at both the normal and CB configurations.
"""
function generate_params_and_vacs(μ::Float64=HIGGS_VEV; tol=1e-8)
    done = false
    while !done
        nvac = generate_normal_vac(HIGGS_VEV)
        cbvac = generate_cb_vac(HIGGS_VEV)
        pars = Params(μ)
        try
            solve_tree_eqns!(nvac, cbvac, pars)
            are_perturbative = check_are_perturbative(pars)
            bounded = is_bounded(pars)
            ngrad = autodiff_gradient(potential_tree, Fields(nvac), pars, (R1, R2, C1, C2, C3, C4, I1, I2))
            cbgrad = autodiff_gradient(potential_tree, Fields(cbvac), pars, (R1, R2, C1, C2, C3, C4, I1, I2))
            maxngrad = maximum(abs.(ngrad))
            maxcbgrad = maximum(abs.(cbgrad))
            if are_perturbative && bounded && maxngrad < tol && maxcbgrad < tol
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
