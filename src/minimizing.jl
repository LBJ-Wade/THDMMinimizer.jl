"""
    minimize!(vac::Vacuum, params::Params{Float64})

Minimize the effective potential starting at the vacuum `vac` given the
parameters `params`. The vacuum will be modified in-place to reflect the
solutions and the method will return a boolean specifying if the minimization
was successful.
"""
function minimize!(vac::Vacuum, params::Params{Float64}; method=Optim.NewtonTrustRegion())


    tup = (R1, R2, C1)
    f(x) = potential_eff(Fields(x), params)
    function g!(G, x)
        G[:] = autodiff_gradient(potential_eff, Fields(x), params, tup)
    end
    function h!(H, x)
        H[:,:] = autodiff_hessian(potential_eff, Fields(x), params, tup, tup)
    end

    res = Optim.optimize(
        f, g!, h!,
        [vac.v1, vac.v2, vac.α],
        method,
        Optim.Options(g_tol=1e-8, f_tol=-1.0, x_tol=-1.0)
    )
    vac.v1, vac.v2, vac.α = Optim.minimizer(res)
    vac.potential = potential_eff(Fields(vac), params)
    res
end
