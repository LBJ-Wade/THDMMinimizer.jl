struct SolveRootEqnsMaxIter <: Exception end
Base.showerror(io::IO, e::SolveRootEqnsMaxIter) = print(io, "Too many iterations in `solve_root_eqns`")

"""
    solve_root_eqns([;μ=HIGGS_VEV, maxiter=100)

Try to find parameters and vacuua (normal and CB) such that the effective
potential is extremized at the normal and CB vacuua. Additionally, we require
that the gradients are less than 1e-8, the potential is bounded and vacuua
contain the correct number of Goldstones.

We randomly choose vacuua and parameters such that the tree-level potential is
approximately extremized at the vacuaa. Then we try to find new parameters
that extremize the effective potential. We allow for `maxiter` number of tries.
"""
function solve_root_eqns(;μ::Float64 = HIGGS_VEV, maxiter::Int = 100, tol=1e-8)
    iter::Int = 1
    while iter < maxiter
        nvac, cbvac, pars = generate_params_and_vacs(μ)
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
    throw(SolveRootEqnsMaxIter())
end

"""
    find_new_minima(params::Params[;num_vacs=50, method=Optim.NewtonTrustRegion()])

Find new minima of the effective potential given the parameters `params` by
picking `num_vacs` new vacuua and performing a minimization. Additionally, we
perform checks on the gradients and masses to validate the new vacuua.
"""
function find_new_minima(params::Params{Float64}; num_vacs::Int=50, method=Optim.NewtonTrustRegion(),tol=1e-8)
    vacs = [generate_cb_vac(params.μ) for _ in 1:num_vacs]
    goodvacs = Array{Vacuum}(undef, 0)
    for vac in vacs
        minimize!(vac, params; method=method)
        categorize_vacuum!(vac, params; tol=tol)
        # perform checks
        good_grad = check_gradients(vac, params;tol=tol)
        type = abs(vac.α) < 1e-5 ? :normal : :cb
        good_masses = check_one_loop_masses(vac, params, type; tol=tol)
        if good_grad && good_masses && vac.extremum_type == Minimum
            push!(goodvacs, vac)
        end
    end
    goodvacs
end

"""
    catagorize_results(vacs::Vector{Vacuum})

Given a set of vacuua, determine if the set is of type: A1, A2, B or C. These
types are such that:
    A1 : normal + CB minima with CB the deepest => returns :typea1
    A2 : normal + CB minima with normal the deepest => returns :typea2
    B : only a normal minimum => returns :typeb
    C : only a CB minimum => returns :typec
If no minima exists, returns :undefined.
"""
function catagorize_results(vacs::Vector{Vacuum})
    minima = filter(vac->vac.extremum_type == Minimum, vacs)
    if length(minima) > 0
        cb_mins = filter(vac->abs(vac.α) > 1.0, minima)
        n_mins = filter(vac->abs(vac.α) < 1e-5, minima)
        if length(cb_mins) > 0 && length(n_mins) > 0 # Type A1,A2
            deepestcb = first(cb_mins)
            for vac in cb_mins
                if vac.potential < deepestcb.potential
                    deepestcb = vac
                end
            end
            deepestn = first(n_mins)
            for vac in n_mins
                if vac.potential < deepestn.potential
                    deepestn = vac
                end
            end
            if deepestn.potential < deepestcb.potential
                return :typea2
            else
                return :typea1
            end
        elseif length(n_mins) > 0 && length(cb_mins) == 0 # Type B
            return :typeb
        elseif length(cb_mins) > 0 && length(n_mins) == 0 # Type C
            return :typec
        end
    else # do something
        return :undefined
    end
end

"""
    find_deepest_normal_min(vacs::Vector{Vacuum})

Find the deepest normal minimum out of a set of vacuua.
"""
function find_deepest_normal_min(vacs::Vector{Vacuum})
    normal_vacs = filter(vac->abs(vac.α) < 1e-5, vacs)
    normal_mins = filter(vac->vac.extremum_type == Minimum, normal_vacs)
    deepest_n = first(normal_mins)
    for normal_min in normal_mins
        if normal_min.potential < deepest_n.potential
            deepest_n = normal_min
        end
    end
    deepest_n
end

"""
    find_deepest_cb_min(vacs::Vector{Vacuum})

Find the deepest cb minimum out of a set of vacuua.
"""
function find_deepest_cb_min(vacs::Vector{Vacuum})
    cb_vacs = filter(vac->abs(vac.α) > 1.0, vacs)
    cb_mins = filter(vac->vac.extremum_type == Minimum, cb_vacs)
    deepest_cb = first(cb_mins)
    for cb_min in cb_mins
        if cb_min.potential < deepest_cb.potential
            deepest_cb = cb_min
        end
    end
    deepest_cb
end

"""
    write_results_to_file(filename, nvac, cbvac, params)

Write the parameters and vacuua (normal and CB) to the file `filename`.
"""
function write_results_to_file(filename, nvac, cbvac, params)
    res = [params.m112 params.m122 params.m222 params.λ1 params.λ2 params.λ3 params.λ4 params.λ5 params.μ params.yt params.gp params.g nvac.v1 nvac.v2 nvac.α cbvac.v1 cbvac.v2 cbvac.α]
    open(filename, "a") do io
        writedlm(io, res, ',')
    end
end

"""
    write_results_to_file(filename, nvac, cbvac, params)

Write the parameters and the single vacuum to the file `filename`.
"""
function write_results_to_file(filename, vac, params)
    res = [params.m112 params.m122 params.m222 params.λ1 params.λ2 params.λ3 params.λ4 params.λ5 params.μ params.yt params.gp params.g vac.v1 vac.v2 vac.α]
    open(filename, "a") do io
        writedlm(io, res, ',')
    end
end
