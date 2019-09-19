"""
    check_gradients(vacuum, params [,tol=1e-8])

Check that the gradients of the effective potential given a vacuum
configuration `vacuum` and parameters `params` are zero within tolerance `tol`.
"""
function check_gradients(vac::Vacuum, params::Params{Float64}; tol=1e-8)
    flds = Fields(vac)
    fields_idxs = (R1, R2, C1, C2, C3, C4, I1, I2)
    derivs = autodiff_gradient(potential_eff, flds, params, fields_idxs)
    maximum(abs.(derivs)) < tol
end

"""
    check_one_loop_masses(masses, type[,tol=1e-8])

Check that the correct number of goldstone bosons are present in the array
`masses`. If `type` is `:cb`, there there should be 4 goldstones. If the type is
normal, there should be 3 golstones. A mass is consider zero if it is less than
`tol` in absolute value.
"""
function check_one_loop_masses(masses::Array{Float64,1}, type::Symbol; tol=1e-8)
    ms = sort(abs.(masses))
    good::Bool = true
    good = good && ms[1] < tol
    good = good && ms[2] < tol
    good = good && ms[3] < tol
    good = good && (type == :cb ? ms[4] < tol : true)
    good
end

"""
    check_one_loop_masses(vac::Vacuum, params::Params, type[,tol=1e-8])

Check that the correct number of goldstone bosons are present given the vacuum
`vac` and parameters `params`. If `type` is `:cb`, there there should be 4
goldstones. If the type is normal, there should be 3 golstones. A mass is
consider zero if it is less than `tol` in absolute value.
"""
function check_one_loop_masses(vac::Vacuum, params::Params, type::Symbol; tol=1e-8)
    masses = one_loop_masses(vac, params)
    check_one_loop_masses(masses, type; tol=tol)
end

"""
    check_are_perturbative(params)

Check that the parameters are perturbative. We consider parameters to be
perturbative if: the masses are less than 10 times the renormalization scale
(to avoid large logs) and if the dimensionless parameters are much less than 4π
(we take the cutoff to be 10).
"""
function check_are_perturbative(params::Params{Float64})
    good = true
    good = good && sqrt(abs(params.m112)) < (10params.μ)^2
    good = good && sqrt(abs(params.m122)) < (10params.μ)^2
    good = good && sqrt(abs(params.m222)) < (10params.μ)^2
    good = good && abs(params.λ1) < 10.0
    good = good && abs(params.λ2) < 10.0
    good = good && abs(params.λ3) < 10.0
    good = good && abs(params.λ4) < 10.0
    good = good && abs(params.λ5) < 10.0
    good
end

"""
    categorize_vacuum!(vac::Vacuum, masses::Array{Float64,1}; tol=1e-8)

Categorize the vacuum `vac` as a `Minimum`, `Maximum`, `Saddle` or `Undefined`
given the one-loop masses `masses`. The catagorization is as follows:
    - all m >= 0                 => minimum
    - all m <= 0                 => maximum
    - some m <=0 and some m >= 0 => saddle
    - all m == 0                 => undefined
"""
function categorize_vacuum!(vac::Vacuum, masses::Array{Float64,1}; tol=1e-8)
    n_neg_masses = count(m -> (abs(m) > tol && m < 0.0), masses)
    n_pos_masses = count(m -> (abs(m) > tol && m > 0.0), masses)
    if n_neg_masses == 0 && n_pos_masses > 0 # minimum
        vac.extremum_type = Minimum
    elseif n_neg_masses > 0 && n_pos_masses == 0 # maximum
        vac.extremum_type = Maximum
    elseif n_neg_masses > 0 && n_pos_masses > 0 # saddle
        vac.extremum_type = Saddle
    else
        vac.extremum_type = Undefined
    end
end

"""
    categorize_vacuum!(vac::Vacuum, params::Params; tol=1e-8)

Categorize the vacuum `vac` as a `Minimum`, `Maximum`, `Saddle` or `Undefined`
given the parameters `params`. Same as `categorize_vacuum!(vac, masses)` but
the masses are computed.
"""
function categorize_vacuum!(vac::Vacuum, params::Params{Float64}; tol=1e-8)
    masses = one_loop_masses(vac, params)
    n_neg_masses = count(m -> (abs(m) > tol && m < 0.0), masses)
    n_pos_masses = count(m -> (abs(m) > tol && m > 0.0), masses)
    if n_neg_masses == 0 && n_pos_masses > 0 # minimum
        vac.extremum_type = Minimum
    elseif n_neg_masses > 0 && n_pos_masses == 0 # maximum
        vac.extremum_type = Maximum
    elseif n_neg_masses > 0 && n_pos_masses > 0 # saddle
        vac.extremum_type = Saddle
    else
        vac.extremum_type = Undefined
    end
end
