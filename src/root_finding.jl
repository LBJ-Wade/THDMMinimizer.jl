"""
  try_find_root_eff!(nvac::Vacuum, cbvac::Vacuum, params::Params{Float64})

Tries to find a set of parameters such the the effective potential is
extremized at both the normal and CB vacuua `nvac` and `cbvac`. The parameters
are modified in-place.
"""
function try_find_root_eff!(
  nvac::Vacuum,
  cbvac::Vacuum,
  params::Params{Float64}
)
  # Choose 5 parameters at random to solve the root equations
  allpars = [M112, M122, M222, Λ1, Λ2, Λ3, Λ4, Λ5]
  parsample = Tuple(sample(allpars, 5; replace = false))

  flds_n = Fields(nvac)
  flds_cb = Fields(cbvac)

  # gradiants of the effective potential evaluated at normal and cb vacuua
  function f!(F, x)
    for (i, par) in enumerate(parsample)
      params[par] = x[i]
    end
    derivs_n = autodiff_gradient(potential_eff, flds_n, params, (R1, R2))
    derivs_cb = autodiff_gradient(potential_eff, flds_cb, params, (R1, R2, C1))
    F[1] = derivs_n[1]
    F[2] = derivs_n[2]
    F[3] = derivs_cb[1]
    F[4] = derivs_cb[2]
    F[5] = derivs_cb[3]
  end
  # jacobian of the effective potential evaluated at normal and cb vacuua
  function j!(J, x)
    for (i, par) in enumerate(parsample)
      params[par] = x[i]
    end
    hessian_n = autodiff_hessian(
      potential_eff,
      flds_n,
      params,
      (R1, R2),
      parsample
    )
    hessian_cb = autodiff_hessian(
      potential_eff,
      flds_cb,
      params,
      (R1, R2, C1),
      parsample
    )
    for i = 1:5
      J[1, i] = hessian_n[1, i]
      J[2, i] = hessian_n[2, i]
      J[3, i] = hessian_cb[1, i]
      J[4, i] = hessian_cb[2, i]
      J[5, i] = hessian_cb[3, i]
    end
  end

  nlsolve(f!, j!, params[parsample]; xtol = 0.0, ftol = 1e-8)
end

"""
    solve_tree_eqns!(nvac::Vacuum, cbvac::Vacuum, params::Params)

Finds parameters such that the tree-level potential is extremized at the normal
and charge-breaking vacuua `nvac` and `cbvac`. This is done via least squares.
i.e., if system is singluar, we find the best-fit parameters.
"""
function solve_tree_eqns!(nvac::Vacuum, cbvac::Vacuum, params::Params{Float64})
  nv1, nv2 = nvac.v1, nvac.v2
  cv1, cv2, cv3 = cbvac.v1, cbvac.v2, cbvac.α

  coefs1 = [
    nv1,
    -nv2,
    0,
    nv1^3 / 2.,
    0,
    (nv1 * nv2^2) / 2.,
    (nv1 * nv2^2) / 2.,
    (nv1 * nv2^2) / 2.
  ]
  coefs2 = [
    0,
    -nv1,
    nv2,
    0,
    nv2^3 / 2.,
    (nv1^2 * nv2) / 2.,
    (nv1^2 * nv2) / 2.,
    (nv1^2 * nv2) / 2.
  ]
  coefs3 = [
    cv1,
    -cv2,
    0,
    (4 * cv1^3 + 4 * cv1 * cv3^2) / 8.,
    0,
    (cv1 * cv2^2) / 2.,
    (cv1 * cv2^2) / 2.,
    (cv1 * cv2^2) / 2.
  ]
  coefs4 = [
    0,
    -cv1,
    cv2,
    0,
    cv2^3 / 2.,
    (4 * cv1^2 * cv2 + 4 * cv2 * cv3^2) / 8.,
    (cv1^2 * cv2) / 2.,
    (cv1^2 * cv2) / 2.
  ]
  coefs5 = [
    cv3,
    0,
    0,
    (4 * cv1^2 * cv3 + 4 * cv3^3) / 8.,
    0,
    (cv2^2 * cv3) / 2.,
    0,
    0
  ]

  paridxs = [M112, M122, M222, Λ1, Λ2, Λ3, Λ4, Λ5]
  parsample = sample(paridxs, 5; replace = false)
  parunused = [par for par in paridxs if !(par in parsample)]

  lhs1 = [coefs1[Int(par)] for par in parsample]
  lhs2 = [coefs2[Int(par)] for par in parsample]
  lhs3 = [coefs3[Int(par)] for par in parsample]
  lhs4 = [coefs4[Int(par)] for par in parsample]
  lhs5 = [coefs5[Int(par)] for par in parsample]

  rhs1 = -sum([coefs1[Int(par)] * params[par] for par in parunused])
  rhs2 = -sum([coefs2[Int(par)] * params[par] for par in parunused])
  rhs3 = -sum([coefs3[Int(par)] * params[par] for par in parunused])
  rhs4 = -sum([coefs4[Int(par)] * params[par] for par in parunused])
  rhs5 = -sum([coefs5[Int(par)] * params[par] for par in parunused])

  M = transpose(hcat(lhs1, lhs2, lhs3, lhs4, lhs5))
  b = [rhs1, rhs2, rhs3, rhs4, rhs5]
  sols = pinv(M' * M) * M' * b
  for (par, sol) in zip(parsample, sols)
    params[par] = sol
  end
end

"""
  find_tree_vacuua(params::Params{Float64})

Find all of the vacuua of the tree-level potential using homotopy continuation.
"""
function find_tree_vacuua(params::Params{Float64}; tol = 1e-8)
  m112, m122, m222 = params.m112, params.m122, params.m222
  λ1, λ2, λ3, λ4, λ5 = params.λ1, params.λ2, params.λ3, params.λ4, params.λ5

  @polyvar v1 v2 v3

  result = HomotopyContinuation.solve([
    ((2m112 * v1 + v2 * (-2m122 + (λ3 + λ4 + λ5) * v1 * v2) +
      λ1 * v1 * (v1^2 + v3^2)) / 2),
    (-(m122 * v1) +
     (v2 * (2m222 + (λ3 + λ4 + λ5) * v1^2 + λ2 * v2^2 + λ3 * v3^2)) / 2),
    ((v3 * (2m112 + λ3 * v2^2 + λ1 * (v1^2 + v3^2))) / 2)
  ])

  sols = real_solutions(result)
  vacuua = Array{Vacuum,1}(undef, 0)
  for sol in sols
    vac = Vacuum(sol[1], sol[2], sol[3], 0.0, NotSet)
    fldidxs = (R1, R2, C1, C2, C3, C4, I1, I2)
    grad = autodiff_gradient(potential_tree, Fields(vac), params, fldidxs)
    if maximum(abs.(grad)) < tol
      push!(vacuua, vac)
    end
  end
  vacuua
end
