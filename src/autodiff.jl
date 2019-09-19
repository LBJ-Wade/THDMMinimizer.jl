"""
    autodiff_gradient(f, fields, params, fld_indices)

Compute the gradient matrix of the function `f` w.r.t. fields with indices
`fld_indices`. i.e. for each field-index in `fld_indices` we compute ∂f/∂ϕᵢ.

The acceptable values for `fld_indices` are:
    R1, R2, C1, C2, C3, C4, I1 or I2.

# Example
```julia
julia> flds = Fields([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
julia> pars = Params([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
julia> autodiff_gradient(potential_tree, flds, pars, (R1, R2))
[8.0, 8.0]
```
"""
function autodiff_gradient(
    f::Function,
    fields::Fields{T},
    params::Params{T},
    syms::Tuple{Vararg{FieldIndex}}
) where T<:Real
    _fields = to_dual(fields, syms...)
    _results = f(_fields, params)
    _poss = [Int(sym) for sym in syms]
    gradients = [partials(_res, _pos) for _pos in _poss, _res in _results]
    length(gradients) == 1 ? gradients[1] : gradients
end

"""
    autodiff_hessian(f, fields, params, fld_indices1, fld_indices2)

Compute the hessian matrix of the function `f` w.r.t. fields with indices
`fld_indices1` and `fld_indices2`. i.e. for each field-index in fld_indices1
and fld_indices2, we compute ∂²f/∂ϕᵢ∂ϕⱼ.

The acceptable values for `fld_indices1` and `fld_indices2` are:
    R1, R2, C1, C2, C3, C4, I1 or I2.

# Example
```julia
julia> f(fields, params) = potential_tree(fields, params);
julia> flds = Fields([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
julia> pars = Params([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
julia> hessian = autodiff_hessian(f, flds, pars, (R1, R2), (R1, R2))
[7.0 5.0; 5.0 7.0]
# We can check that this gives the same answer as the mass matrix:
julia> scalar_squared_mass_matrix(flds, pars)[1:2, 1:2]
[7.0 5.0; 5.0 7.0]
```
"""
function autodiff_hessian(
    f::Function,
    fields::Fields{T},
    params::Params{T},
    fld_indices1::Tuple{Vararg{FieldIndex}},
    fld_indices2::Tuple{Vararg{FieldIndex}}
) where T<:Real
    _fields = to_dual(to_dual(fields, fld_indices1...), fld_indices2...)
    _params = to_dual(to_dual(params))

    _pos1s = [Int(sym) for sym in fld_indices1]
    _pos2s = [Int(sym) for sym in fld_indices2]

    _results = f(_fields, _params)
    hessians = [[partials(partials(_res, _pos2), _pos1)
                 for _pos1 in _pos1s, _pos2 in _pos2s]
                for _res in _results]
    length(hessians) == 1 ? hessians[1] : hessians
end

"""
    autodiff_hessian(f, fields, params, fld_indices, par_indices)

Compute the hessian matrix of the function `f` w.r.t. fields and parameter
indices `fld_indices` and `par_indices`. i.e. for each field-index in
`fld_indices` and parameter index in `par_indices`, we compute ∂²f/∂ϕᵢ∂pⱼ.

The acceptable values for `fld_indices` are:
    R1, R2, C1, C2, C3, C4, I1 or I2.

The acceptable values for `par_indices` are:
    M112, M122, M222, Λ1, Λ2, Λ3, Λ4, Λ5.

# Example
```julia
julia> f(fields, params) = potential_tree(fields, params);
julia> flds = Fields([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
julia> pars = Params([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
julia> hessian = autodiff_hessian(f, flds, pars, (R1, R2), (M112, M122))
[1.0 -1.0; 0.0 -1.0]
```
"""
function autodiff_hessian(
    f::Function,
    fields::Fields{T},
    params::Params{T},
    sym1s::Tuple{Vararg{FieldIndex}},
    sym2s::Tuple{Vararg{ParamIndex}}
) where T<:Real
    _fields = to_dual(to_dual(fields, sym1s...))
    _params = to_dual(to_dual(params), sym2s...)

    _pos1s = [Int(sym) for sym in sym1s]
    _pos2s = [Int(sym) for sym in sym2s]

    _vals = f(_fields, _params)
    hessians = [[partials(partials(_val, _pos2), _pos1)
                 for _pos1 in _pos1s, _pos2 in _pos2s]
                for _val in _vals]
    length(hessians) == 1 ? hessians[1] : hessians
end
