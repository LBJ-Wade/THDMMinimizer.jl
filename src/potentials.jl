@inline function _x1(fields::Fields{T}) where T <: Real
  (fields.c1^2 + fields.c2^2 + fields.i1^2 + fields.r1^2) / 2
end;

@inline function _x2(fields::Fields{T}) where T <: Real
  (fields.c3^2 + fields.c4^2 + fields.i2^2 + fields.r2^2) / 2
end;

@inline function _x3(fields::Fields{T}) where T <: Real
  (fields.c1 * fields.c3 + fields.c2 * fields.c4 + fields.i1 * fields.i2 +
   fields.r1 * fields.r2) / 2
end;

@inline function _x4(fields::Fields{T}) where T <: Real
  (-(fields.c2 * fields.c3) + fields.c1 * fields.c4 + fields.i2 * fields.r1 -
   fields.i1 * fields.r2) / 2
end;


"""
    potential_tree(fields::Fields, params::Params)

Computes the THDM potential given field and parameters values
`fields` and `params`.
"""
function potential_tree(
  fields::Fields{T},
  params::Params{U}
) where {T<:Real,U<:Real}
  r1 = fields.r1
  r2 = fields.r2
  c1 = fields.c1
  c2 = fields.c2
  c3 = fields.c3
  c4 = fields.c4
  i1 = fields.i1
  i2 = fields.i2

  m112 = params.m112
  m122 = params.m122
  m222 = params.m222
  λ1 = params.λ1
  λ2 = params.λ2
  λ3 = params.λ3
  λ4 = params.λ4
  λ5 = params.λ5

  return @horner(
    r1,
    @horner(
      r2,
      @horner(
        c1,
        @horner(
          c2,
          @horner(
            c3,
            @horner(
              c4,
              @horner(
                i1,
                @horner(i2, 0, 0, m222 / 2., 0, λ2 / 8.),
                @horner(i2, 0, -m122),
                @horner(i2, m112 / 2., 0, (2 * λ3 + 2 * λ4 + 2 * λ5) / 8.),
                0,
                λ1 / 8.
              ),
              0,
              @horner(i1, @horner(i2, m222 / 2., 0, λ2 / 4.), 0, λ3 / 4.),
              0,
              λ2 / 8.
            ),
            0,
            @horner(
              c4,
              @horner(i1, @horner(i2, m222 / 2., 0, λ2 / 4.), 0, λ3 / 4.),
              0,
              λ2 / 4.
            ),
            0,
            λ2 / 8.
          ),
          @horner(
            c4,
            0,
            @horner(i1, -m122, @horner(i2, 0, (4 * λ4 + 4 * λ5) / 8.))
          ),
          @horner(
            c3,
            @horner(
              c4,
              @horner(i1, @horner(i2, m112 / 2., 0, λ3 / 4.), 0, λ1 / 4.),
              0,
              (2 * λ3 + 2 * λ4 + 2 * λ5) / 8.
            ),
            0,
            (2 * λ3 + 2 * λ4 - 2 * λ5) / 8.
          ),
          0,
          λ1 / 8.
        ),
        @horner(
          c2,
          @horner(
            c3,
            0,
            @horner(i1, -m122, @horner(i2, 0, (4 * λ4 + 4 * λ5) / 8.))
          ),
          @horner(c3, 0, @horner(c4, 0, λ5))
        ),
        @horner(
          c2,
          @horner(
            c3,
            @horner(
              c4,
              @horner(i1, @horner(i2, m112 / 2., 0, λ3 / 4.), 0, λ1 / 4.),
              0,
              (2 * λ3 + 2 * λ4 - 2 * λ5) / 8.
            ),
            0,
            (2 * λ3 + 2 * λ4 + 2 * λ5) / 8.
          ),
          0,
          λ1 / 4.
        ),
        0,
        λ1 / 8.
      ),
      @horner(
        c1,
        @horner(c2, 0, @horner(c3, 0, @horner(i1, 0, (4 * λ4 - 4 * λ5) / 8.))),
        @horner(c4, 0, @horner(i1, 0, (-4 * λ4 + 4 * λ5) / 8.))
      ),
      @horner(
        c1,
        @horner(
          c2,
          @horner(
            c3,
            @horner(
              c4,
              @horner(
                i1,
                @horner(i2, m222 / 2., 0, λ2 / 4.),
                0,
                (2 * λ3 + 2 * λ4 - 2 * λ5) / 8.
              ),
              0,
              λ2 / 4.
            ),
            0,
            λ2 / 4.
          ),
          0,
          λ3 / 4.
        ),
        0,
        λ3 / 4.
      ),
      0,
      λ2 / 8.
    ),
    @horner(
      r2,
      @horner(
        c1,
        @horner(c2, 0, @horner(c3, 0, @horner(i2, 0, (-4 * λ4 + 4 * λ5) / 8.))),
        @horner(c4, 0, @horner(i2, 0, (4 * λ4 - 4 * λ5) / 8.))
      ),
      @horner(
        c1,
        @horner(
          c2,
          @horner(i1, -m122, @horner(i2, 0, λ5)),
          @horner(c4, 0, (4 * λ4 + 4 * λ5) / 8.)
        ),
        @horner(c3, 0, (4 * λ4 + 4 * λ5) / 8.)
      )
    ),
    @horner(
      r2,
      @horner(
        c1,
        @horner(
          c2,
          @horner(
            c3,
            @horner(
              c4,
              @horner(
                i1,
                @horner(i2, m112 / 2., 0, (2 * λ3 + 2 * λ4 - 2 * λ5) / 8.),
                0,
                λ1 / 4.
              ),
              0,
              λ3 / 4.
            ),
            0,
            λ3 / 4.
          ),
          0,
          λ1 / 4.
        ),
        0,
        λ1 / 4.
      ),
      0,
      (2 * λ3 + 2 * λ4 + 2 * λ5) / 8.
    ),
    0,
    λ1 / 8.
  )
end

"""
    potential_one_loop_generic(msqrds, c, μ, dof)

Returns the generic form of the one-loop effective potential. This form is
  ∑ᵢ dof * mᵢ⁴/64π² * (log(mᵢ²/μ²) - c)
where the sum is over all the masses mᵢ, `μ` is the renormalization scale,
`c` is 3/2 (5/6) for scalars (fermions) and `dof` is the total number of
degrees of freedom in the fields.
"""
function potential_one_loop_generic(msqrds, c, μ, dof)
  loop = zero(eltype(msqrds))
  for msqrd in msqrds
    loop += abs(msqrd) > 0 ? msqrd^2 * (0.5log(msqrd^2 / μ^4) - c) : 0
  end
  dof * loop / (64π^2)
end;

"""
    potential_one_loop(fields::Fields, params::Params)

Computes the one-loop correction to the THDM scalar potential.
"""
function potential_one_loop(
  fields::Fields{T},
  params::Params{U}
) where {T<:Real,U<:Real}

  scalar_msqrds = scalar_squared_masses(fields, params)
  gauge_msqrds = gauge_squared_masses(fields, params)
  ghost_msqrds = params.ξ .* gauge_msqrds
  fermion_msqrds = fermion_squared_masses(fields, params)

  return (potential_one_loop_generic(scalar_msqrds, 3 / 2, params.μ, 1) +
          # 3 for polarizations
          potential_one_loop_generic(gauge_msqrds, 5 / 6, params.μ, 3) +
          # Ghost contributions
          potential_one_loop_generic(ghost_msqrds, 3 / 2, params.μ, -1) +
          # minus for fermion and 12 for 2 spins, 2 charges and 3 colors
          potential_one_loop_generic(fermion_msqrds, 3 / 2, params.μ, -12))
end;

"""
    potential_eff(fields::Array{T,1}, params::Array{U,1})

Compute the effective potential evaluated at the field values `fields` with
parameters `params`.

# Example
julia> nvac = Vacuum(227.11, 94.56, 0.0, 0.0, NotSet)
julia> params = Params([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
julia> set_top_yukawa!(nvac, params)
julia> potential_eff(Fields(nvac), params)
6.776513836866474e8
"""
function potential_eff(
  fields::Fields{T},
  params::Params{U}
) where {T<:Real,U<:Real}
  potential_tree(fields, params) + potential_one_loop(fields, params)
end;

"""
  one_loop_masses(vacuum, params)

Compute the one-loop scalar masses.

# Example
julia> vac = Vacuum(1.0, 1.0, 1.0, 0.0, NotSet)
julia> params = Params([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
julia> one_loop_masses(vac, params)
[1.23, 1.95, 0.52, 1.27, 1.48, 0.52, 0.52, 1.23]
"""
function one_loop_masses(vac::Vacuum, params::Params{Float64})
  hessian = autodiff_hessian(
    potential_eff,
    Fields(vac),
    params,
    (R1, R2, C1, C2, C3, C4, I1, I2),
    (R1, R2, C1, C2, C3, C4, I1, I2)
  )
  jacobi(hessian)
end;
