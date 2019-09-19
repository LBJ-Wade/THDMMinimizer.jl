@inline function _x(fields::Fields{T}) where T<:Real
  fields.r1^2 + fields.r2^2 + fields.i1^2 + fields.i2^2
end

@inline function _y(fields::Fields{T}) where T<:Real
  fields.c1^2 + fields.c2^2 + fields.c3^2 + fields.c4^2
end

@inline function _z(fields::Fields{T}) where T<:Real
  fields.c2 * fields.i1 + fields.c4 * fields.i2 + fields.c1 * fields.r1 +
  fields.c3 * fields.r2
end

@inline function _w(fields::Fields{T}) where T<:Real
  fields.c1 * fields.i1 + fields.c3 * fields.i2 - fields.c2 * fields.r1 -
  fields.c4 * fields.r2
end

"""
    gauge_sqaured_mass_matrix(fields::Fields, params::Params)

Returns the gauge squared-mass matrix given the fields values `fields` and
parameters `params`.
"""
function gauge_sqaured_mass_matrix(
  fields::Fields{T},
  params::Params{U}
) where {T <: Real, U<:Real}

  x = _x(fields)
  y = _y(fields)
  z = _z(fields)
  w = _w(fields)
  tw = params.gp / params.g

  params.g^2 / 4. *
  [
   x + y zero(T) zero(T) 2tw * z;
   zero(T) x + y zero(T) 2tw * w;
   zero(T) zero(T) x + y tw * (y - x);
   2tw * z 2tw * w tw * (y - x) tw^2 * (x + y)
  ]
end

"""
    gauge_sqaured_mass_matrix(fields::Fields, params::Params)

Return an array of the gauge squared masses in GeV.
"""
function gauge_squared_masses(
  fields::Fields{T},
  params::Params{U}
) where {T <: Real, U<:Real}

  x = _x(fields)
  y = _y(fields)
  z = _z(fields)
  w = _w(fields)
  tw = params.gp / params.g

  sw = params.gp / sqrt(params.gp^2 + params.g^2);
  cw = params.g / sqrt(params.gp^2 + params.g^2);

  sqrt_fac = sqrt((x + y)^2 + 16sw^2 * cw^2 * (w^2 - x * y + z^2));

  mW = params.g^2 / 4 * (x + y);
  mZ = params.g^2 / 8 * (1 + tw^2) * (x + y + sqrt_fac);
  mA = params.g^2 / 8 * (1 + tw^2) * (x + y - sqrt_fac);
  [mW, mW, mZ, mA]
end
