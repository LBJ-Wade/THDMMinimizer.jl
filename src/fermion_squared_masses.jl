"""
    fermion_squared_masses(fields::Fields, params::Params)

Return an array of all the squared fermion masses.
"""
@inline function fermion_squared_masses(
    fields::Fields{T},
    params::Params{U}
) where {T <: Real, U<:Real}
    [(fields.c3^2 + fields.c4^2 + fields.i2^2 + fields.r2^2) * params.yt^2 / 2]
end
