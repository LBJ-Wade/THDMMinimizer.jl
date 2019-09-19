"""
Indices for the THDM fields. i.e. the field `r1` is always in location of
`Int(R1)`, etc.
"""
@enum FieldIndex begin
    R1 = 1
    R2 = 2
    C1 = 3
    C2 = 4
    C3 = 5
    C4 = 6
    I1 = 7
    I2 = 8
end

"Fields struct for the 8 THDM fields"
mutable struct Fields{T<:Real}
    "real component of lower component of first doublet: Re(ϕ₁[2])"
    r1::T
    "real component of lower component of second doublet: Re(ϕ₂[2])"
    r2::T
    "real component of upper component of first doublet: Re(ϕ₁[1])"
    c1::T
    "imaginary component of upper component of first doublet: Im(ϕ₁[1])"
    c2::T
    "real component of upper component of second doublet: Re(ϕ₂[1])"
    c3::T
    "imaginary component of upper component of second doublet: Im(ϕ₂[1])"
    c4::T
    "imaginary component of lower component of first doublet: Im(ϕ₁[2])"
    i1::T
    "imaginary component of lower component of second doublet: Im(ϕ₂[2])"
    i2::T
end

"""
    Fields(arr::Vector)

Create a fields struct from a vector of fields. The vector can be of length 3
(to specify r1, r2 and c1) or of length 8 (to specify all fields.)
"""
function Fields(arr::Vector{T}) where T <: Real
    if length(arr) == 8
        return Fields(arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8])
    elseif length(arr) == 3
        Fields(arr[1], arr[2], arr[3], zero(T), zero(T), zero(T), zero(T), zero(T))
    else
        throw("Fields can only be initialized with array of length 3 or 8.")
    end
end

"""
    Fields(arr::Vector)

Create a fields struct from the values for `r1`, `r2` and `c1`. All other
fields are set to zero.
"""
function Fields(r1::T, r2::T, c1::T) where T <: Real
    Fields(r1, r2, c1, zero(T), zero(T), zero(T), zero(T), zero(T))
end

function Base.getindex(fields::Fields{T}, fld::FieldIndex) where T<:Real
    if fld == R1
        return fields.r1
    elseif fld == R2
        return fields.r2
    elseif fld == C1
        return fields.c1
    elseif fld == C2
        return fields.c2
    elseif fld == C3
        return fields.c3
    elseif fld == C4
        return fields.c4
    elseif fld == I1
        return fields.i1
    else
        return fields.i2
    end
end

function Base.getindex(fields::Fields{T}, flds::Array{FieldIndex,1}) where T<:Real
    [getindex(fields, fld) for fld in flds]
end

function Base.getindex(fields::Fields{T}, flds::Tuple{Vararg{FieldIndex}}) where T<:Real
    [getindex(fields, fld) for fld in flds]
end

function Base.setindex!(fields::Fields{T}, val::U, flds::FieldIndex) where {T<:Real, U<:Real}
    if fld == R1
        fields.r1 = convert(T, val)
    elseif fld == R2
        fields.r2 = convert(T, val)
    elseif fld == C1
        fields.c1 = convert(T, val)
    elseif fld == C2
        fields.c2 = convert(T, val)
    elseif fld == C3
        fields.c3 = convert(T, val)
    elseif fld == C4
        fields.c4 = convert(T, val)
    elseif fld == I1
        fields.i1 = convert(T, val)
    else
        fields.i2 = convert(T, val)
    end
end

"""
    to_dual(fields::Fields{T}, flds::Symbol...) where T<:Real

Convert a fields object `fields` of type `T` into a fields object of
type `ForwardDiff.Dual{T}` with one's inserted in the locations of `flds`. i.e., if
`flds` = (:r1, :c1), one's are inserted into the first and third partials of
the Dual object. If a function is then called with the new fields object, we
will get derivatives of the function w.r.t. r1 and r2.

# Example
```julia
# Create a `Fields` object with random entries
julia> flds = Fields(rand(8))
# Convert to a dual object which can be used to take derivs. w.r.t r1
julia> to_dual(flds, :λ4)
# Convert to a dual object which can be used to take derivs. w.r.t c1 and i1
julia> to_dual(flds, :c1, :i1)
"""
function to_dual(fields::Fields{T}, flds::FieldIndex...) where T <: Real
    if length(unique(flds)) != length(flds)
        throw("'flds'=$flds passed to 'to_dual' contains a non-unique value.")
    end
    length(flds) <= 8 || throw("too many field symbols passed to 'to_dual'")

    ϵ = eps(T)

    Fields(
        R1 in flds ? Dual{T}(fields.r1 + ϵ, 1, 0, 0, 0, 0, 0, 0, 0) :
        Dual{T}(fields.r1, 0, 0, 0, 0, 0, 0, 0, 0),
        R2 in flds ? Dual{T}(fields.r2 + ϵ, 0, 1, 0, 0, 0, 0, 0, 0) :
        Dual{T}(fields.r2, 0, 0, 0, 0, 0, 0, 0, 0),
        C1 in flds ? Dual{T}(fields.c1 + ϵ, 0, 0, 1, 0, 0, 0, 0, 0) :
        Dual{T}(fields.c1, 0, 0, 0, 0, 0, 0, 0, 0),
        C2 in flds ? Dual{T}(fields.c2 + ϵ, 0, 0, 0, 1, 0, 0, 0, 0) :
        Dual{T}(fields.c2, 0, 0, 0, 0, 0, 0, 0, 0),
        C3 in flds ? Dual{T}(fields.c3 + ϵ, 0, 0, 0, 0, 1, 0, 0, 0) :
        Dual{T}(fields.c3, 0, 0, 0, 0, 0, 0, 0, 0),
        C4 in flds ? Dual{T}(fields.c4 + ϵ, 0, 0, 0, 0, 0, 1, 0, 0) :
        Dual{T}(fields.c4, 0, 0, 0, 0, 0, 0, 0, 0),
        I1 in flds ? Dual{T}(fields.i1 + ϵ, 0, 0, 0, 0, 0, 0, 1, 0) :
        Dual{T}(fields.i1, 0, 0, 0, 0, 0, 0, 0, 0),
        I2 in flds ? Dual{T}(fields.i2 + ϵ, 0, 0, 0, 0, 0, 0, 0, 1) :
        Dual{T}(fields.i2, 0, 0, 0, 0, 0, 0, 0, 0)
    )
end

"""
Indices for the THDM parameters. i.e. the parameter `λ1` is always in location
of `Int(Λ1)`, etc.
"""
@enum ParamIndex begin
    M112 = 1
    M122 = 2
    M222 = 3
    Λ1 = 4
    Λ2 = 5
    Λ3 = 6
    Λ4 = 7
    Λ5 = 8
end

"Parameters struct for the THDM parameters"
mutable struct Params{T<:Real}
    "coefficient of Φ̄₁Φ₁"
    m112::T
    "coefficient of Φ̄₁Φ₂ + Φ̄₂Φ₁"
    m122::T
    "coefficient of Φ̄₂Φ₂"
    m222::T
    "coefficient of (Φ̄₁Φ₁)²"
    λ1::T
    "coefficient of (Φ̄₂Φ₂)²"
    λ2::T
    "coefficient of Φ̄₁Φ₁Φ̄₂Φ₂"
    λ3::T
    "coefficient of Φ̄₂Φ₁Φ̄₁Φ₂"
    λ4::T
    "coefficient of (Φ̄₂Φ₁)² + (Φ̄₁Φ₂)²"
    λ5::T
    "renormalization scale"
    μ::T
    "top quark yukawa"
    yt::T # Top Yukawa
    "U(1) gauge coupling"
    gp::T # U(1)_Y gauge coupling
    "SU(2) gauge coupling"
    g::T # SU(2)_L gauge coupling
end

"""
    Params(params::Array{T,1}) where T<:Real

Construct a `Params` object from a vector of length 8
(for m112, m122, m222, λ1, λ2, λ3, λ4, λ5). Note, μ is set to the Higgs VEV.
"""
function Params(ps::Array{T,1}) where T<:Real
    if length(ps) == 8
        return Params(
            ps[1],
            ps[2],
            ps[3],
            ps[4],
            ps[5],
            ps[6],
            ps[7],
            ps[8],
            convert(T, HIGGS_VEV),
            zero(T),
            convert(T, U1Y_COUP),
            convert(T, SU2_COUP)
        )
    elseif length(ps) == 12
        Params(
            ps[1],
            ps[2],
            ps[3],
            ps[4],
            ps[5],
            ps[6],
            ps[7],
            ps[8],
            ps[9],
            ps[10],
            ps[11],
            ps[12]
        )
    else
        throw("length of input vector must be 8 or 12 when initializing params")
    end
end

"""
    Params{T}() where T<:Real

Construct an empty `Params` object that contains parameters of the THDM.
"""
function Params{T}() where T <: Real
    Params(
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        convert(T, U1Y_COUP),
        convert(T, SU2_COUP)
    )
end

function Base.getindex(params::Params{T}, par::ParamIndex) where T<:Real
    if par == M112
        return params.m112
    elseif par == M122
        return params.m122
    elseif par == M222
        return params.m222
    elseif par == Λ1
        return params.λ1
    elseif par == Λ2
        return params.λ2
    elseif par == Λ3
        return params.λ3
    elseif par == Λ4
        return params.λ4
    else
        return params.λ5
    end
end

function Base.getindex(params::Params{T}, pars::Array{ParamIndex,1}) where T<:Real
    [getindex(params, par) for par in pars]
end

function Base.getindex(params::Params{T}, pars::Tuple{Vararg{ParamIndex}}) where T<:Real
    [getindex(params, par) for par in pars]
end

function Base.setindex!(params::Params{T}, val::U, par::ParamIndex) where {T<:Real, U<:Real}
    if par == M112
        params.m112 = convert(T, val)
    elseif par == M122
        params.m122 = convert(T, val)
    elseif par == M222
        params.m222 = convert(T, val)
    elseif par == Λ1
        params.λ1 = convert(T, val)
    elseif par == Λ2
        params.λ2 = convert(T, val)
    elseif par == Λ3
        params.λ3 = convert(T, val)
    elseif par == Λ4
        params.λ4 = convert(T, val)
    else
        params.λ5 = convert(T, val)
    end
end

"""
    Params([μ=246.0,])

Construct a `Params` object that contains parameters of the THDM which are
bounded from below.
"""
function Params(μ::Float64=HIGGS_VEV)
    _is_bounded::Bool = false
    params::Params{Float64} = Params{Float64}()
    params.μ = μ
    while !_is_bounded
        params.m112 = 2μ^2 * (rand() - 0.5)
        params.m122 = 2μ^2 * (rand() - 0.5)
        params.m222 = 2μ^2 * (rand() - 0.5)
        params.λ1 = 10rand()
        params.λ2 = 10rand()

        geo_mean::Float64 = sqrt(params.λ1 * params.λ2)

        params.λ3 = (10 + geo_mean) * rand() - geo_mean
        params.λ4 = 2 * (rand() - 0.5)
        params.λ5 = 2 * (rand() - 0.5)

        _is_bounded = is_bounded(params)
    end
    params
end

"""
    is_bounded(params::Params)

Determine if the tree-level THDM potential is bounded from below.
"""
function is_bounded(params::Params{T}) where T<:Real
    bounded::Bool = false;

    if params.λ1 > 0 && params.λ2 > 0
        neg_geo_mean::T = -sqrt(params.λ1 * params.λ2)
        if params.λ3 >= neg_geo_mean
            if params.λ3 + params.λ4 - abs(params.λ5) >= neg_geo_mean
                bounded = true
            end
        end
    end
    bounded
end

"""
    to_dual(params::Params{T}, pars::Symbol...) where T<:Real

Convert a `Params` object `params` of type `T` into a `Params` object of
type `ForwardDiff.Dual{T}` with one's inserted in the locations of `pars`. i.e.,
if `pars` = (M112, Λ2), one's are inserted into the first and fourth partials
of the Dual object. If a function is then called with the new fields object, we
will get derivatives of the function w.r.t. m112 and λ1.

# Example
```julia
# Create a Parameters object with random entries
julia> pars = Params(rand(9))
# Convert to a dual object which can be used to take derivs. w.r.t λ4
julia> to_dual(pars, Λ1)
# Convert to a dual object which can be used to take derivs. w.r.t m112 and λ1
julia> to_dual(pars, M112, Λ4)
```
"""
function to_dual(params::Params{T}, pars::ParamIndex...) where T <: Real
    if length(unique(pars)) != length(pars)
        throw("'pars'=$pars passed to 'to_dual' contains a non-unique value.")
    end
    length(pars) <= 8 || throw("too many parameter symbols passed to 'to_dual'")

    Params(
        M112 in pars ? Dual{T}(params.m112, 1, 0, 0, 0, 0, 0, 0, 0) :
        Dual{T}(params.m112, 0, 0, 0, 0, 0, 0, 0, 0),
        M122 in pars ? Dual{T}(params.m122, 0, 1, 0, 0, 0, 0, 0, 0) :
        Dual{T}(params.m122, 0, 0, 0, 0, 0, 0, 0, 0),
        M222 in pars ? Dual{T}(params.m222, 0, 0, 1, 0, 0, 0, 0, 0) :
        Dual{T}(params.m222, 0, 0, 0, 0, 0, 0, 0, 0),
        Λ1 in pars ? Dual{T}(params.λ1, 0, 0, 0, 1, 0, 0, 0, 0) :
        Dual{T}(params.λ1, 0, 0, 0, 0, 0, 0, 0, 0),
        Λ2 in pars ? Dual{T}(params.λ2, 0, 0, 0, 0, 1, 0, 0, 0) :
        Dual{T}(params.λ2, 0, 0, 0, 0, 0, 0, 0, 0),
        Λ3 in pars ? Dual{T}(params.λ3, 0, 0, 0, 0, 0, 1, 0, 0) :
        Dual{T}(params.λ3, 0, 0, 0, 0, 0, 0, 0, 0),
        Λ4 in pars ? Dual{T}(params.λ4, 0, 0, 0, 0, 0, 0, 1, 0) :
        Dual{T}(params.λ4, 0, 0, 0, 0, 0, 0, 0, 0),
        Λ5 in pars ? Dual{T}(params.λ5, 0, 0, 0, 0, 0, 0, 0, 1) :
        Dual{T}(params.λ5, 0, 0, 0, 0, 0, 0, 0, 0),
        Dual{T}(params.μ, 0, 0, 0, 0, 0, 0, 0, 0),
        Dual{T}(params.yt, 0, 0, 0, 0, 0, 0, 0, 0),
        Dual{T}(params.g, 0, 0, 0, 0, 0, 0, 0, 0),
        Dual{T}(params.gp, 0, 0, 0, 0, 0, 0, 0, 0)
    )
end

"enum for the possible extrema type of a multi-dimensional function"
@enum ExtremaType begin
    Minimum
    Maximum
    Saddle
    Undefined
    NotSet
end

"""
Vacuum struct holding the VEVs of r₁, r₂, c₁, the potential at these VEVs and
the extremum type
"""
mutable struct Vacuum
    "VEV of r₁"
    v1::Float64
    "VEV of r₂"
    v2::Float64
    "VEV of c₁"
    α::Float64
    "value of the potential evaluated at v1, v2, α"
    potential::Float64
    "characterization of the extremum type of potential at VEVs"
    extremum_type::ExtremaType
end

"""
    generate_normal_vac(v::Float64=HIGGS_VEV)

Generate a random normal vacuum such that v₁² + v₂² = `v`^2 and α = 0.
"""
function generate_normal_vac(v::Float64=HIGGS_VEV)
    β::Float64 = π * rand() - π / 2
    v1::Float64 = v * cos(β)
    v2::Float64 = v * sin(β)
    Vacuum(v1, v2, 0.0, Inf, NotSet)
end

"""
    generate_normal_vac(v::Float64=HIGGS_VEV)

Generate a random CB vacuum such that v₁, v₂, α ∈ (v,-v).
"""
function generate_cb_vac(v::Float64=HIGGS_VEV)
    v1::Float64 = 2v * rand() - v
    v2::Float64 = 2v * rand() - v
    α::Float64 = 2v * rand() - v
    Vacuum(v1, v2, α, Inf, NotSet)
end

"""
    set_top_yukawa!(nvac::Vacuum, params::Params{Float64})

Set the value of the top quark Yukawa given the EW vacuum.
"""
function set_top_yukawa!(nvac::Vacuum, params::Params{Float64})
    params.yt = sqrt(2) * M_TOP / abs(nvac.v2);
    return
end

"""
    Fields(vac::Vacuum)

Create a `Fields` object from a vacuum.
"""
function Fields(vac::Vacuum)
    Fields(vac.v1, vac.v2, vac.α, 0.0, 0.0, 0.0, 0.0, 0.0)
end
