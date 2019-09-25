using LinearAlgebra: issymmetric

struct SymmEigMaxIterError <: Exception end
Base.showerror(io::IO, e::SymmEigMaxIterError) = print(io, "Too many iterations in `symmeig`")

"""
	pythag(a, b)

computes sqrt(a^2 + b^2) without destructive underflow or overflow.
"""
function pythag(a::T, b::T) where T <: Real
    absa = abs(a)
    absb = abs(b)
    if absa > absb
        return absa * sqrt(one(T) + (absb / absa)^2)
    else
        if absb == zero(T)
            return zero(T)
        else
            absb * sqrt(one(T) + (absa / absb)^2)
        end
    end
end

pythag(a::T, b::U) where {T<:Real,U<:Real} = pythag(promote(a, b)...)

"""
	samesign(a, b)

Return `a` with the same sign as `b`. Equivalent to abs(a) * sign(b)
"""
@inline samesign(a::T, b::T) where T <: Real =
    b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a)

"""
	householder_reduction!(z, d, e[;yesvecs=true])

Perform Householder reduction on a symmetric matrix `z` to reduce matrix into
a tridiagonal matrix. Values of the diagonal are stored in `d` and values of
off-diagonal are stored in `e`. If `yesvecs` is true, then computation will
proceed with this in mind.

# Arguments
`z::Array{T<:Real, 2}`: symmetric matrix
`d::Array{T<:Real, 1}`: array to store diagonal of tridiagonal matrix
`e::Array{T<:Real, 1}`: array to store off-diagonal of tridiagonal matrix
`yesvec::Bool=true`: flag indicating if eigen vectors are wanted.
"""
function householder_reduction!(
    z::Array{T,2},
    d::Array{T,1},
    e::Array{T,1};
    yesvecs = true
) where T <: Real
    n = size(z)[1]

    @inbounds for i = n:-1:2
        l = i - 1
        h = scale = zero(T)
        if l > 1
            for k = 1:(i-1)
                scale += abs(z[i, k])
            end
            if scale == zero(T)
                e[i] = z[i, l]
            else
                for k = 1:(i-1)
                    z[i, k] /= scale
                    h += z[i, k] * z[i, k]
                end
                f = z[i, l]
                g = (f >= zero(T) ? -sqrt(h) : sqrt(h))
                e[i] = scale * g
                h -= f * g
                z[i, l] = f - g
                f = zero(T)
                for j = 1:(i-1)
                    if (yesvecs)
                        z[j, i] = z[i, j] / h
                    end
                    g = zero(T)
                    for k = 1:j
                        g += z[j, k] * z[i, k]
                    end
                    for k = (j+1):(i-1)
                        g += z[k, j] * z[i, k]
                    end
                    e[j] = g / h
                    f += e[j] * z[i, j]
                end
                hh = f / (h + h)
                for j = 1:(i-1)
                    f = z[i, j]
                    e[j] = g = e[j] - hh * f
                    for k = 1:j
                        z[j, k] -= (f * e[k] + g * z[i, k])
                    end
                end
            end
        else
            e[i] = z[i, l]
        end
        d[i] = h
    end
    @inbounds if yesvecs
        d[1] = zero(T)
    end
    e[1] = zero(T)
    @inbounds for i = 1:n
        if yesvecs
            if d[i] != zero(T)
                for j = 1:(i-1)
                    g = zero(T)
                    for k = 1:(i-1)
                        g += z[i, k] * z[k, j]
                    end
                    for k = 1:(i-1)
                        z[k, j] -= g * z[k, i]
                    end
                end
            end
            d[i] = z[i, i]
            z[i, i] = one(T)
            for j = 1:(i-1)
                z[j, i] = z[i, j] = zero(T)
            end
        else
            d[i] = z[i, i]
        end
    end
end;

"""
	implict_ql!(z, d, e[;yesvecs=true])

Perform an implict QL decomposition on a tridiagonal matrix produced from
`householder_reduction!`. If `yesvecs` is true, the eigen vectors will be
stored in the columns of `z`. The eigenvalues are stored in `d`.

# Arguments
`z::Array{T<:Real, 2}`: symmetric matrix
`d::Array{T<:Real, 1}`: array containing diagonal of tridiagonal matrix
`e::Array{T<:Real, 1}`: array containg off-diagonal of tridiagonal matrix
`yesvec::Bool=true`: flag indicating if eigen vectors are wanted.
"""
function implict_ql!(
    z::Array{T,2},
    d::Array{T,1},
    e::Array{T,1};
    yesvecs = true
) where T <: Real
    n = size(z)[1]
    EPS = eps(T)
    @inbounds for i = 2:n
        e[i-1] = e[i]
    end
    e[n] = zero(T)
    @inbounds for l = 1:n
        iter = 0
        done = false
        while !done
            m = l
            while m <= n - 1
                dd = abs(d[m]) + abs(d[m+1])
                if (abs(e[m]) <= EPS * dd)
                    break
                end
                m += 1
            end
            if m != l
                iter += 1
                if iter == 30
                    throw(SymmEigMaxIterError())
                end
                g = (d[l+1] - d[l]) / 2e[l]
                r = pythag(g, 1)
                g = d[m] - d[l] + e[l] / (g + samesign(r, g))
                s = c = one(T)
                p = zero(T)
                for i = (m - 1):-1:l
                    f = s * e[i]
                    b = c * e[i]
                    e[i+1] = r = pythag(f, g)
                    if r == zero(T)
                        d[i+1] -= p
                        e[m] = zero(T)
                        break
                    end
                    s = f / r
                    c = g / r
                    g = d[i+1] - p
                    r = (d[i] - g) * s + 2c * b
                    p = s * r
                    d[i+1] = g + p
                    g = c * r - b
                    if yesvecs
                        for k = 1:n
                            f = z[k, i+1]
                            z[k, i+1] = s * z[k, i] + c * f
                            z[k, i] = c * z[k, i] - s * f
                        end
                    end
                end
                if r == zero(T) && i >= l
                    continue
                end
                d[l] -= p
                e[l] = g
                e[m] = zero(T)
            end
            if m == l
                done = true
            end
        end
    end
end

"""
	symmeig(M[;compute_evecs=false])

Compute the eigenvalues and, if `compute_evecs` is true, the eigenvectors of
the symmetric matrix `M`. First, the matrix is reduced to tridiagonal form
using householder reduction. Then, an implict QL decomposition is performed.

# Arguments
`M::Array{T<:Real, 2}`: symmetric matrix

# Example
```julia
julia> M = [1 2 3; 2 4 5; 3 5 6];
julia> symmeig(M)
[0.171, -0.516, 11.345]
julia> symmeig(M; compute_evecs=true)
([0.171, -0.516, 11.345], [0.591 -0.737 0.328; -0.737 -0.328 0.591; 0.328 0.591 0.737])
```
"""
function symmeig(M::Array{T,2}; compute_evecs = false) where T <: Real
    @assert issymmetric(M)
    z = deepcopy(M)
    n = size(z)[1]
    d = zeros(T, n)
    e = zeros(T, n)
    householder_reduction!(z, d, e; yesvecs = compute_evecs)
    implict_ql!(z, d, e; yesvecs = compute_evecs)
    compute_evecs ? (d, z) : d
end
