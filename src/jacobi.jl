struct JacobiMaxIterError <: Exception end
Base.showerror(io::IO, e::JacobiMaxIterError) = print(io, "Too many iterations in `jacobi`")


@inline function rot!(
    mat::Array{T,2},
    s::T,
    τ::T,
    i::Int,
    j::Int,
    k::Int,
    l::Int
) where T <: Real
    g::T = mat[i, j]
    h::T = mat[k, l]
    @inbounds mat[i, j] = g - s * (h + g * τ)
    @inbounds mat[k, l] = h + s * (g - h * τ)
end

"""
    jacobi(matrix)

Compute the eigenvalues of a square, symmetric matrix.
"""
function jacobi(mmat::Array{T,2}) where T <: Real
    # Copy the matrix so we don't destroy original
    mat = deepcopy(mmat)
    n::Int = size(mat)[1]

    v::Array{T,2} = Matrix{T}(I, n, n)
    d::Array{T,1} = diag(mat)
    b::Array{T,1} = diag(mat)
    z::Array{T,1} = zeros(T, n)

    ϵ::T = eps(T)
    JACOBI_MAX_ITER = 100
    for i = 1:JACOBI_MAX_ITER
        sm::T = zero(T)
        for ip = 1:(n-1)
            for iq = (ip+1):n
                @inbounds sm += abs(mat[ip, iq])
            end
        end
        if sm == 0
            return d
        end

        thresh::T = i < 4 ? 0.2 * sm / n^2 : zero(T)

        for ip = 1:(n-1)
            for iq = (ip+1):n
                g::T = 100abs(mat[ip, iq])
                if i > 4 && g <= ϵ * abs(d[ip]) && g <= ϵ * abs(d[iq])
                     @inbounds mat[ip, iq] = zero(T)
                elseif abs(mat[ip, iq]) > thresh
                    @inbounds h::T = d[iq] - d[ip]
                    t::T = zero(T)
                    if g <= ϵ * abs(h)
                         @inbounds t = mat[ip, iq] / h
                    else
                        @inbounds θ::T = 0.5 * h / mat[ip, iq]
                        t = 1.0 / (abs(θ) + sqrt(1.0 + θ^2))
                        if θ < 0
                            t = -t
                        end
                    end
                    c::T = 1.0 / sqrt(1 + t^2)
                    s::T = t * c
                    τ::T = s / (1 + c)
                    @inbounds h = t * mat[ip, iq]
                    @inbounds z[ip] -= h
                    @inbounds z[iq] += h
                    @inbounds d[ip] -= h
                    @inbounds  d[iq] += h
                    @inbounds mat[ip, iq] = zero(T)

                    for j = 1:(ip-1)
                        rot!(mat, s, τ, j, ip, j, iq)
                    end
                    for j = (ip+1):(iq-1)
                        rot!(mat, s, τ, ip, j, j, iq)
                    end
                    for j = (iq+1):n
                        rot!(mat, s, τ, ip, j, iq, j)
                    end
                    for j = 1:n
                        rot!(v, s, τ, j, ip, j, iq)
                    end
                end
            end
        end

        for ip = 1:n
            @inbounds b[ip] += z[ip]
            @inbounds d[ip] = b[ip]
            @inbounds z[ip] = zero(T)
        end
    end
    throw(JacobiMaxIterError())
end
