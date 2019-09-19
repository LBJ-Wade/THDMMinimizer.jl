using SpecialFunctions
using ForwardDiff: Dual

function test_matrix1(a::T, b::T) where T <: Real
    [
     exp(a + b) log(a^2 + b^2);
     log(a^2 + b^2) (a * b + a^2 * sin(b))
    ]
end

function test_matrix2(a::T, b::T) where T <: Real
    [
     sinh(sqrt(a * b)) besselk(2, a);
     besselk(2, a) a^2 + b^2 + a + b + one(T)
    ]
end

function test_matrix3(a::T, b::T) where T <: Real
    [
     a b;
     b a
    ]
end

function test_matrix1_evals(a::T, b::T) where T <: Real
    M11 = exp(a + b)
    M12 = log(a^2 + b^2)
    M22 = a * b + a^2 * sin(b)
    eval1 = 1 / 2 * (M11 + M22 - sqrt((M11 - M22)^2 + 4M12^2))
    eval2 = 1 / 2 * (M11 + M22 + sqrt((M11 - M22)^2 + 4M12^2))
    [eval1, eval2]
end

function test_matrix2_evals(a::T, b::T) where T <: Real
    M11 = sinh(sqrt(a * b))
    M12 = besselk(2, a)
    M22 = a^2 + b^2 + a + b + one(T)
    eval1 = 1 / 2 * (M11 + M22 - sqrt((M11 - M22)^2 + 4M12^2))
    eval2 = 1 / 2 * (M11 + M22 + sqrt((M11 - M22)^2 + 4M12^2))
    [eval1, eval2]
end

function test_matrix3_evals(a::T, b::T) where T <: Real
    M11 = a
    M12 = b
    M22 = a
    eval1 = 1 / 2 * (M11 + M22 - sqrt((M11 - M22)^2 + 4M12^2))
    eval2 = 1 / 2 * (M11 + M22 + sqrt((M11 - M22)^2 + 4M12^2))
    [eval1, eval2]
end

function test_matrix1_eval_derivs(a::T, b::T, deriv::Symbol) where T <: Real
    M11 = exp(a + b)
    M12 = log(a^2 + b^2)
    M22 = a * b + a^2 * sin(b)
    if deriv == :a
        dM11 = exp(a + b)
        dM12 = 2a / (a^2 + b^2)
        dM22 = 2a * sin(b) + b
    elseif deriv == :b
        dM11 = exp(a + b)
        dM12 = 2b / (a^2 + b^2)
        dM22 = a^2 * cos(b) + a
    end
    d1 = (dM11 + dM22 +
          ((-dM11 + dM22) * M11 - 4 * dM12 * M12 + (dM11 - dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    d2 = (dM11 + dM22 +
          ((dM11 - dM22) * M11 + 4 * dM12 * M12 + (-dM11 + dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    e1, e2 = test_matrix1_evals(a, b)
    [Dual{Float64}(e1, d1), Dual{Float64}(e2, d2)]
end

function test_matrix2_eval_derivs(a::T, b::T, deriv::Symbol) where T <: Real
    M11 = sinh(sqrt(a * b))
    M12 = besselk(2, a)
    M22 = a^2 + b^2 + a + b + one(T)
    if deriv == :a
        dM11 = b * cosh(sqrt(a * b)) / (2sqrt(a * b))
        dM12 = -(besselk(1, a) + besselk(3, a)) / 2
        dM22 = 2a + one(T)
    elseif deriv == :b
        dM11 = a * cosh(sqrt(a * b)) / (2sqrt(a * b))
        dM12 = zero(T)
        dM22 = 2b + one(T)
    end
    d1 = (dM11 + dM22 +
          ((-dM11 + dM22) * M11 - 4 * dM12 * M12 + (dM11 - dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    d2 = (dM11 + dM22 +
          ((dM11 - dM22) * M11 + 4 * dM12 * M12 + (-dM11 + dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    e1, e2 = test_matrix2_evals(a, b)
    [Dual{Float64}(e1, d1), Dual{Float64}(e2, d2)]
end

function test_matrix3_eval_derivs(a::T, b::T, deriv::Symbol) where T <: Real
    M11 = a
    M12 = b
    M22 = a
    if deriv == :a
        dM11 = one(T)
        dM12 = zero(T)
        dM22 = one(T)
    elseif deriv == :b
        dM11 = zero(T)
        dM12 = one(T)
        dM22 = zero(T)
    end
    d1 = (dM11 + dM22 +
          ((-dM11 + dM22) * M11 - 4 * dM12 * M12 + (dM11 - dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    d2 = (dM11 + dM22 +
          ((dM11 - dM22) * M11 + 4 * dM12 * M12 + (-dM11 + dM22) * M22) /
          sqrt(M11^2 + 4 * M12^2 - 2 * M11 * M22 + M22^2)) / 2.
    e1, e2 = test_matrix3_evals(a, b)
    [Dual{Float64}(e1, d1), Dual{Float64}(e2, d2)]
end
