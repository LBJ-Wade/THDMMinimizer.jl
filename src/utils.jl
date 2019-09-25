"""
  horner(x, p...)

Convert an expression of the form a₀ + a₁x + a₂x² + a₃x³ + ... into Horner form:
  a₀ + a₁x + a₂x² -> a₀ + x(a₁ + x(a₂ + x(a₃ + ⋅⋅⋅)
and combines all expressions of the form `ax + y` operatations into `muladd`.

# Example
horner(x, 1, 2, 3) = 1 + 2x + 3x² -> muladd(x, muladd(3,x,2), 1)
"""
macro horner(x, p...)
  ex = esc(p[end])
  for i = length(p) - 1:-1:1
    ex = :(muladd(t, $ex, $(esc(p[i]))))
  end
  Expr(:block, :(t = $(esc(x))), ex)
end
