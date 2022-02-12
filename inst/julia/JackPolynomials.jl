module JackPolynomials

import DynamicPolynomials
export JackPolynomial
export ZonalPolynomial
export SchurPolynomial
export Zonal
export Schur
export Jack

function isPartition(lambda::Vector{<:Integer})
  return all(diff(lambda) .<= 0) && all(lambda .>= 0)
end

function dualPartition(lambda::Vector{T}) where {T<:Integer}
  out = T[]
  if !isempty(lambda)
    for i = 1:lambda[1]
      push!(out, sum(lambda .>= i))
    end
  end
  return out
end

function betaratio(
  kappa::Vector{I},
  mu::Vector{I},
  k::I,
  alpha::T,
) where {T<:Real,I<:Integer}
  muk = mu[k]
  t = k - alpha * muk
  u = map(i -> t + 1 - i + alpha * kappa[i], 1:k)
  v = map(i -> t - i + alpha * mu[i], 1:(k-1))
  muPrime = dualPartition(mu)
  w = map(i -> muPrime[i] - t - alpha * i, 1:(muk-1))
  prod1 = prod(u ./ (u .+ alpha .- 1))
  prod2 = prod((v .+ alpha) ./ v)
  prod3 = prod((w .+ alpha) ./ w)
  return alpha * prod1 * prod2 * prod3
end

function _N(lambda::Vector{I}, mu::Vector{I}) where {I<:Integer}
  prods = map(i -> prod(Iterators.drop(lambda .+ 1, i)), 1:length(lambda))
  return sum(mu .* prods)
end

# ------------------------------------------------------------------------------
#~~ Jack polynomial ~~~~##
# ------------------------------------------------------------------------------
"""
    Jack(x, lambda, alpha)

Evaluates a Jack polynomial.

# Arguments
- `x`: vector of real or complex numbers
- `lambda`: partition of an integer
- `alpha`: alpha parameter
"""
function Jack(
  x::Vector{C},
  lambda::Vector{I},
  alpha::T,
) where {T<:Real,I<:Integer,C<:Number}
  if !isPartition(lambda)
    error("`lambda` must be a partition of an integer")
  end
  if alpha <= 0
    error("`alpha` must be positive")
  end
  function jac(m::I, k::I, mu::Vector{I}, nu::Vector{I}, beta::Real)
    if isempty(nu) || nu[1] == 0 || m == 0
      return C(1)
    end
    if length(nu) > m && nu[m+1] > 0
      return C(0)
    end
    if m == 1
      return x[1]^(nu[1]) * prod(1 .+ alpha .* collect(1:(nu[1]-1)))
    end
    v = S[_N(lambda, nu), m]
    if k == 0 && !ismissing(v)
      return v
    end
    i = max(1, k)
    s = jac(m - 1, 0, nu, nu, T(1)) * beta * x[m]^(sum(mu) - sum(nu))
    while length(nu) >= i && nu[i] > 0
      if length(nu) == i || nu[i] > nu[i+1]
        nuPrime = copy(nu)
        nuPrime[i, 1] -= 1
        gamma = beta * betaratio(mu, nu, i, alpha)
        if nu[i] > 1
          s = s + jac(m, i, mu, nuPrime, gamma)
        else
          s =
            s +
            jac(m - 1, 0, nuPrime, nuPrime, T(1)) * gamma * x[m]^(sum(mu) - sum(nuPrime))
        end
      end
      i += 1
    end
    if k == 0
      S[_N(lambda, nu), m] = s
    end
    return s
  end # end jac --------------------------------------------------------------
  lx = length(x)
  S = Array{Union{Missing,C}}(missing, _N(lambda, lambda), lx)
  jac(lx, 0, lambda, lambda, T(1))
end

# ------------------------------------------------------------------------------
#~~ Symbolic Jack polynomial ~~~~##
# ------------------------------------------------------------------------------
function JackPolynomial0(
  m::I, lambda::Vector{I}, alpha::T
) where {T<:Real,I<:Integer}
  function jac(m::I, k::I, mu::Vector{I}, nu::Vector{I}, beta::Real)
    if isempty(nu) || nu[1] == 0 || m == 0
      return T(1)
    end
    lnu = length(nu)
    if lnu > m && nu[m+1] > 0
      return T(0)
    end
    if m == 1
      return x[1]^(nu[1]) * prod(1 .+ alpha .* collect(1:(nu[1]-1)))
    end
    v = S[_N(lambda, nu), m]
    if k == 0 && !ismissing(v)
      return v
    end
    i = max(1, k)
    s = jac(m - 1, 0, nu, nu, T(1)) * beta * x[m]^(sum(mu) - sum(nu))
    while lnu >= i && nu[i] > 0
      if lnu == i || nu[i] > nu[i+1]
        nuPrime = copy(nu)
        nuPrime[i, 1] -= 1
        gamma = beta * betaratio(mu, nu, i, alpha)
        if nu[i] > 1
          s += jac(m, i, mu, nuPrime, gamma)
        else
          s += jac(m - 1, 0, nuPrime, nuPrime, T(1)) * gamma * x[m]^(sum(mu) - sum(nuPrime))
        end
      end
      i += 1
    end
    if k == 0
      S[_N(lambda, nu), m] = s
    end
    return s
  end
  DynamicPolynomials.@polyvar x[1:m]
  S = Array{Union{Missing,DynamicPolynomials.Polynomial{true,T}}}(
    missing,
    _N(lambda, lambda),
    m,
  )
  jac(m, 0, lambda, lambda, T(1))
end

"""
    JackPolynomial(x, lambda, alpha)

Symbolic Jack polynomial. The coefficients of the polynomial will have the
same type as `alpha`.

# Arguments
- `m`: integer, the number of variables
- `lambda`: partition of an integer
- `alpha`: alpha parameter
"""
function JackPolynomial(
  m::I, lambda::Vector{I}, alpha::T, R::Bool=false
) where {T<:Real,I<:Integer}
  if !isPartition(lambda)
    error("`lambda` must be a partition of an integer")
  end
  if alpha <= 0
    error("`alpha` must be positive")
  end
  jack = JackPolynomial0(m, lambda, alpha)
  if(typeof(jack) == T)
    DynamicPolynomials.@polyvar x[1:m]
    jack = sum(T(0) * x) + jack
  end
  if R
    return (
      coefficients = jack.a,
      powers = jack.x.Z
    )
  end
  return jack
end

# ------------------------------------------------------------------------------
#~~ Zonal polynomial ~~~~##
# ------------------------------------------------------------------------------
function _i(lambda::Vector{T}) where {T<:Integer}
  out = T[]
  for i in 1:length(lambda)
    out = vcat(out, repeat([i], lambda[i]))
  end
  return out
end

function _j(lambda::Vector{T}) where {T<:Integer}
  out = T[]
  mu = filter(x -> x > 0, lambda)
  for m in mu
    out = vcat(out, 1:m)
  end
  return out
end

function hookLengths(lambda::Vector{<:Integer}, alpha::Real)
  i = _i(lambda)
  j = _j(lambda)
  lambdap = dualPartition(lambda)
  upper = lambdap[j] .- i .+ alpha .* (lambda[i] .- j .+ 1)
  lower = lambdap[j] .- i .+ 1 .+ alpha .* (lambda[i] .- j)
  return vcat(upper, lower)
end

"""
    Zonal(x, lambda)

Evaluates a zonal polynomial.

# Arguments
- `x`: vector of real or complex numbers
- `lambda`: partition of an integer
"""
function Zonal(
  x::Vector{<:Union{R,Complex{R}}},
  lambda::Vector{<:Integer},
) where {R<:Real}
  jack = Jack(x, lambda, R(2))
  jlambda = prod(hookLengths(lambda, R(2)))
  n = sum(lambda)
  return jack * 2^n * factorial(n) / jlambda
end

# ------------------------------------------------------------------------------
#~~ Symbolic zonal polynomial ~~~~##
# ------------------------------------------------------------------------------
"""
    ZonalPolynomial(m, lambda[, type])

Symbolic zonal polynomial.

# Arguments
- `m`: integer, the number of variables
- `lambda`: partition of an integer
- `type`: the type of the coefficients of the polynomial; default `Rational`
"""
function ZonalPolynomial(
  m::I,
  lambda::Vector{I},
  type::Type = Real
) where {I<:Integer}
  jack = JackPolynomial(m, lambda, 2//1)
  jlambda = prod(hookLengths(lambda, 2//1))
  n = sum(lambda)
  poly = jack * 2^n * factorial(n) / jlambda
  return (
    qcoefficients = poly.a,
    coefficients = convert(Vector{Float64}, poly.a),
    powers = poly.x.Z
  )
end


# ------------------------------------------------------------------------------
#~~ Schur polynomial ~~~~##
# ------------------------------------------------------------------------------
"""
    Schur(x, lambda)

Evaluates a Schur polynomial.

# Arguments
- `x`: vector of real or complex numbers
- `lambda`: partition of an integer
"""
function Schur(x::Vector{T}, lambda::Vector{I}) where {T<:Number,I<:Integer}
  if !isPartition(lambda)
    error("`lambda` must be a partition of an integer")
  end
  function sch(m::Int64, k::Int64, nu::Vector{I})
    if isempty(nu) || nu[1] == 0 || m == 0
      return T(1)
    end
    if length(nu) > m && nu[m+1] > 0
      return T(0)
    end
    if m == 1
      return x[1]^nu[1]
    end
    v = S[_N(lambda, nu), m]
    if !ismissing(v)
      return v
    end
    s = sch(m - 1, 1, nu)
    lnu = length(nu)
    i = k
    while lnu >= i && nu[i] > 0
      if lnu == i || nu[i] > nu[i+1]
        nup = copy(nu)
        nup[i] -= 1
        if nu[i] > 1
          s += x[m] * sch(m, i, nup)
        else
          s += x[m] * sch(m - 1, 1, nup)
        end
      end
      i += 1
    end
    if k == 1
      S[_N(lambda, lambda), m] = s
    end
    return s
  end # end sch --------------------------------------------------------------
  lx = length(x)
  S = Array{Union{Missing,T}}(missing, _N(lambda, lambda), lx)
  sch(lx, 1, lambda)
end

# ------------------------------------------------------------------------------
#~~ Symbolic Schur polynomial ~~~~##
# ------------------------------------------------------------------------------
function SchurPolynomial0(
  m::I,
  lambda::Vector{I},
  T::Type
) where {I<:Integer}
  function sch(m::Int64, k::Int64, nu::Vector{I})
    if isempty(nu) || nu[1] == 0 || m == 0
      return T(1)
    end
    if length(nu) > m && nu[m+1] > 0
      return T(0)
    end
    if m == 1
      return x[1]^nu[1]
    end
    v = S[_N(lambda, nu), m]
    if !ismissing(v)
      return v
    end
    s = sch(m - 1, 1, nu)
    lnu = length(nu)
    i = k
    while lnu >= i && nu[i] > 0
      if lnu == i || nu[i] > nu[i+1]
        nup = copy(nu)
        nup[i] -= 1
        if nu[i] > 1
          s += x[m] * sch(m, i, nup)
        else
          s += x[m] * sch(m - 1, 1, nup)
        end
      end
      i += 1
    end
    if k == 1
      S[_N(lambda, lambda), m] = s
    end
    return s
  end # end sch --------------------------------------------------------------
  DynamicPolynomials.@polyvar x[1:m]
  S = Array{Union{Missing,DynamicPolynomials.Polynomial{true,T}}}(
    missing,
    _N(lambda, lambda),
    m
  )
  sch(m, 1, lambda)
end

"""
    SchurPolynomial(m, lambda[, type])

Symbolic Schur polynomial.

# Arguments
- `m`: integer, the number of variables
- `lambda`: partition of an integer
- `type`: the type of the coefficients of the polynomial; default `Int64`
"""
function SchurPolynomial(
  m::I,
  lambda::Vector{I},
  type::Type = Int64
) where {I<:Integer}
  if !isPartition(lambda)
    error("`lambda` must be a partition of an integer")
  end
  schur = SchurPolynomial0(m, lambda, type)
  if(typeof(schur) == type)
    DynamicPolynomials.@polyvar x[1:m]
    schur = sum(type(0) * x) + schur
  end
  return (
    coefficients = schur.a,
    powers = schur.x.Z
  )
end

end # module
