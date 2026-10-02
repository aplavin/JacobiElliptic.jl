export am

# Max AGM steps: Float64 needs ≤ 9 even at m = prevfloat(1.0), Float32 ≤ 6; wider types more
_am_maxsteps(::Type{<:Union{Float16,Float32,Float64}}) = 10
_am_maxsteps(::Type) = 32

# A&S 16.4: descending AGM for 0 ≤ m ≤ 1, absolute error O(eps·u) uniformly in m, so it is used up
# to m → 1. No small-(1-m) asymptotic (A&S 16.15.4): its error grows like (1-m)·cosh(u), wrong for u ~ K(m).
function _am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)

    # Use an immutable local buffer so GPU backends can keep this in registers.
    nmax = _am_maxsteps(T)
    _ambuf = zero(StaticArrays.SVector{nmax,T})
    u == 0 && return zero(T)

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4 (O(m²) remainder)
        return u - m * (u - sin(2 * u) / 2) / 4
    end
    m1 = one(T) - m
    # m = 1: am(u|1) = gd(u)
    iszero(m1) && return asin(tanh(u))

    a = one(T)
    b = sqrt(m1)
    c = sqrt(m)
    n = 0
    scale = one(T)  # 2^n, accumulated exactly (no integer power: GPU-safe)
    while abs(c) > tol && n < nmax
        a, b, c, n = ((a + b) / 2, sqrt(a * b), (a - b) / 2, n + 1)
        scale *= 2
        @inbounds _ambuf = Base.setindex(_ambuf, c / a, n)
    end
    !(abs(c) <= tol) && return T(NaN)  # also NaN m, for which the loop above never runs

    phi = a * u * scale
    for i = n:-1:1
        @inbounds phi = (phi + asin(_ambuf[i] * sin(phi))) / 2
    end
    phi
end
function _am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return _am(u, m, eps(T))
end

"""
    am(u::Real, m::Real)

Returns amplitude, φ, such that u = F(φ | m)
"""
function am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    if m < 0
        mu1 = inv(1 - m)
        mu = -m * mu1
        sqrtmu1 = sqrt(mu1)
        v = u / sqrtmu1
        phi = _am(v, mu)
        s = sin(phi)
        t = floor((phi + A(π / 2)) / A(π))

        return t * π + ((-1)^t) * asin(sqrtmu1 * s / sqrt(1 - mu * s^2))
    elseif m <= 1 # 0 <= m <= 1
        return _am(u, m, eps(T))
    else # m > 1
        # Reciprocal Modulus Transformation
        # P.F. Byrd & M.D. Friedman: Handbook of Elliptic Integrals for Engineers and Scientists, 1971, p. 38, eq. 162.01
        # am(u,k) = arcsin(sin(am(u,k))) = arcsin(sn(u,k)) = arcsin(1 / k * sn(k * u, 1 / k))
        k = sqrt(m)
        return asin(inv(k) * sn(k * u, inv(m)))
    end
end

@inline function custom_atanh(a::T) where {T}
    oneT = one(T)
    arg1 = abs(oneT + a)
    arg2 = abs(oneT - a)

    ans = (log(arg1 / arg2)) / 2
    return ans
end

@inline function FukushimaT(t::A, h::B) where {A,B}
    T = promote_type(A, B)
    if h > zero(T)
        sqrt_h = √h
        return atan(t * sqrt_h) / sqrt_h
    elseif h == zero(T)
        return t
    else
        sqrt_neg_h = √(-h)
        arg = t * sqrt_neg_h
        ans = abs(arg) < one(T) ? atanh(arg) : custom_atanh(arg)
        return ans / sqrt_neg_h
    end
end

#https://link-springer-com.ezp-prod1.hul.harvard.edu/article/T(10).1007/BF02165405
function cel(kc::A, p::B, a::C, b::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    #ca = T(1e-6)
    ca = eps(T)
    oneT = one(T)
    twoT = T(2)
    pi_over_2 = T(π / 2)
    kc = abs(kc)
    e = kc
    m = oneT

    f, g, q = T(0), T(0), T(0)
    if p > T(0)
        p = √p
        b = b / p
    else
        f = kc^2
        q = oneT - f
        g = oneT - p
        f = f - p
        q = (b - a * p) * q
        p = √(f / g)
        a = (a - b) / g
        b = -q * (g^2 * p) + a * p
    end
    count = 0
    while count < 1000
        count+=1
        f = a
        invp = inv(p)
        a = muladd(invp, b, a)
        g = e * invp
        b = twoT * muladd(f, g, b)
        p = g + p
        g = m
        m = kc + m
        !(abs(g - kc) >= g * ca) && break  # `!(≥)`: NaN terminates
        kc = twoT * √e
        e = kc * m
    end
    return pi_over_2 * muladd(a, m, b) / (m * (m + p))
end
