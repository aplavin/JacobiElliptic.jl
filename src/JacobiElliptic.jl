module JacobiElliptic
using StaticArrays
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         """


export Fukushima, Carlson
export am, sn, cn, dn, nn, sd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns # Do not export cd to avoid name clash with Base.cd
export ellipj, ellipj_add, ellipke
export E, F, K, Pi, J
include(joinpath(@__DIR__, "Fukushima.jl"))
include(joinpath(@__DIR__, "Carlson.jl"))
include(joinpath(@__DIR__, "AGM.jl"))

abstract type AbstractAlgorithm end
struct Fukushima <: AbstractAlgorithm end
struct Carlson <: AbstractAlgorithm end
struct AGM <: AbstractAlgorithm end


func_syms = [
    :E,
    :F,
    :K,
    :Pi,
    :J,
    :sn,
    :cn,
    :dn,
    :nn,
    :sd,
    :dd,
    :nd,
    :sc,
    :cc,
    :dc,
    :nc,
    :ss,
    :cs,
    :ds,
    :ns,
    :cd,
]

# Creates function func(::alg, args...) = algAlg.sym(args...)
for alg in [:Fukushima, :Carlson]
    for func in func_syms
        f = Symbol(alg, :Alg)
        @eval $func(::$alg, args...) = $f.$func(args...)
    end
end

K(alg::AGM, m) = ArithmeticGeometricMeanAlg.K(m)
E(alg::AGM, m) = ArithmeticGeometricMeanAlg.E(m)


# Sets default functions to Carlson Alg
# func(args...) = CarlsonAlg.func(args...)
alg = :Carlson
for func in func_syms
    f = Symbol(alg, :Alg)
    @eval $func(args...) = $f.$func(args...)
end

K(m) = ArithmeticGeometricMeanAlg.K(m)
E(m) = ArithmeticGeometricMeanAlg.E(m)

ellipj(u, m) = CarlsonAlg.ellipj(u, m)
ellipj(::Carlson, u, m) = CarlsonAlg.ellipj(u, m)
ellipj(u, m, K) = CarlsonAlg.ellipj(u, m, K)

"""
    ellipj_add((sn₁, cn₁, dn₁), (sn₂, cn₂, dn₂), m) -> (sn, cn, dn)

Jacobi elliptic functions at `u₁ + u₂` for `m ≤ 1`, given `ellipj(u₁, m)` and `ellipj(u₂, m)`, by the
addition theorem (A&S 16.17.1–2):

    D  = 1 - m sn₁² sn₂²
    sn = (sn₁ cn₂ dn₂ + sn₂ cn₁ dn₁) / D
    cn = (cn₁ cn₂ - sn₁ sn₂ dn₁ dn₂) / D

Intended for repeated stepping `u → u + Δu` with a fixed `ellipj(Δu, m)`. To keep such chains
from drifting, `(sn, cn)` is renormalized onto `sn² + cn² = 1` by one Newton step (a norm residual
`e` becomes `O(e²)`), and `dn = √(cn² + (1-m) sn²)` (the positive root, no cancellation near
`m = 1`). For `m > 1`, where `dn` changes sign, `dn` is `NaN`. Phase error along a chain is not
corrected.
"""
@inline function ellipj_add(j1::NTuple{3,T}, j2::NTuple{3,T}, m::T) where {T<:Real}
    s1, c1, d1 = j1
    s2, c2, d2 = j2
    invD = inv(one(T) - m * s1 * s1 * s2 * s2)
    s = (s1 * c2 * d2 + s2 * c1 * d1) * invD
    c = (c1 * c2 - s1 * s2 * d1 * d2) * invD
    # one Newton step for the inverse norm: residual e = s² + c² - 1 → -3e²/4 + e³/4
    v = (T(3) - s * s - c * c) / T(2)
    s *= v
    c *= v
    d = sqrt(c * c + (one(T) - m) * s * s)
    return (s, c, ifelse(m > one(T), T(NaN), d))
end

asn = FukushimaAlg.asn
acn = FukushimaAlg.acn
am = CarlsonAlg.am
end
