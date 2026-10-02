import StaticArrays
include(joinpath(@__DIR__, "common.jl"))

export sn, cn, dn, nn, sd, cd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns
sn(u, m) = ellipj(u, m)[1]
cn(u, m) = ellipj(u, m)[2]
dn(u, m) = ellipj(u, m)[3]

# pq(u, m) = p(u, m) / q(u, m) with n ≡ 1, from a single ellipj evaluation
xn = ((:s, 1), (:c, 2), (:d, 3), (:n, 0))
for (p, i) in xn, (q, j) in xn
    f = Symbol(p, q)
    if p == q
        @eval ($f)(::A, ::B) where {A,B} = one(promote_type(A, B))
    elseif q != :n
        @eval ($f)(u, m) = (scd = ellipj(u, m); $(i == 0 ? :(inv(scd[$j])) : :(scd[$i] / scd[$j])))
    end
end

