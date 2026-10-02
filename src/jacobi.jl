import StaticArrays
include(joinpath(@__DIR__, "common.jl"))

export sn, cn, dn, nn, sd, cd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns
sn(u, m) = ellipj(u, m)[1]
cn(u, m) = ellipj(u, m)[2]
dn(u, m) = ellipj(u, m)[3]

xn = ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, :(1)))
for (p, num) in xn, (q, den) in xn
    f = Symbol(p, q)
    #@eval begin
    #    """
    #        $($f)(u::Real, m::Real)

    #    Compute the Jacobi elliptic function $($f)(u | m)
    #    """
    #    ($f)(u::Real, m::Real) = ($f)(Float64(u), Float64(m))
    #end

    if (p == q)
        @eval ($f)(::A, ::B) where {A,B} = begin
            one(promote_type(A, B))
        end
    elseif (q != :n)
        @eval ($f)(u::A, m::B) where {A,B} = begin
            ($num) / ($den)
        end
    end
end

#end # module
