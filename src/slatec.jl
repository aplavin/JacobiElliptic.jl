#module SLATEC

export DRC, DRD, DRF, DRJ

#***BEGIN PROLOGUE  DRF
#***PURPOSE  Compute the incomplete or complete elliptic integral of the
#            1st kind.  For X, Y, and Z non-negative and at most one of
#            them zero, RF(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -1/2
#                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X, Y or Z is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

@inline function DRF(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)

    ERRTOL = (4 * eps(T) / 2)^T(1 / 6)
    LOLIM = 5floatmin(T)
    UPLIM = floatmax(T) / 5
    oneT = one(T)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)
    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)
    C0 = T(1 / 10)

    ans = zero(T)
    min(X, Y, Z) < zero(T) && return (ans, 1)
    max(X, Y, Z) > UPLIM && return (ans, 3)
    min(X + Y, X + Z, Y + Z) < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    ZN = Z
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    while true
        MU = (XN + YN + ZN) * inv3
        ninvMU = -inv(MU)
        XNDEV = muladd(ninvMU, (MU + XN), 2)
        YNDEV = muladd(ninvMU, (MU + YN), 2)
        ZNDEV = muladd(ninvMU, (MU + ZN), 2)
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        (EPSLON < ERRTOL) && break
        XNROOT = _sqrt(XN)
        YNROOT = _sqrt(YN)
        ZNROOT = _sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, (YNROOT + ZNROOT), YNROOTZNROOT)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
    end
    XNDEVYNDEV = XNDEV * YNDEV
    E2 = muladd(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV * ZNDEV
    S = oneT + muladd(E2, muladd(-C2, E3, muladd(C1, E2, -C0)), C3 * E3)
    ans = S / _sqrt(MU)

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRD
#***PURPOSE  Compute the incomplete or complete elliptic integral of
#            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
#            DRD(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -3/2
#                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X or Y is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

@inline function DRD(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)

    ERRTOL = (eps(T) / 6)^T(1 / 6)
    LOLIM = 2 / (floatmax(T))^T(2 / 3)
    TUPLIM = floatmin(T)^T(1 / 3)
    TUPLIM = (ERRTOL / 10)^T(1 / 3) / TUPLIM
    UPLIM = TUPLIM^2
    oneT = one(T)
    twoT = T(2)
    threeT = T(3)
    sixT = T(6)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
    C1 = T(3 / 14)
    C2 = T(1 / 6)
    C3 = T(9 / 22)
    C4 = T(3 / 26)

    ans = zero(T)
    min(X, Y) < 0 && return (ans, 1)
    max(X, Y, Z) > UPLIM && return (ans, 3)
    min(X + Y, Z) < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    ZN = Z
    SIGMA = zero(T)
    POWER4 = one(T)
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    while true
        MU = (XN + YN + threeT * ZN) * inv5
        invMU = inv(MU)
        XNDEV = (MU - XN) * invMU
        YNDEV = (MU - YN) * invMU
        ZNDEV = (MU - ZN) * invMU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        (EPSLON < ERRTOL) && break
        XNROOT = _sqrt(XN)
        YNROOT = _sqrt(YN)
        ZNROOT = _sqrt(ZN)
        ynrootznroot = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, ynrootznroot)
        SIGMA = muladd(POWER4, inv(ZNROOT * (ZN + LAMDA)), SIGMA)
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
    end

    EA = XNDEV * YNDEV
    EB = ZNDEV * ZNDEV
    EC = EA - EB
    ED = EA - sixT * EB
    EF = ED + EC + EC
    S1 = ED * (-C1 + C3 * ED * inv4 - threeT * C4 * ZNDEV * EF / twoT)
    S2 = ZNDEV * (C2 * EF + ZNDEV * (-C3 * EC + ZNDEV * C4 * EA))
    ans = threeT * SIGMA + POWER4 * (oneT + S1 + S2) / (MU * _sqrt(MU))

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRC
#***PURPOSE  Calculate a double precision approximation to
#             DRC(X,Y) = Integral from zero to infinity of
#                              -1/2     -1
#                    (1/2)(t+X)    (t+Y)  dt,
#            where X is nonnegative and Y is positive.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RC-S, DRC-D)
#***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
#             ELLIPTIC INTEGRAL, TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

@inline function DRC(X::A, Y::B) where {A,B}
    T = promote_type(A, B)

    ERRTOL = (eps(T) / 32)^T(1 / 6)
    LOLIM = T(5) * floatmin(T)
    UPLIM = T(floatmax(T) / 5)
    oneT = one(T)
    twoT = T(2)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)
    C1 = T(1 / 7)
    C2 = T(9 / 22)
    C3 = T(3 / 8)
    C4 = T(3 / 10)

    ans = zero(T)

    X < zero(T) || Y <= zero(T) && return (ans, 1)
    max(X, Y) > UPLIM && return (ans, 3)
    X + Y < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    MU = zero(T)
    SN = zero(T)

    while true
        MU = (XN + YN + YN) * inv3
        invMU = inv(MU)
        SN = muladd(invMU, YN + MU, -twoT)
        abs(SN) < ERRTOL && break
        LAMDA = muladd(twoT * _sqrt(XN), _sqrt(YN), YN)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
    end

    S = SN^2 * muladd(SN, muladd(SN, muladd(SN, C2, C3), C1), C4)
    ans = (oneT + S) / _sqrt(MU)

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRJ
#***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero)
#            elliptic integral of the 3rd kind.  For X, Y, and Z non-
#            negative, at most one of them zero, and P positive,
#             RJ(X,Y,Z,P) = Integral from zero to infinity of
#                              -1/2     -1/2     -1/2     -1
#                    (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RJ-S, DRJ-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRJ(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)

    ERRTOL = (eps(T) / 6)^T(1 / 6)
    LOLIM = (5floatmin(T))^T(1 / 3)
    UPLIM = T(3 / 10) * (floatmax(T) / 5)^T(1 / 3)
    twoT = T(2)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)

    C1 = T(3 / 14)
    C2 = T(1 / 3)
    C3 = T(3 / 22)
    C4 = T(3 / 26)

    ans = zero(T)

    min(X, Y, Z) < zero(T) && return (ans, 1)
    max(X, Y, Z, P) > UPLIM && return (ans, 3)
    min(X + Y, X + Z, Y + Z, P) < LOLIM && return (ans, 2)

    IER = 0
    XN = X
    YN = Y
    ZN = Z
    PN = P
    SIGMA = zero(T)
    POWER4 = one(T)
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)
    PNDEV = zero(T)

    while true
        XNYNZN = XN + YN + ZN
        MU = (XNYNZN + twoT * PN) * inv5
        invMU = inv(MU)
        XNDEV = (MU - XN) * invMU
        YNDEV = (MU - YN) * invMU
        ZNDEV = (MU - ZN) * invMU
        PNDEV = (MU - PN) * invMU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))
        EPSLON < ERRTOL && break
        XNROOT = _sqrt(XN)
        YNROOT = _sqrt(YN)
        ZNROOT = _sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        rootsum = XNROOT + YNROOT + ZNROOT
        LAMDA = muladd(XNROOT, (YNROOT + ZNROOT), YNROOTZNROOT)
        pn_plus_lamda = PN + LAMDA
        alpha_base = muladd(PN, rootsum, XNROOT * YNROOTZNROOT)
        ALFA = alpha_base * alpha_base
        BETA = PN * pn_plus_lamda * pn_plus_lamda
        drc, IER = DRC(ALFA, BETA)
        SIGMA = muladd(POWER4, drc, SIGMA)
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
        PN = pn_plus_lamda * inv4
    end
    YNDEVZNDEV = YNDEV * ZNDEV
    EA = muladd(XNDEV, (YNDEV + ZNDEV), YNDEVZNDEV)
    EB = XNDEV * YNDEVZNDEV
    EC = PNDEV * PNDEV
    E2 = EA - 3 * EC
    E3 = EB + 2 * PNDEV * (EA - EC)
    S1 = 1 + E2 * (-C1 + 3C3 / 4 * E2 - 3C4 / 2 * E3)
    S2 = EB * (C2 / 2 + PNDEV * (-C3 - C3 + PNDEV * C4))
    S3 = PNDEV * EA * (C2 - PNDEV * C3) - C2 * PNDEV * EC
    ans = 3SIGMA + POWER4 * (S1 + S2 + S3) / (MU * _sqrt(MU))

    return (ans, IER)
end

"""
    DRJ(X, Y, Z, Ps::NTuple{N}) -> NTuple{N, Tuple{T, Int}}

Batched Carlson `RJ` for several `P` sharing the same `(X, Y, Z)`. For valid input each element is
bit-identical to the scalar `DRJ(X, Y, Z, P)`, but the shared `(X, Y, Z)` duplication step (three
square roots, `LAMDA`) runs once per iteration for all `P`. Each `P` stops iterating at its own
convergence. Invalid input gives `(NaN, ier)` for the affected `P`, with the scalar method's error
code `ier`. Unrolled over `N`, allocation-free, GPU-compatible.
"""
@inline DRJ(X, Y, Z, ::Tuple{}) = ()
@inline function DRJ(X::A, Y::B, Z::C, Ps::NTuple{N,D}) where {A,B,C,D,N}
    T = promote_type(A, B, C, D)

    ERRTOL = (eps(T) / 6)^T(1 / 6)
    LOLIM = (5floatmin(T))^T(1 / 3)
    UPLIM = T(3 / 10) * (floatmax(T) / 5)^T(1 / 3)
    twoT = T(2)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)

    C1 = T(3 / 14)
    C2 = T(1 / 3)
    C3 = T(3 / 22)
    C4 = T(3 / 26)

    XN = T(X)
    YN = T(Y)
    ZN = T(Z)
    # closures below use `z` instead of capturing `T`: a captured DataType defeats inference
    z = zero(T)

    # Error codes of the scalar DRJ, per P; an invalid P is frozen from the start and returns NaN.
    # Inside the valid region DRC never errors (the SLATEC bounds are chosen for that).
    neg = min(XN, YN, ZN) < z
    maxXYZ = max(XN, YN, ZN)
    minsum = min(XN + YN, XN + ZN, YN + ZN)
    ier = map(Ps) do P
        PT = oftype(z, P)
        ifelse(neg, 1, ifelse(max(maxXYZ, PT) > UPLIM, 3, ifelse(min(minsum, PT) < LOLIM, 2, 0)))
    end

    # Per-P state (all-float, so it stays in registers). A converged P is frozen, `done`.
    st = map(P -> (PN=oftype(z, P), SIGMA=z, MU=z, XNDEV=z, YNDEV=z, ZNDEV=z, PNDEV=z, POWER4=z), Ps)
    done = map(!iszero, ier)
    POWER4 = one(T)

    while true
        XNYNZN = XN + YN + ZN
        # `let`: closures must not capture the loop-reassigned XN/YN/ZN/POWER4 (they would be boxed);
        # in-body `@inline`: without it the closures are not inlined for N ≥ 2 (stack copies, ~2x slower)
        st_done = let XN = XN, YN = YN, ZN = ZN, XNYNZN = XNYNZN, POWER4 = POWER4
            map(st, done) do s, d
                @inline
                d && return (s, d)
                MU = (XNYNZN + twoT * s.PN) * inv5
                invMU = inv(MU)
                XNDEV = (MU - XN) * invMU
                YNDEV = (MU - YN) * invMU
                ZNDEV = (MU - ZN) * invMU
                PNDEV = (MU - s.PN) * invMU
                EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))
                # `!(≥)` rather than `<`: a NaN input terminates (→ NaN) instead of looping forever
                ((; s..., MU, XNDEV, YNDEV, ZNDEV, PNDEV, POWER4), !(EPSLON >= ERRTOL))
            end
        end
        st = map(first, st_done)
        done = map(last, st_done)
        all(done) && break

        XNROOT = _sqrt(XN)
        YNROOT = _sqrt(YN)
        ZNROOT = _sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        rootsum = XNROOT + YNROOT + ZNROOT
        LAMDA = muladd(XNROOT, (YNROOT + ZNROOT), YNROOTZNROOT)
        st = let LAMDA = LAMDA, POWER4 = POWER4
            map(st, done) do s, d
                @inline
                d && return s
                pn_plus_lamda = s.PN + LAMDA
                alpha_base = muladd(s.PN, rootsum, XNROOT * YNROOTZNROOT)
                ALFA = alpha_base * alpha_base
                BETA = s.PN * pn_plus_lamda * pn_plus_lamda
                drc, _ = DRC(ALFA, BETA)
                (; s..., SIGMA=muladd(POWER4, drc, s.SIGMA), PN=pn_plus_lamda * inv4)
            end
        end
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
    end

    return map(st, ier) do s, e
        @inline
        (; MU, XNDEV, YNDEV, ZNDEV, PNDEV, SIGMA) = s
        YNDEVZNDEV = YNDEV * ZNDEV
        EA = muladd(XNDEV, (YNDEV + ZNDEV), YNDEVZNDEV)
        EB = XNDEV * YNDEVZNDEV
        EC = PNDEV * PNDEV
        E2 = EA - 3 * EC
        E3 = EB + 2 * PNDEV * (EA - EC)
        S1 = 1 + E2 * (-C1 + 3C3 / 4 * E2 - 3C4 / 2 * E3)
        S2 = EB * (C2 / 2 + PNDEV * (-C3 - C3 + PNDEV * C4))
        S3 = PNDEV * EA * (C2 - PNDEV * C3) - C2 * PNDEV * EC
        (ifelse(iszero(e), 3SIGMA + s.POWER4 * (S1 + S2 + S3) / (MU * _sqrt(MU)), oftype(z, NaN)), e)
    end
end

@inline function DRFJ(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)

    ERRTOLJ = (eps(T) / 6)^T(1 / 6)
    LOLIMJ = (5floatmin(T))^T(1 / 3)
    UPLIMJ = T(3 / 10) * (floatmax(T) / 5)^T(1 / 3)

    ERRTOL = (4 * eps(T) / 2)^T(1 / 6)
    LOLIM = 5floatmin(T)
    UPLIM = floatmax(T) / 5
    oneT = one(T)
    twoT = T(2)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
    c0 = T(1 / 10)

    ans = zero(T)
    ansJ = zero(T)

    ierr1 = 0
    ierr2 = 0
    min(X, Y, Z) < zero(T) && return (ans, 1, ansJ, 1)
    ierr1 = max(X, Y, Z) > UPLIM ? 3 : 0
    ierr2 = max(X, Y, Z, P) > UPLIMJ ? 3 : 0
    if ierr1 != 0 || ierr2 != 0
        return (ans, ierr1, ansJ, ierr2)
    end
    ierr1 = min(X + Y, X + Z, Y + Z) < LOLIM ? 2 : 0
    ierr2 = min(X + Y, X + Z, Y + Z, P) < LOLIMJ ? 2 : 0
    if ierr1 != 0 || ierr2 != 0
        return (ans, ierr1, ansJ, ierr2)
    end


    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)

    C1J = T(3 / 14)
    C2J = T(1 / 3)
    C3J = T(3 / 22)
    C4J = T(3 / 26)


    XN = X
    YN = Y
    ZN = Z
    PN = P
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    IER = 0
    SIGMA = zero(T)
    POWER4 = one(T)
    MUJ = zero(T)
    XNDEVJ = zero(T)
    YNDEVJ = zero(T)
    ZNDEVJ = zero(T)
    PNDEVJ = zero(T)

    fflag = true
    jflag = true
    while true
        XNYNZN = XN + YN + ZN
        if fflag
            MU = XNYNZN * inv3
            ninvMU = -inv(MU)
            XNDEV = muladd(ninvMU, (MU + XN), 2)
            YNDEV = muladd(ninvMU, (MU + YN), 2)
            ZNDEV = muladd(ninvMU, (MU + ZN), 2)
            EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
            fflag = EPSLON ≥ ERRTOL
        end

        if jflag
            MUJ = (XNYNZN + twoT * PN) * inv5
            invMUJ = inv(MUJ)
            XNDEVJ = (MUJ - XN) * invMUJ
            YNDEVJ = (MUJ - YN) * invMUJ
            ZNDEVJ = (MUJ - ZN) * invMUJ
            PNDEVJ = (MUJ - PN) * invMUJ
            EPSLONJ = max(abs(XNDEVJ), abs(YNDEVJ), abs(ZNDEVJ), abs(PNDEVJ))
            jflag = EPSLONJ ≥ ERRTOLJ
        end
        fflag || jflag || break

        XNROOT = _sqrt(XN)
        YNROOT = _sqrt(YN)
        ZNROOT = _sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        rootsum = XNROOT + YNROOT + ZNROOT
        LAMDA = muladd(XNROOT, (YNROOT + ZNROOT), YNROOTZNROOT)
        pn_plus_lamda = PN + LAMDA
        ALFA = muladd(PN, rootsum, XNROOT * YNROOTZNROOT)^2
        BETA = PN * pn_plus_lamda^2
        drc, IER = DRC(ALFA, BETA)
        SIGMA = muladd(POWER4, drc, SIGMA)
        POWER4 *= inv4

        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
        PN = pn_plus_lamda * inv4
    end
    XNDEVYNDEV = XNDEV * YNDEV
    E2 = muladd(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV * ZNDEV
    S = oneT + muladd(E2, muladd(-C2, E3, muladd(C1, E2, -c0)), C3 * E3)
    ans = S / _sqrt(MU)

    YNDEVZNDEVJ = YNDEVJ * ZNDEVJ
    EAJ = muladd(XNDEVJ, (YNDEVJ + ZNDEVJ), YNDEVZNDEVJ)
    EBJ = XNDEVJ * YNDEVZNDEVJ
    ECJ = PNDEVJ * PNDEVJ
    E2J = EAJ - 3 * ECJ
    E3J = muladd(twoT * PNDEVJ, EAJ - ECJ, EBJ)
    S1J = oneT + E2J * (-C1J + 3C3J / 4 * E2J - 3C4J / 2 * E3J)
    S2J = EBJ * (C2J / twoT + PNDEVJ * (-C3J - C3J + PNDEVJ * C4J))
    S3J = PNDEVJ * EAJ * (C2J - PNDEVJ * C3J) - C2J * PNDEVJ * ECJ
    ansJ = 3SIGMA + POWER4 * (S1J + S2J + S3J) / (MUJ * _sqrt(MUJ))


    return (ans, 0, ansJ, IER)
end


#end # module
