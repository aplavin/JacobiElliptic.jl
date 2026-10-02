#
# Tests for JacobiElliptic Functions
#

@testset "JacobiElliptic" begin
    @testset "Abramowitz & Stegun Table 16.1" begin
        dataloc = "data/"
        # table 16.1
        t161, _ = readdlm(joinpath(dataloc, "table_16_1.csv"), ',', header = true)
        # table 17.2
        t172, _ = readdlm(joinpath(dataloc, "table_17_2.csv"), ',', header = true)

        K_lut = Dict(zip(t172[:, 1], t172[:, 2]))

        # data vary first by ϵ ∈ 0:5:90 then α ∈ 0:5:85
        αs = 0:5:85
        ϵs = 0:5:90
        θss = reshape(t161[:, 3], length(ϵs), length(αs))
        θns = reshape(t161[:, 4], length(ϵs), length(αs))

        @testset "α = $α" for (i, α) in enumerate(αs)
            K = K_lut[α]
            m = sind(α)^2
            denom = sqrt(secd(α))

            @testset "ϵ = $ϵ" for (j, ϵ) in enumerate(ϵs)
                j₁ = length(ϵs) - j + 1
                ϵ₁ = ϵs[j₁]
                u = ϵ * K / 90

                θs = θss[j, i]
                θn = θns[j, i]
                θc = θss[j₁, i] / denom
                θd = θns[j₁, i] / denom

                @test JacobiElliptic.sn(u, m) ≈ θs / θn atol = 2.5e-9
                @test JacobiElliptic.cn(u, m) ≈ θc / θn atol = 2.5e-9
                @test JacobiElliptic.dn(u, m) ≈ θd / θn atol = 2.5e-9
                @test JacobiElliptic.nn(u, m) == 1.0

                @test JacobiElliptic.sd(u, m) ≈ θs / θd atol = 2.5e-9
                @test JacobiElliptic.cd(u, m) ≈ θc / θd atol = 2.5e-9
                @test JacobiElliptic.dd(u, m) == 1.0
                @test JacobiElliptic.nd(u, m) ≈ θn / θd atol = 2.5e-9

                @test JacobiElliptic.cc(u, m) == 1.0
                if ϵ != 90
                    # very sensitive around u = K,
                    # estimate of K(0) = pi/2 + 4e-9, so cosine causes errors
                    # also, errors build up in ϵ ≥75°, so lower tolerences for that
                    @test JacobiElliptic.sc(u, m) ≈ θs / θc atol = 1e-7
                    @test JacobiElliptic.dc(u, m) ≈ θd / θc atol = 1e-7
                    @test JacobiElliptic.nc(u, m) ≈ θn / θc atol = 1e-7
                end

                @test JacobiElliptic.ss(u, m) == 1.0
                @test JacobiElliptic.cs(u, m) ≈ θc / θs atol = 5.5e-8
                @test JacobiElliptic.ds(u, m) ≈ θd / θs atol = 5.5e-8
                @test JacobiElliptic.ns(u, m) ≈ θn / θs atol = 5.5e-8

                # ellipj
                s, c, d = JacobiElliptic.CarlsonAlg.ellipj(u, m)
                @test s ≈ θs / θn atol = 1e-9
                @test c ≈ θc / θn atol = 1e-9
                @test d ≈ θd / θn atol = 1e-9
            end
        end
    end

    @testset "u = $u" for u = -1.0:0.21:1.0
        @test JacobiElliptic.am(u, 0.0) ≈ u
        @test JacobiElliptic.sn(u, 0.0) ≈ sin(u)
        @test JacobiElliptic.cn(u, 0.0) ≈ cos(u)
        @test JacobiElliptic.dn(u, 0.0) ≈ 1.0
        @test JacobiElliptic.nn(u, 0.0) ≈ 1.0
        @test JacobiElliptic.cd(u, 0.0) ≈ cos(u)
        @test JacobiElliptic.sd(u, 0.0) ≈ sin(u)
        @test JacobiElliptic.nd(u, 0.0) ≈ 1.0
        @test JacobiElliptic.dd(u, 0.0) ≈ 1.0
        @test JacobiElliptic.dc(u, 0.0) ≈ sec(u)
        @test JacobiElliptic.nc(u, 0.0) ≈ sec(u)
        @test JacobiElliptic.sc(u, 0.0) ≈ tan(u)
        @test JacobiElliptic.cc(u, 0.0) ≈ 1.0
        @test JacobiElliptic.ns(u, 0.0) ≈ csc(u)
        @test JacobiElliptic.ds(u, 0.0) ≈ csc(u)
        @test JacobiElliptic.cs(u, 0.0) ≈ cot(u)
        @test JacobiElliptic.ss(u, 0.0) ≈ 1.0

        @test JacobiElliptic.am(u, 1.0) ≈ atan(sinh(u))
        @test JacobiElliptic.sn(u, 1.0) ≈ tanh(u)
        @test JacobiElliptic.cn(u, 1.0) ≈ sech(u)
        @test JacobiElliptic.dn(u, 1.0) ≈ sech(u)
        @test JacobiElliptic.nn(u, 1.0) ≈ 1.0
        @test JacobiElliptic.cd(u, 1.0) ≈ 1.0
        @test JacobiElliptic.sd(u, 1.0) ≈ sinh(u)
        @test JacobiElliptic.nd(u, 1.0) ≈ cosh(u)
        @test JacobiElliptic.dd(u, 1.0) ≈ 1.0
        @test JacobiElliptic.dc(u, 1.0) ≈ 1.0
        @test JacobiElliptic.nc(u, 1.0) ≈ cosh(u)
        @test JacobiElliptic.sc(u, 1.0) ≈ sinh(u)
        @test JacobiElliptic.cc(u, 1.0) ≈ 1.0
        @test JacobiElliptic.ns(u, 1.0) ≈ coth(u)
        @test JacobiElliptic.ds(u, 1.0) ≈ csch(u)
        @test JacobiElliptic.cs(u, 1.0) ≈ csch(u)
        @test JacobiElliptic.ss(u, 1.0) ≈ 1.0
    end

    @testset "negative modulus for Jacobi am" begin
        @test JacobiElliptic.am(0.5, -10.0) ≈ 0.717286208837285
        @test JacobiElliptic.am(10.0, -10.0) ≈ 19.58906782107541
    end

    @testset "modulus > 1 for Jacobi am" begin
        @test JacobiElliptic.am(0.5, 1.1) ≈ 0.4784409095408496
        @test JacobiElliptic.am(10.0, 1.1) ≈ -0.0038260319603258524
        @test JacobiElliptic.am(10.0, 10.0) ≈ -0.18538447091861052
    end

    @testset "am near m = 1 vs 256-bit reference" begin
        # am(u|m) from bisection + Newton on ArbNumerics.elliptic_f
        for ((u, m), ref) in (
            (20.0, 1 - 1e-10) => 1.576860869527043,
            (15.0, 1 - 1e-9) => 1.571204342143267,
            (12.0, 1 - 1e-8) => 1.5709874818544491,
            (25.0, 0.999) => 7.8818237562111815,
            (3.0, 0.3) => 2.7175782549216283,
        )
            @test JacobiElliptic.am(u, m) ≈ ref rtol = 4u * eps()
        end
        for ((u, m), ref) in ((10f0, 1 - 1f-5) => 1.5982737166354217, (6f0, 1 - 1f-4) => 1.5708812708993438)
            @test JacobiElliptic.am(u, m) ≈ ref rtol = 4u * eps(1f0)
        end
        @test JacobiElliptic.am(big"0.5", big"0.99") isa BigFloat
    end

    @testset "ellipj vs 256-bit reference" begin
        # (u, m) => (sn, cn, dn): inverse of ArbNumerics.elliptic_f by Newton; m > 1 via A&S 16.11
        refs = (
            (10.0, 1 - 1e-5) => (0.999623546614056, -0.027436564193360935, 0.027618065230675056),
            (5.0, 1 - 1e-12) => (0.9999092042628449, 0.013475282202772691, 0.013475282239870108),
            (12.0, 1 - 1e-8) => (0.9999999817298717, -0.00019115505838839932, 0.0002157319077751768),
            (25.0, 0.999) => (0.999612433151943, -0.027838525249234494, 0.04212135449583993),
            (-7.0, 0.9) => (0.9686123629830704, -0.24857612571675636, 0.39447570425657064),
            (3.0, 0.3) => (0.4114226606340007, -0.9114446743038438, 0.9742789222266137),
            (1.2, -3.0) => (0.9709312730622157, -0.23935844039804655, 1.9565588698096064),
            (0.8, 4.0) => (0.4986210795513635, 0.866820061504712, 0.07421641400749103),
        )
        for ((u, m), ref) in refs
            @test all(map((x, r) -> isapprox(x, r; rtol = 1e-13), ellipj(u, m), ref))
            @test (sn(u, m), cn(u, m), dn(u, m)) === ellipj(u, m)
        end
        @test all(map((x, r) -> isapprox(x, r; rtol = 32eps(1f0)), ellipj(10f0, 1 - 1f-5),
                      (0.9996225202746573, -0.027473932367654354, 0.027655433130660895)))
        @inferred ellipj(0.3f0, 0.5f0)
        @inferred ellipj(1, 0.5)
    end
end
