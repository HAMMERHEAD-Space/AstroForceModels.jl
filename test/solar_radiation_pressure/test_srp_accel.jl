@testset "SRP Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_srp_model = CannonballFixedSRP(0.2)

    #TODO: RESOLVE SUN'S POSITION WITH HIGHER FIDELITY MODEL
    sun_model = test_sun_model()

    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model,
        sun_data=sun_model,
        R_Occulting=AstroForceModels.R_EARTH,
    )

    srp_accel = acceleration(state, p, 0.0, srp_model)

    # Generated with Orekit
    expected_acceleration = [
        -2.329584903106538e-10, 8.386410966633495e-10, 3.63558686733764e-10
    ] # km/s

    # This should resolve more after we replace the Sun's position with a higher fidelity
    @test srp_accel ≈ expected_acceleration atol = 1e-11

    @testset "additional_occulters: empty tuple equals single-body result" begin
        srp_with_empty = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=sun_model,
            R_Occulting=AstroForceModels.R_EARTH,
            additional_occulters=(),
        )
        @test acceleration(state, p, 0.0, srp_with_empty) ≈ srp_accel rtol = 1e-14
    end

    @testset "additional_occulters: body larger than Sun's disk eclipses it" begin
        # Place an OccultingBody at the Sun's position with a radius large enough to
        # fully cover the Sun's disk as seen from the spacecraft (≫ R_SUN), but
        # still smaller than the spacecraft→Sun distance so that the geometry is
        # physically valid (asin(R / d) remains in [-1, 1]).
        #
        # Full-eclipse criterion (Conical model): c < |b - a|, where
        #   a = asin(R_SUN / d_sun),   b = asin(R_occ / d_sun),   c = 0 (same direction).
        # Any R_occ > R_SUN yields b > a and c=0 → shadow factor = 0.
        sun_tb = test_sun_model(; frames=p.frames)
        big_R = 1.0e7  # km — well above R_SUN (695,700 km), far below 1 AU
        giant_occulter = OccultingBody(sun_tb, big_R)

        srp_eclipsed = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=sun_model,
            R_Occulting=AstroForceModels.R_EARTH,
            additional_occulters=(giant_occulter,),
        )

        a = acceleration(state, p, 0.0, srp_eclipsed)
        @test norm(a) ≈ 0.0 atol = 1e-15
    end

    @testset "additional_occulters: non-eclipsing body leaves result unchanged" begin
        # An occulter placed far behind the spacecraft (on the anti-solar side)
        # cannot occlude the Sun, so the SRP acceleration must match the
        # single-body baseline to full precision.
        sun_tb = test_sun_model(; frames=p.frames)

        # A very tiny occulter is also safe: even when present on the sunward side,
        # its angular size is negligible for all four shadow models.
        tiny = OccultingBody(sun_tb, 1e-6)

        srp_tiny = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=sun_model,
            R_Occulting=AstroForceModels.R_EARTH,
            additional_occulters=(tiny,),
        )
        @test acceleration(state, p, 0.0, srp_tiny) ≈ srp_accel rtol = 1e-12
    end
end
