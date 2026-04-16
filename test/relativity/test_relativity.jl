@testset "Relativity Acceleration" begin
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

    sun_model = test_sun_model()
    earth_model = ThirdBodyModel(
        body=EarthBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=399, axes=:ICRF),
    )

    satellite_lense_thirring_model = RelativityModel(;
        central_body=earth_model,
        sun_body=sun_model,
        J=SVector{3}(0.0, 0.0, AstroForceModels.EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS),
        schwarzschild_effect=false,
        lense_thirring_effect=true,
        de_Sitter_effect=false,
    )

    lense_thirring_accel = acceleration(state, p, 0.0, satellite_lense_thirring_model)

    # Generated with Orekit
    expected_acceleration = [
        -9.318853003321375e-14, -2.495516635286037e-13, -3.817283343045848e-14
    ] # km/s

    @test lense_thirring_accel ≈ expected_acceleration atol = 1e-13

    satellite_de_sitter_model = RelativityModel(;
        central_body=earth_model,
        sun_body=sun_model,
        J=SVector{3}(0.0, 0.0, AstroForceModels.EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS),
        schwarzschild_effect=false,
        lense_thirring_effect=false,
        de_Sitter_effect=true,
    )

    de_sitter_accel = acceleration(state, p, 0.0, satellite_de_sitter_model)

    # Generated with Orekit
    expected_acceleration = [
        1.992497087381168e-14, 5.674171816594502e-14, -1.2737072865828191e-15
    ] # km/s

    #TODO: THIS ISN'T A GREAT TEST BUT THE SUN POSITION MODEL NEEDS TO BE IMPROVED FIRST
    @test de_sitter_accel ≈ expected_acceleration atol = 1E-13

    satellite_schwarzschild_model = RelativityModel(;
        central_body=earth_model,
        sun_body=sun_model,
        J=SVector{3}(0.0, 0.0, AstroForceModels.EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS),
        schwarzschild_effect=true,
        lense_thirring_effect=false,
        de_Sitter_effect=false,
    )

    schwarzschild_accel = acceleration(state, p, 0.0, satellite_schwarzschild_model)

    # Generated with Orekit
    expected_acceleration = [
        4.591078254044568e-12, -1.4642193807865207e-11, -1.4370450162013855e-12
    ] # km/s

    @test schwarzschild_accel ≈ expected_acceleration rtol = 1E-15
end
