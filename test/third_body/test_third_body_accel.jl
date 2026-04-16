@testset "Third Body Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)

    sun_third_body = test_sun_model()
    moon_third_body = test_moon_model()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_accel = acceleration(state, p, 0.0, sun_third_body)
    moon_accel = acceleration(state, p, 0.0, moon_third_body)

    expected_solar_accel =
        [2.2643603803606135e-07, -3.715524933257638e-07, -2.6950132428989296e-07] ./ 1E3
    expected_lunar_accel =
        [-7.507182245381024e-07, 1.2454140260475816e-07, -1.7206122388643457e-07] ./ 1E3

    # Slightly wider tolerance due to Vallado analytical ephemeris (no EOP) vs reference (with EOP)
    @test sun_accel ≈ expected_solar_accel rtol = 5E-3
    @test moon_accel ≈ expected_lunar_accel rtol = 5E-3
end
