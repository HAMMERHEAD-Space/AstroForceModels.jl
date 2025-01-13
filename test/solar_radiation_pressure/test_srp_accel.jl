@testset "SRP Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()

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
    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model, sun_data=sun_model, eop_data=eop_data
    )

    srp_accel = acceleration(state, p, 0.0, srp_model)

    # Generated with Orekit
    expected_acceleration = [
        -2.329584903106538e-10, 8.386410966633495e-10, 3.63558686733764e-10
    ] # km/s

    # This should resolve more after we replace the Sun's position with a higher fidelity
    @test srp_accel ≈ expected_acceleration atol = 1e-11
end
