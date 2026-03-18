@testset "Dynamics Builder" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, body_fixed_frame=:ITRF, propagation_frame=:ICRF, order=36, degree=36
    )
    sun_third_body = test_sun_model()
    moon_third_body = test_moon_model()

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model, sun_data=sun_third_body, R_Occulting=AstroForceModels.R_EARTH
    )

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=JB2008(),
        eop_data=eop_data,
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    lt_model = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))

    t = 0.0
    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model, lt_model)
    )

    total_accel = build_dynamics_model(state, p, t, model_list)

    total_accel_summed =
        acceleration(state, p, t, grav_model) +
        acceleration(state, p, t, moon_third_body) +
        acceleration(state, p, t, sun_third_body) +
        acceleration(state, p, t, srp_model) +
        acceleration(state, p, t, drag_model) +
        acceleration(state, p, t, lt_model)

    @test total_accel_summed ≈ total_accel rtol=1e-14
end
