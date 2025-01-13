@testset "Dynamics Builder" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model, sun_data=sun_third_body, eop_data=eop_data
    )

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(satellite_drag_model, JB2008(), eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    t = 0.0
    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    total_accel = build_dynamics_model(state, p, t, model_list)

    total_accel_summed =
        acceleration(state, p, t, grav_model) +
        acceleration(state, p, t, moon_third_body) +
        acceleration(state, p, t, sun_third_body) +
        acceleration(state, p, t, srp_model) +
        acceleration(state, p, t, drag_model)

    @test total_accel_summed == total_accel
end
