@testset "Drag Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_drag_model = CannonballFixedDrag(0.2)

    drag_model = DragAstroModel(satellite_drag_model, ExpAtmo(), eop_data)

    @check_allocs drag_accel(state, p, t, drag_model) = acceleration(
        state, p, t, drag_model
    )
    drag_accel(state, p, t, drag_model)
end

@testset "Gravitational Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    @check_allocs zonal_accel(state, p, t, grav_model) = acceleration(
        state, p, t, grav_model
    )
    zonal_accel(state, p, t, grav_model)
end

@testset "Relativity Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_lense_thirring_model = RelativityModel(;
        schwartzchild_effect=false, lense_thirring_effect=true, de_Sitter_effect=false
    )

    @check_allocs lense_thirring_accel(state, p, t, satellite_lense_thirring_model) = acceleration(
        state, p, t, satellite_lense_thirring_model
    )

    lense_thirring_accel(state, p, t, satellite_lense_thirring_model)

    satellite_de_sitter_model = RelativityModel(;
        schwartzchild_effect=false, lense_thirring_effect=false, de_Sitter_effect=true
    )

    @check_allocs de_sitter_accel(state, p, t, satellite_de_sitter_model) = acceleration(
        state, p, t, satellite_de_sitter_model
    )
    de_sitter_accel(state, p, t, satellite_de_sitter_model)

    satellite_schwartzchild_model = RelativityModel(;
        schwartzchild_effect=true, lense_thirring_effect=false, de_Sitter_effect=false
    )

    @check_allocs schwartzchild_accel(state, p, t, satellite_schwartzchild_model) = acceleration(
        state, p, t, satellite_schwartzchild_model
    )
    schwartzchild_accel(state, p, t, satellite_schwartzchild_model)
end

@testset "SRP Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0
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

    srp_model = SRPAstroModel(satellite_srp_model, sun_model, eop_data, Conical())

    @check_allocs srp_accel(state, p, t, srp_model) = acceleration(state, p, t, srp_model)

    srp_accel(state, p, t, srp_model)
end

@testset "Third Body Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    @check_allocs sun_accel(state, p, t, sun_third_body) = acceleration(
        state, p, t, sun_third_body
    )
    @check_allocs moon_accel(state, p, t, moon_third_body) = acceleration(
        state, p, t, moon_third_body
    )

    sun_accel(state, p, t, sun_third_body)
    moon_accel(state, p, t, moon_third_body)
end

@testset "Dynamics Builder Allocations" begin
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
    srp_model = SRPAstroModel(satellite_srp_model, sun_third_body, eop_data, Conical())

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

    @check_allocs accel(u, p, t, models) = build_dynamics_model(u, p, t, models)
    accel(state, p, t, model_list)
end
