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

    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=ExpAtmo(),
        eop_data=eop_data,
    )

    @check_allocs dg_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test dg_accel(state, p, t, drag_model) isa SVector
end

@testset "Gravitational Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        eop_data=eop_data,
        order=36,
        degree=36,
        P=MMatrix{37,37,Float64}(zeros(37, 37)),
        dP=MMatrix{37,37,Float64}(zeros(37, 37)),
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    @check_allocs zon_accel(state, p, t, grav_model) = acceleration(state, p, t, grav_model)

    @test zon_accel(state, p, t, grav_model) isa SVector
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

    @check_allocs lense_thirr_accel(state, p, t, satellite_lense_thirring_model) = acceleration(
        state, p, t, satellite_lense_thirring_model
    )

    @test lense_thirr_accel(state, p, t, satellite_lense_thirring_model) isa SVector

    satellite_de_sitter_model = RelativityModel(;
        schwartzchild_effect=false, lense_thirring_effect=false, de_Sitter_effect=true
    )

    @check_allocs de_sitt_accel(state, p, t, satellite_de_sitter_model) = acceleration(
        state, p, t, satellite_de_sitter_model
    )
    @test de_sitt_accel(state, p, t, satellite_de_sitter_model) isa SVector

    satellite_schwartzchild_model = RelativityModel(;
        schwartzchild_effect=true, lense_thirring_effect=false, de_Sitter_effect=false
    )

    @check_allocs schwartz_accel(state, p, t, satellite_schwartzchild_model) = acceleration(
        state, p, t, satellite_schwartzchild_model
    )
    @test schwartz_accel(state, p, t, satellite_schwartzchild_model) isa SVector
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

    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model, sun_data=sun_model, eop_data=eop_data
    )
    @check_allocs sr_accel(state, p, t, srp_model) = acceleration(state, p, t, srp_model)

    @test sr_accel(state, p, t, srp_model) isa SVector
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

    @check_allocs sun_3body_accel(state, p, t, sun_third_body) = acceleration(
        state, p, t, sun_third_body
    )
    @check_allocs moon_3body_accel(state, p, t, moon_third_body) = acceleration(
        state, p, t, moon_third_body
    )

    @test sun_3body_accel(state, p, t, sun_third_body) isa SVector
    @test moon_3body_accel(state, p, t, moon_third_body) isa SVector
end

@testset "Dynamics Builder Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        eop_data=eop_data,
        order=36,
        degree=36,
        P=MMatrix{37,37,Float64}(zeros(37, 37)),
        dP=MMatrix{37,37,Float64}(zeros(37, 37)),
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model, sun_data=sun_third_body, eop_data=eop_data
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

    t = 0.0
    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    @check_allocs accel(u, p, t, models) = build_dynamics_model(u, p, t, models)
    @test accel(state, p, t, model_list) isa SVector
end
