@testset "Drag Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    t = 0.0

    SpaceIndices.init()
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

    satellite_drag_model = CannonballFixedDrag(0.2)

    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=ExpAtmo(),
        frames=p.frames,
    )

    @check_allocs dg_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test dg_accel(state, p, t, drag_model) isa SVector
end

@testset "Gravitational Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    t = 0.0

    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        body_fixed_frame=:ITRF,
        propagation_frame=:ICRF,
        order=36,
        degree=36,
        P=MMatrix{37,37,Float64}(zeros(37, 37)),
        dP=MMatrix{37,37,Float64}(zeros(37, 37)),
        frames=p.frames,
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

    @check_allocs pot_accel(state, p, t, grav_model) = potential(state, p, t, grav_model)

    @test pot_accel(state, p, t, grav_model) isa Number

    @check_allocs pot_time_accel(state, p, t, grav_model) = potential_time_derivative(
        state, p, t, grav_model
    )

    @test pot_time_accel(state, p, t, grav_model) isa Number
end

@testset "Relativity Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    t = 0.0
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

    earth_body_model = ThirdBodyModel(
        body=EarthBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=399, axes=:ICRF),
    )
    sun_model = test_sun_model(; frames=p.frames)

    satellite_lense_thirring_model = RelativityModel(;
        central_body=earth_body_model,
        sun_body=sun_model,
        J=SVector{3}(0.0, 0.0, AstroForceModels.EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS),
        schwarzschild_effect=false,
        lense_thirring_effect=true,
        de_Sitter_effect=false,
    )

    @check_allocs lense_thirr_accel(state, p, t, satellite_lense_thirring_model) = acceleration(
        state, p, t, satellite_lense_thirring_model
    )

    @test lense_thirr_accel(state, p, t, satellite_lense_thirring_model) isa SVector

    satellite_de_sitter_model = RelativityModel(;
        central_body=earth_body_model,
        sun_body=sun_model,
        schwarzschild_effect=false,
        lense_thirring_effect=false,
        de_Sitter_effect=true,
    )

    @check_allocs de_sitt_accel(state, p, t, satellite_de_sitter_model) = acceleration(
        state, p, t, satellite_de_sitter_model
    )
    @test de_sitt_accel(state, p, t, satellite_de_sitter_model) isa SVector

    satellite_schwarzschild_model = RelativityModel(;
        central_body=earth_body_model,
        sun_body=sun_model,
        schwarzschild_effect=true,
        lense_thirring_effect=false,
        de_Sitter_effect=false,
    )

    @check_allocs schwartz_accel(state, p, t, satellite_schwarzschild_model) = acceleration(
        state, p, t, satellite_schwarzschild_model
    )
    @test schwartz_accel(state, p, t, satellite_schwarzschild_model) isa SVector
end

@testset "SRP Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    t = 0.0
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
    sun_model = test_sun_model(; frames=p.frames)

    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model,
        sun_data=sun_model,
        R_Occulting=AstroForceModels.R_EARTH,
    )
    @check_allocs sr_accel(state, p, t, srp_model) = acceleration(state, p, t, srp_model)

    @test sr_accel(state, p, t, srp_model) isa SVector

    # Multi-occulter path must also be allocation-free. The tuple-recursive
    # `_resolve_occulters` / `_shadow_prod` helpers should specialize away.
    moon_tb = test_moon_model(; frames=p.frames)
    moon_occulter = OccultingBody(moon_tb, AstroForceModels.R_MOON)

    srp_multi = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model,
        sun_data=sun_model,
        R_Occulting=AstroForceModels.R_EARTH,
        additional_occulters=(moon_occulter,),
    )
    @check_allocs sr_multi_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test sr_multi_accel(state, p, t, srp_multi) isa SVector
end

@testset "Third Body Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    sun_third_body = test_sun_model(; frames=p.frames)
    moon_third_body = test_moon_model(; frames=p.frames)

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

@testset "Low Thrust Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    cartesian_model = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(1e-7, 2e-7, 3e-7)
    )
    @check_allocs cart_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test cart_accel(state, p, t, cartesian_model) isa SVector

    tangential_model = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))
    @check_allocs tang_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test tang_accel(state, p, t, tangential_model) isa SVector

    rtn_model = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
    )
    @check_allocs rtn_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test rtn_accel(state, p, t, rtn_model) isa SVector

    vnb_model = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=VNBFrame()
    )
    @check_allocs vnb_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test vnb_accel(state, p, t, vnb_model) isa SVector

    pw_model = LowThrustAstroModel(;
        thrust_model=PiecewiseConstantThrust(
            [0.0, 3600.0, 7200.0],
            [
                SVector{3}(1e-7, 0.0, 0.0),
                SVector{3}(0.0, 1e-7, 0.0),
                SVector{3}(-1e-7, 0.0, 0.0),
            ],
        ),
        frame=RTNFrame(),
    )
    @check_allocs pw_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test pw_accel(state, p, t, pw_model) isa SVector
end

@testset "Plasma Drag Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    t = 0.0
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

    satellite_plasma_drag_model = CannonballFixedPlasmaDrag(0.025)

    plasma_drag_model = PlasmaDragAstroModel(;
        satellite_plasma_drag_model=satellite_plasma_drag_model,
        ionosphere_model=ChapmanIonosphere(),
        frames=p.frames,
    )

    @check_allocs pd_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test pd_accel(state, p, t, plasma_drag_model) isa SVector

    plasma_drag_const = PlasmaDragAstroModel(;
        satellite_plasma_drag_model=satellite_plasma_drag_model,
        ionosphere_model=ConstantIonosphere(; rho_i=1e-17),
        frames=p.frames,
    )

    @check_allocs pd_const_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test pd_const_accel(state, p, t, plasma_drag_const) isa SVector
end

@testset "Solid Earth Tides Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    tides_model = SolidBodyTidesModel(
        tide_raising_bodies=(
            test_sun_model(; frames=p.frames), test_moon_model(; frames=p.frames)
        ),
        R_e=AstroForceModels.R_EARTH,
    )

    @check_allocs tides_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test tides_accel(state, p, t, tides_model) isa SVector
end

@testset "Thermal Emission Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_model = test_sun_model(; frames=p.frames)

    thermal_sat = FixedThermalEmission(0.01)
    thermal_model = ThermalEmissionAstroModel(;
        satellite_thermal_model=thermal_sat,
        sun_data=sun_model,
        R_Occulting=AstroForceModels.R_EARTH,
    )

    @check_allocs thm_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test thm_accel(state, p, t, thermal_model) isa SVector

    moon_tb = test_moon_model(; frames=p.frames)
    moon_occulter = OccultingBody(moon_tb, AstroForceModels.R_MOON)
    thermal_multi = ThermalEmissionAstroModel(;
        satellite_thermal_model=thermal_sat,
        sun_data=sun_model,
        R_Occulting=AstroForceModels.R_EARTH,
        additional_occulters=(moon_occulter,),
    )
    @check_allocs thm_multi_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test thm_multi_accel(state, p, t, thermal_multi) isa SVector
end

@testset "Magnetic Field Dipole Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    mag_model = MagneticFieldAstroModel(;
        spacecraft_charge_model=FixedChargeMassRatio(1e-3),
        geomagnetic_field_model=DipoleMagneticField(),
        frames=p.frames,
    )

    @check_allocs mag_accel(state, p, t, model) = acceleration(state, p, t, model)
    @test mag_accel(state, p, t, mag_model) isa SVector
end

@testset "Albedo Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    t = 0.0

    sun_third_body = test_sun_model(; frames=p.frames)

    satellite_shape_model = CannonballFixedSRP(0.2)

    albedo_model = AlbedoAstroModel(;
        satellite_shape_model=satellite_shape_model,
        sun_data=sun_third_body,
        body_albedo_model=UniformAlbedoModel(0.3, 0.7),
        body_fixed_frame=:ITRF,
        propagation_frame=:ICRF,
        frames=p.frames,
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    @check_allocs albedo_alloc_accel(state, p, t, albedo_model) = acceleration(
        state, p, t, albedo_model
    )

    @test albedo_alloc_accel(state, p, t, albedo_model) isa SVector
end

@testset "Dynamics Builder Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        body_fixed_frame=:ITRF,
        propagation_frame=:ICRF,
        order=36,
        degree=36,
        P=MMatrix{37,37,Float64}(zeros(37, 37)),
        dP=MMatrix{37,37,Float64}(zeros(37, 37)),
        frames=p.frames,
    )
    sun_third_body = test_sun_model(; frames=p.frames)
    moon_third_body = test_moon_model(; frames=p.frames)

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model,
        sun_data=sun_third_body,
        R_Occulting=AstroForceModels.R_EARTH,
    )

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=JB2008(),
        frames=p.frames,
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
