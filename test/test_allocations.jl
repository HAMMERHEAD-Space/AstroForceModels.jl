# AllocCheck reports spurious `jl_get_pgcstack_static` "allocating runtime
# call"s on macOS aarch64 with Julia 1.12+. These are not real heap allocations:
# the analyzed code is allocation-free on every other platform/version (Linux,
# Windows, and macOS on Julia 1.10/1.11). This is a known AllocCheck/Julia
# limitation, so the checks are skipped on the affected platform.
# Ref: https://github.com/SciML/SciMLStructures.jl/issues/59
const _SKIP_ALLOCCHECK = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION >= v"1.12"

if _SKIP_ALLOCCHECK
    @info "Skipping AllocCheck allocation tests (spurious jl_get_pgcstack_static reports on macOS aarch64 + Julia 1.12+; see SciML/SciMLStructures.jl#59)."
end

# Wrapper around `check_allocs` that honors the platform skip and, when real
# allocations are detected, dumps the full vector (with backtraces) to stdout so
# CI logs reveal exactly what is allocating.
function checked_allocs(f, types)
    _SKIP_ALLOCCHECK && return ()
    allocs = check_allocs(f, types)
    if !isempty(allocs)
        printstyled(stdout, "\n[ALLOC] "; color=:red, bold=true)
        println(stdout, f, " with ", types, " => ", length(allocs), " allocation(s)")
        for (i, a) in enumerate(allocs)
            println(stdout, "  ──────── allocation ", i, " ────────")
            show(stdout, MIME"text/plain"(), a)
            println(stdout)
        end
        flush(stdout)
    end
    return allocs
end

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

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(drag_model))
        ),
    ) == 0
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

    grav_types = (typeof(state), typeof(p), typeof(t), typeof(grav_model))

    @test length(checked_allocs(acceleration, grav_types)) == 0
    @test length(checked_allocs(potential, grav_types)) == 0
    @test length(checked_allocs(potential_time_derivative, grav_types)) == 0
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
        schwarzschild_effect=false, lense_thirring_effect=true, de_Sitter_effect=false
    )

    @test length(
        checked_allocs(
            acceleration,
            (typeof(state), typeof(p), typeof(t), typeof(satellite_lense_thirring_model)),
        ),
    ) == 0

    satellite_de_sitter_model = RelativityModel(;
        schwarzschild_effect=false, lense_thirring_effect=false, de_Sitter_effect=true
    )

    @test length(
        checked_allocs(
            acceleration,
            (typeof(state), typeof(p), typeof(t), typeof(satellite_de_sitter_model)),
        ),
    ) == 0

    satellite_schwarzschild_model = RelativityModel(;
        schwarzschild_effect=true, lense_thirring_effect=false, de_Sitter_effect=false
    )

    @test length(
        checked_allocs(
            acceleration,
            (typeof(state), typeof(p), typeof(t), typeof(satellite_schwarzschild_model)),
        ),
    ) == 0
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
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(srp_model))
        ),
    ) == 0
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

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(sun_third_body))
        ),
    ) == 0

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(moon_third_body))
        ),
    ) == 0
end

@testset "Low Thrust Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
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
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(cartesian_model))
        ),
    ) == 0

    tangential_model = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(tangential_model))
        ),
    ) == 0

    rtn_model = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
    )
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(rtn_model))
        ),
    ) == 0

    vnb_model = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=VNBFrame()
    )
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(vnb_model))
        ),
    ) == 0

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
    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(pw_model))
        ),
    ) == 0
end

@testset "Plasma Drag Allocations" begin
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

    satellite_plasma_drag_model = CannonballFixedPlasmaDrag(0.025)

    plasma_drag_model = PlasmaDragAstroModel(;
        satellite_plasma_drag_model=satellite_plasma_drag_model,
        ionosphere_model=ChapmanIonosphere(),
        eop_data=eop_data,
    )

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(plasma_drag_model))
        ),
    ) == 0

    plasma_drag_const = PlasmaDragAstroModel(;
        satellite_plasma_drag_model=satellite_plasma_drag_model,
        ionosphere_model=ConstantIonosphere(; rho_i=1e-17),
        eop_data=eop_data,
    )

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(plasma_drag_const))
        ),
    ) == 0
end

@testset "Solid Earth Tides Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    tides_model = SolidBodyTidesModel(eop_data)

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(tides_model))
        ),
    ) == 0
end

@testset "Thermal Emission Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    thermal_sat = FixedThermalEmission(0.01)
    thermal_model = ThermalEmissionAstroModel(;
        satellite_thermal_model=thermal_sat, sun_data=sun_model, eop_data=eop_data
    )

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(thermal_model))
        ),
    ) == 0
end

@testset "Magnetic Field Dipole Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
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
        eop_data=eop_data,
    )

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(mag_model))
        ),
    ) == 0
end

@testset "Albedo Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    satellite_shape_model = CannonballFixedSRP(0.2)

    albedo_model = AlbedoAstroModel(;
        satellite_shape_model=satellite_shape_model,
        sun_data=sun_third_body,
        body_albedo_model=UniformAlbedoModel(0.3, 0.7),
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

    @test length(
        checked_allocs(
            acceleration, (typeof(state), typeof(p), typeof(t), typeof(albedo_model))
        ),
    ) == 0
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

    @test length(
        checked_allocs(
            build_dynamics_model, (typeof(state), typeof(p), typeof(t), typeof(model_list))
        ),
    ) == 0
end
