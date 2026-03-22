@testset "setup_inertial_frames" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)

    @testset "default (Sun + Moon included)" begin
        p = setup_inertial_frames(epoch)

        @test p isa FrameAwareParams
        @test p.propagation_frame === :ICRF
        @test p.epoch === epoch
        @test p.JD ≈ JD atol = 1e-6

        # ICRF axes, Earth, Sun, Moon must all be registered
        frames = p.frames
        @test has_axes(frames, :ICRF)
        @test haskey(points_alias(frames), :Earth)
        @test haskey(points_alias(frames), :Sun)
        @test haskey(points_alias(frames), :Moon)
    end

    @testset "include_sun=false, include_moon=false" begin
        p = setup_inertial_frames(epoch; include_sun=false, include_moon=false)

        @test p isa FrameAwareParams
        @test has_axes(p.frames, :ICRF)
        @test haskey(points_alias(p.frames), :Earth)
        @test !haskey(points_alias(p.frames), :Sun)
        @test !haskey(points_alias(p.frames), :Moon)
    end

    @testset "include_sun=true, include_moon=false" begin
        p = setup_inertial_frames(epoch; include_sun=true, include_moon=false)
        @test haskey(points_alias(p.frames), :Sun)
        @test !haskey(points_alias(p.frames), :Moon)
    end

    @testset "custom order and numtype" begin
        p = setup_inertial_frames(epoch; order=1, numtype=Float32)
        @test p.frames isa FrameSystem{1,Float32}
    end

    @testset "Sun/Moon ephemeris callable" begin
        p = setup_inertial_frames(epoch)
        t_ft = ft_time(p, 0.0)

        sun_pos = get_position(
            FrameEphemeris(; center_point=399, target_point=10, axes=:ICRF),
            SunBody(),
            p.frames,
            t_ft,
        )
        moon_pos = get_position(
            FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF),
            MoonBody(),
            p.frames,
            t_ft,
        )

        # Sanity check: Sun ~1 AU from Earth, Moon ~384,400 km
        @test norm(sun_pos) > 1e8       # > 10^8 km
        @test norm(moon_pos) > 3.5e5    # > 350,000 km
    end

    @testset "usable with KeplerianGravityAstroModel" begin
        p = setup_inertial_frames(epoch; include_sun=false, include_moon=false)
        grav = KeplerianGravityAstroModel(; μ=AstroForceModels.μ_EARTH)
        dynamics = CentralBodyDynamicsModel(grav, ())

        state = [
            -1076.225324679696
            -6765.896364327722
            -332.3087833503755
            9.356857417032581
            -3.3123476319597557
            -1.1880157328553503
        ]

        accel = build_dynamics_model(state, p, 0.0, dynamics)
        @test accel isa SVector{3}
        @test norm(accel) > 0
    end

    @testset "usable with ThirdBodyModel" begin
        p = setup_inertial_frames(epoch)
        grav = KeplerianGravityAstroModel(; μ=AstroForceModels.μ_EARTH)
        sun = ThirdBodyModel(;
            body=SunBody(),
            ephem_type=FrameEphemeris(; center_point=399, target_point=10, axes=:ICRF),
            frames=p.frames,
        )
        moon = ThirdBodyModel(;
            body=MoonBody(),
            ephem_type=FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF),
            frames=p.frames,
        )
        dynamics = CentralBodyDynamicsModel(grav, (sun, moon))

        state = [
            -1076.225324679696
            -6765.896364327722
            -332.3087833503755
            9.356857417032581
            -3.3123476319597557
            -1.1880157328553503
        ]

        accel = build_dynamics_model(state, p, 0.0, dynamics)
        @test accel isa SVector{3}
        @test norm(accel) > 0
    end
end
