@testset "Low Thrust Acceleration" begin
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

    @testset "LowThrustAstroModel with ConstantCartesianThrust" begin
        thrust_model = ConstantCartesianThrust(1e-7, 2e-7, 3e-7)
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        accel = acceleration(state, p, t, lt_model)

        @test accel isa SVector{3}
        @test accel ≈ SVector{3}(1e-7, 2e-7, 3e-7)
    end

    @testset "LowThrustAstroModel with ConstantTangentialThrust" begin
        thrust_model = ConstantTangentialThrust(1e-7)
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        accel = acceleration(state, p, t, lt_model)

        v = SVector{3}(state[4], state[5], state[6])
        v_hat = normalize(v)

        @test accel isa SVector{3}
        @test norm(accel) ≈ 1e-7
        @test normalize(accel) ≈ v_hat
    end

    @testset "LowThrustAstroModel with StateThrustModel" begin
        f = (u, p, t) -> SVector{3}(1e-7, 0.0, 0.0)
        thrust_model = StateThrustModel(f)
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        accel = acceleration(state, p, t, lt_model)

        @test accel isa SVector{3}
        @test accel ≈ SVector{3}(1e-7, 0.0, 0.0)
    end

    @testset "LowThrustAstroModel with PiecewiseConstantThrust" begin
        thrust_model = PiecewiseConstantThrust(
            [0.0, 3600.0], [SVector{3}(1e-7, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0)]
        )
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        accel = acceleration(state, p, 0.0, lt_model)
        @test accel ≈ SVector{3}(1e-7, 0.0, 0.0)

        accel_coast = acceleration(state, p, 3600.0, lt_model)
        @test accel_coast ≈ SVector{3}(0.0, 0.0, 0.0)
    end

    @testset "LowThrustAstroModel InertialFrame matches thrust_acceleration" begin
        thrust_model = ConstantTangentialThrust(0.5, 500.0)
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        accel_lt = acceleration(state, p, t, lt_model)
        accel_tm = thrust_acceleration(state, p, t, thrust_model)

        @test accel_lt ≈ accel_tm rtol = 1e-15
    end

    @testset "LowThrustAstroModel RTNFrame" begin
        r = SVector{3}(state[1], state[2], state[3])
        v = SVector{3}(state[4], state[5], state[6])
        R̂ = normalize(r)
        h = cross(r, v)
        N̂ = normalize(h)
        T̂ = cross(N̂, R̂)

        # Pure radial thrust in RTN → should be along R̂ in inertial
        lt_radial = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=RTNFrame()
        )
        accel_R = acceleration(state, p, t, lt_radial)
        @test accel_R ≈ 1e-7 * R̂

        # Pure transverse thrust in RTN → should be along T̂ in inertial
        lt_transverse = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
        )
        accel_T = acceleration(state, p, t, lt_transverse)
        @test accel_T ≈ 1e-7 * T̂

        # Pure normal thrust in RTN → should be along N̂ in inertial
        lt_normal = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(0.0, 0.0, 1e-7), frame=RTNFrame()
        )
        accel_N = acceleration(state, p, t, lt_normal)
        @test accel_N ≈ 1e-7 * N̂

        # Combined RTN thrust: verify linearity
        lt_combined = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(1e-7, 2e-7, 3e-7), frame=RTNFrame()
        )
        accel_combined = acceleration(state, p, t, lt_combined)
        expected = 1e-7 * R̂ + 2e-7 * T̂ + 3e-7 * N̂
        @test accel_combined ≈ expected
    end

    @testset "LowThrustAstroModel VNBFrame" begin
        r = SVector{3}(state[1], state[2], state[3])
        v = SVector{3}(state[4], state[5], state[6])
        V̂ = normalize(v)
        h = cross(r, v)
        B̂ = normalize(h)
        N̂ = cross(B̂, V̂)

        # Pure velocity thrust in VNB → should be along V̂ in inertial
        lt_velocity = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=VNBFrame()
        )
        accel_V = acceleration(state, p, t, lt_velocity)
        @test accel_V ≈ 1e-7 * V̂

        # Pure normal thrust in VNB → should be along N̂ in inertial
        lt_normal = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=VNBFrame()
        )
        accel_N = acceleration(state, p, t, lt_normal)
        @test accel_N ≈ 1e-7 * N̂

        # Pure binormal thrust in VNB → should be along B̂ in inertial
        lt_binormal = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(0.0, 0.0, 1e-7), frame=VNBFrame()
        )
        accel_B = acceleration(state, p, t, lt_binormal)
        @test accel_B ≈ 1e-7 * B̂

        # Combined VNB thrust: verify linearity
        lt_combined = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(1e-7, 2e-7, 3e-7), frame=VNBFrame()
        )
        accel_combined = acceleration(state, p, t, lt_combined)
        expected = 1e-7 * V̂ + 2e-7 * N̂ + 3e-7 * B̂
        @test accel_combined ≈ expected
    end

    @testset "VNBFrame velocity-aligned matches ConstantTangentialThrust" begin
        # VNBFrame with [magnitude, 0, 0] should give the same result as
        # ConstantTangentialThrust with the same magnitude
        mag = 1e-7
        lt_vnb = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(mag, 0.0, 0.0), frame=VNBFrame()
        )
        lt_tang = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(mag))

        accel_vnb = acceleration(state, p, t, lt_vnb)
        accel_tang = acceleration(state, p, t, lt_tang)

        @test accel_vnb ≈ accel_tang rtol = 1e-14
    end

    @testset "LowThrustAstroModel in CentralBodyDynamicsModel" begin
        # Verify low thrust can be composed with gravity in the dynamics builder
        keplerian = KeplerianGravityAstroModel()
        thrust_model = ConstantCartesianThrust(1e-7, 0.0, 0.0)
        lt_model = LowThrustAstroModel(; thrust_model=thrust_model)

        dynamics = CentralBodyDynamicsModel(keplerian, (lt_model,))
        total_accel = build_dynamics_model(state, p, t, dynamics)

        grav_accel = acceleration(state, p, t, keplerian)
        lt_accel = acceleration(state, p, t, lt_model)

        @test total_accel ≈ grav_accel + lt_accel rtol = 1e-14
    end

    @testset "LowThrustAstroModel RTN in CentralBodyDynamicsModel" begin
        keplerian = KeplerianGravityAstroModel()
        lt_rtn = LowThrustAstroModel(;
            thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
        )

        dynamics = CentralBodyDynamicsModel(keplerian, (lt_rtn,))
        total_accel = build_dynamics_model(state, p, t, dynamics)

        grav_accel = acceleration(state, p, t, keplerian)
        lt_accel = acceleration(state, p, t, lt_rtn)

        @test total_accel ≈ grav_accel + lt_accel rtol = 1e-14
    end
end
