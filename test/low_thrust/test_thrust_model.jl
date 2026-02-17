@testset "Thrust Model" begin
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

    @testset "ConstantCartesianThrust" begin
        # Direct acceleration specification
        model = ConstantCartesianThrust(1e-7, 2e-7, 3e-7)
        accel = thrust_acceleration(state, p, t, model)

        @test accel ≈ SVector{3}(1e-7, 2e-7, 3e-7)
        @test accel isa SVector{3}

        # From physical parameters: direction, thrust [N], mass [kg]
        # 1 N thrust, 1000 kg mass → a = 1/(1000*1000) = 1e-6 km/s²
        direction = [1.0, 0.0, 0.0]
        model_phys = ConstantCartesianThrust(direction, 1.0, 1000.0)
        accel_phys = thrust_acceleration(state, p, t, model_phys)

        @test accel_phys ≈ SVector{3}(1e-6, 0.0, 0.0)

        # Diagonal direction should be normalized
        direction_diag = [1.0, 1.0, 1.0]
        model_diag = ConstantCartesianThrust(direction_diag, 1.0, 1000.0)
        accel_diag = thrust_acceleration(state, p, t, model_diag)

        expected_mag = 1.0 / (1000.0 * 1000.0)
        @test norm(accel_diag) ≈ expected_mag
        @test accel_diag[1] ≈ accel_diag[2]
        @test accel_diag[2] ≈ accel_diag[3]
    end

    @testset "ConstantCartesianThrust Validation" begin
        @test_throws ArgumentError ConstantCartesianThrust([1.0, 0.0, 0.0], -1.0, 100.0)
        @test_throws ArgumentError ConstantCartesianThrust([1.0, 0.0, 0.0], 1.0, 0.0)
        @test_throws ArgumentError ConstantCartesianThrust([1.0, 0.0, 0.0], 1.0, -1.0)
    end

    @testset "ConstantTangentialThrust" begin
        # Direct acceleration specification (orbit raising)
        model = ConstantTangentialThrust(1e-7)
        accel = thrust_acceleration(state, p, t, model)

        v = SVector{3}(state[4], state[5], state[6])
        v_hat = normalize(v)

        @test norm(accel) ≈ 1e-7
        @test normalize(accel) ≈ v_hat

        # Negative magnitude (braking)
        model_brake = ConstantTangentialThrust(-1e-7)
        accel_brake = thrust_acceleration(state, p, t, model_brake)

        @test norm(accel_brake) ≈ 1e-7
        @test normalize(accel_brake) ≈ -v_hat

        # From physical parameters: thrust [N], mass [kg]
        model_phys = ConstantTangentialThrust(1.0, 1000.0)
        accel_phys = thrust_acceleration(state, p, t, model_phys)

        @test norm(accel_phys) ≈ 1.0 / (1000.0 * 1000.0)
        @test normalize(accel_phys) ≈ v_hat
    end

    @testset "ConstantTangentialThrust Validation" begin
        @test_throws ArgumentError ConstantTangentialThrust(1.0, 0.0)
        @test_throws ArgumentError ConstantTangentialThrust(1.0, -1.0)
    end

    @testset "StateThrustModel" begin
        # Simple constant function returning fixed acceleration
        f_const = (u, p, t) -> SVector{3}(1e-7, 2e-7, 3e-7)
        model = StateThrustModel(f_const)
        accel = thrust_acceleration(state, p, t, model)

        @test accel ≈ SVector{3}(1e-7, 2e-7, 3e-7)

        # State-dependent: radial thrust
        f_radial = (u, p, t) -> begin
            r = SVector{3}(u[1], u[2], u[3])
            return 1e-7 * normalize(r)
        end
        model_radial = StateThrustModel(f_radial)
        accel_radial = thrust_acceleration(state, p, t, model_radial)

        r = SVector{3}(state[1], state[2], state[3])
        expected = 1e-7 * normalize(r)
        @test accel_radial ≈ expected

        # Time-dependent: thrust only during first hour
        f_scheduled = (u, p, t) -> begin
            if t < 3600.0
                return SVector{3}(1e-7, 0.0, 0.0)
            else
                return SVector{3}(0.0, 0.0, 0.0)
            end
        end
        model_sched = StateThrustModel(f_scheduled)

        accel_on = thrust_acceleration(state, p, 0.0, model_sched)
        @test accel_on ≈ SVector{3}(1e-7, 0.0, 0.0)

        accel_off = thrust_acceleration(state, p, 7200.0, model_sched)
        @test accel_off ≈ SVector{3}(0.0, 0.0, 0.0)
    end

    @testset "PiecewiseConstantThrust" begin
        # Single arc
        model_single = PiecewiseConstantThrust([0.0], [SVector{3}(1e-7, 0.0, 0.0)])
        accel = thrust_acceleration(state, p, 100.0, model_single)
        @test accel ≈ SVector{3}(1e-7, 0.0, 0.0)

        # Two arcs: thrust then coast
        model_two = PiecewiseConstantThrust(
            [0.0, 3600.0], [SVector{3}(1e-7, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0)]
        )
        accel_on = thrust_acceleration(state, p, 0.0, model_two)
        @test accel_on ≈ SVector{3}(1e-7, 0.0, 0.0)

        accel_mid = thrust_acceleration(state, p, 1800.0, model_two)
        @test accel_mid ≈ SVector{3}(1e-7, 0.0, 0.0)

        accel_off = thrust_acceleration(state, p, 3600.0, model_two)
        @test accel_off ≈ SVector{3}(0.0, 0.0, 0.0)

        accel_off2 = thrust_acceleration(state, p, 7200.0, model_two)
        @test accel_off2 ≈ SVector{3}(0.0, 0.0, 0.0)

        # Three arcs: direction changes 
        model_three = PiecewiseConstantThrust(
            [0.0, 3600.0, 7200.0],
            [
                SVector{3}(1e-7, 0.0, 0.0),
                SVector{3}(0.0, 1e-7, 0.0),
                SVector{3}(0.0, 0.0, 1e-7),
            ],
        )
        @test thrust_acceleration(state, p, 0.0, model_three) ≈ SVector{3}(1e-7, 0.0, 0.0)
        @test thrust_acceleration(state, p, 3600.0, model_three) ≈
            SVector{3}(0.0, 1e-7, 0.0)
        @test thrust_acceleration(state, p, 7200.0, model_three) ≈
            SVector{3}(0.0, 0.0, 1e-7)
        @test thrust_acceleration(state, p, 5000.0, model_three) ≈
            SVector{3}(0.0, 1e-7, 0.0)

        # Before first arc returns first arc
        @test thrust_acceleration(state, p, -100.0, model_three) ≈
            SVector{3}(1e-7, 0.0, 0.0)

        # Constructor validation: mismatched lengths
        @test_throws ArgumentError PiecewiseConstantThrust(
            [0.0, 3600.0], [SVector{3}(1e-7, 0.0, 0.0)]
        )

        # Constructor validation: non-increasing times
        @test_throws ArgumentError PiecewiseConstantThrust(
            [3600.0, 0.0], [SVector{3}(1e-7, 0.0, 0.0), SVector{3}(0.0, 1e-7, 0.0)]
        )

        # Constructor validation: empty
        @test_throws ArgumentError PiecewiseConstantThrust(Float64[], SVector{3,Float64}[])
    end
end
