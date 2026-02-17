@testset "Low Thrust Cartesian Differentiability State" begin
    _lt_cartesian = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(1e-7, 2e-7, 3e-7)
    )

    for backend in _BACKENDS
        testname = "Low Thrust Cartesian Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _lt_cartesian), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _lt_cartesian)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Low Thrust Tangential Differentiability State" begin
    _lt_tangential = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))

    for backend in _BACKENDS
        testname = "Low Thrust Tangential Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _lt_tangential), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _lt_tangential)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-4
        end
    end
end

@testset "Low Thrust Tangential Differentiability Time" begin
    _lt_tangential = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))

    for backend in _BACKENDS
        testname = "Low Thrust Tangential Time Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _lt_tangential), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _lt_tangential)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Low Thrust RTN Differentiability State" begin
    _lt_rtn = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
    )

    for backend in _BACKENDS
        testname = "Low Thrust RTN Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _lt_rtn), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _lt_rtn)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-4
        end
    end
end

@testset "Low Thrust VNB Differentiability State" begin
    _lt_vnb = LowThrustAstroModel(;
        thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=VNBFrame()
    )

    for backend in _BACKENDS
        testname = "Low Thrust VNB Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _lt_vnb), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _lt_vnb)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-4
        end
    end
end

@testset "Low Thrust Piecewise Differentiability State" begin
    _lt_pw = LowThrustAstroModel(;
        thrust_model=PiecewiseConstantThrust(
            [0.0, 7200.0], [SVector{3}(0.0, 1e-7, 0.0), SVector{3}(1e-7, 0.0, 0.0)]
        ),
        frame=RTNFrame(),
    )

    for backend in _BACKENDS
        testname = "Low Thrust Piecewise State Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _lt_pw), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _lt_pw)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-4
        end
    end
end
