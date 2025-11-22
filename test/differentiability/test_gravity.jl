@testset "Harmonics Differentiability State" begin
    for backend in _BACKENDS
        testname = "Gravity Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _grav_model), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _grav_model)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 2e-1

            f_fd2, df_fd2 = value_and_gradient(
                (x) -> potential(x, _p, _t, _grav_model), AutoFiniteDiff(), _state
            )

            f_ad2, df_ad2 = value_and_gradient(
                (x) -> potential(x, _p, _t, _grav_model), backend[2], _state
            )

            @test f_fd2 ≈ f_ad2
            @test df_fd2 ≈ df_ad2 rtol = 2e-1

            f_fd3, df_fd3 = value_and_gradient(
                (x) -> potential_time_derivative(x, _p, _t, _grav_model),
                AutoFiniteDiff(),
                _state,
            )

            f_ad3, df_ad3 = value_and_gradient(
                (x) -> potential_time_derivative(x, _p, _t, _grav_model),
                backend[2],
                _state,
            )

            @test f_fd3 ≈ f_ad3
            @test df_fd3 ≈ df_ad3 rtol = 2e-1
        end
    end
end

@testset "Harmonics Differentiability Time" begin
    for backend in _BACKENDS
        testname = "Gravity Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _grav_model), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _grav_model)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10

            f_fd2, df_fd2 = value_and_derivative(
                (x) -> potential(_state, _p, x, _grav_model), AutoFiniteDiff(), _t
            )

            f_ad2, df_ad2 = value_and_derivative(
                (x) -> potential(_state, _p, x, _grav_model), backend[2], _t
            )

            @test f_fd2 ≈ f_ad2
            @test df_fd2 ≈ df_ad2 atol = 1e-6

            f_fd3, df_fd3 = value_and_derivative(
                (x) -> potential_time_derivative(_state, _p, x, _grav_model),
                AutoFiniteDiff(),
                _t,
            )

            f_ad3, df_ad3 = value_and_derivative(
                (x) -> potential_time_derivative(_state, _p, x, _grav_model),
                backend[2],
                _t,
            )
        end
    end
end
