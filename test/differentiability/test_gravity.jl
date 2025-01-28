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
        end
    end
end
