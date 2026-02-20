@testset "Solid Body Tides Differentiability State" begin
    for backend in _BACKENDS
        testname = "Solid Body Tides State " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _tides_model), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _tides_model)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 2e-1
        end
    end
end

@testset "Solid Body Tides Differentiability Time" begin
    for backend in _BACKENDS
        testname = "Solid Body Tides Time " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _tides_model), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _tides_model)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Solid Body Tides Differentiability Parameters" begin
    for backend in _BACKENDS
        testname = "Solid Body Tides k2 " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(
                    _state,
                    _p,
                    _t,
                    SolidBodyTidesModel(;
                        tide_raising_bodies=_tides_model.tide_raising_bodies, k2=x
                    ),
                ),
                AutoFiniteDiff(),
                _tides_model.k2,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(
                    acceleration(
                        _state,
                        _p,
                        _t,
                        SolidBodyTidesModel(;
                            tide_raising_bodies=_tides_model.tide_raising_bodies, k2=x
                        ),
                    ),
                ),
                backend[2],
                _tides_model.k2,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end

        testname = "Solid Body Tides k3 " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(
                    _state,
                    _p,
                    _t,
                    SolidBodyTidesModel(;
                        tide_raising_bodies=_tides_model.tide_raising_bodies, k3=x
                    ),
                ),
                AutoFiniteDiff(),
                _tides_model.k3,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(
                    acceleration(
                        _state,
                        _p,
                        _t,
                        SolidBodyTidesModel(;
                            tide_raising_bodies=_tides_model.tide_raising_bodies, k3=x
                        ),
                    ),
                ),
                backend[2],
                _tides_model.k3,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end
