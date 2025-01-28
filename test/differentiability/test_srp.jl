@testset "SRP Differentiability State" begin
    for backend in _BACKENDS
        testname = "SRP Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _srp_model), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _srp_model)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "SRP Differentiability Time" begin
    for backend in _BACKENDS
        testname = "SRP Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _srp_model), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _srp_model)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "SRP Differentiability SRP Parameters" begin
    for backend in _BACKENDS
        testname = "SRP Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(
                    _state,
                    _p,
                    _t,
                    SRPAstroModel(;
                        satellite_srp_model=CannonballFixedSRP(x),
                        sun_data=_sun_model,
                        eop_data=_eop_data,
                        shadow_model=Conical(),
                    ),
                ),
                AutoFiniteDiff(),
                _RC,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(
                    acceleration(
                        _state,
                        _p,
                        _t,
                        SRPAstroModel(;
                            satellite_srp_model=CannonballFixedSRP(x),
                            sun_data=_sun_model,
                            eop_data=_eop_data,
                            shadow_model=Conical(),
                        ),
                    ),
                ),
                backend[2],
                _RC,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end
