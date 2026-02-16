@testset "SRP Differentiability State" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            testname = "SRP State " * shadow[1] * " " * backend[1]
            @testset "$testname" begin
                srp_model = SRPAstroModel(;
                    satellite_srp_model=_satellite_srp_model,
                    sun_data=_sun_model,
                    eop_data=_eop_data,
                    shadow_model=shadow[2],
                )

                f_fd, df_fd = value_and_jacobian(
                    (x) -> acceleration(x, _p, _t, srp_model), AutoFiniteDiff(), _state
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(acceleration(x, _p, _t, srp_model)), backend[2], _state
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-10
            end
        end
    end
end

@testset "SRP Differentiability Time" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            testname = "SRP Time " * shadow[1] * " " * backend[1]
            @testset "$testname" begin
                srp_model = SRPAstroModel(;
                    satellite_srp_model=_satellite_srp_model,
                    sun_data=_sun_model,
                    eop_data=_eop_data,
                    shadow_model=shadow[2],
                )

                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(_state, _p, x, srp_model), AutoFiniteDiff(), _t
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(acceleration(_state, _p, x, srp_model)), backend[2], _t
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-10
            end
        end
    end
end

@testset "SRP Differentiability SRP Parameters" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            testname = "SRP Parameters " * shadow[1] * " " * backend[1]
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
                            shadow_model=shadow[2],
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
                                shadow_model=shadow[2],
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
end
