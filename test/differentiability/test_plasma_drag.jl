@testset "Plasma Drag Differentiability State" begin
    for backend in _BACKENDS
        testname = "Plasma Drag State " * backend[1]
        @testset "$testname" begin
            for iono in _IONOSPHERE_MODELS
                iononame = iono[1]
                @testset "$iononame" begin
                    model = PlasmaDragAstroModel(;
                        satellite_plasma_drag_model=_satellite_plasma_drag_model,
                        ionosphere_model=iono[2],
                        eop_data=_eop_data,
                    )

                    f_fd, df_fd = value_and_jacobian(
                        (x) -> acceleration(x, _p, _t, model), AutoFiniteDiff(), _state
                    )

                    f_ad, df_ad = value_and_jacobian(
                        (x) -> Array(acceleration(x, _p, _t, model)), backend[2], _state
                    )

                    @test f_fd ≈ f_ad
                    @test df_fd ≈ df_ad atol = 1e-5
                end
            end
        end
    end
end

@testset "Plasma Drag Differentiability Time" begin
    for backend in _BACKENDS
        testname = "Plasma Drag Time " * backend[1]
        @testset "$testname" begin
            for iono in _IONOSPHERE_MODELS
                iononame = iono[1]
                if backend[1] == "Zygote" && iononame == "NoIonosphere"
                    continue
                end
                @testset "$iononame" begin
                    model = PlasmaDragAstroModel(;
                        satellite_plasma_drag_model=_satellite_plasma_drag_model,
                        ionosphere_model=iono[2],
                        eop_data=_eop_data,
                    )

                    f_fd, df_fd = value_and_derivative(
                        (x) -> acceleration(_state, _p, x, model), AutoFiniteDiff(), _t
                    )

                    f_ad, df_ad = value_and_derivative(
                        (x) -> Array(acceleration(_state, _p, x, model)), backend[2], _t
                    )

                    @test f_fd ≈ f_ad
                    @test df_fd ≈ df_ad atol = 2e-1
                end
            end
        end
    end
end

@testset "Plasma Drag Differentiability Model Parameters" begin
    for backend in _BACKENDS
        testname = "Plasma Drag Parameters " * backend[1]
        @testset "$testname" begin
            for iono in _IONOSPHERE_MODELS
                iononame = iono[1]
                @testset "$iononame" begin
                    f_fd, df_fd = value_and_derivative(
                        (x) -> acceleration(
                            _state,
                            _p,
                            _t,
                            PlasmaDragAstroModel(;
                                satellite_plasma_drag_model=CannonballFixedPlasmaDrag(x),
                                ionosphere_model=iono[2],
                                eop_data=_eop_data,
                            ),
                        ),
                        AutoFiniteDiff(),
                        _BC_i,
                    )

                    f_ad, df_ad = value_and_derivative(
                        (x) -> Array(
                            acceleration(
                                _state,
                                _p,
                                _t,
                                PlasmaDragAstroModel(;
                                    satellite_plasma_drag_model=CannonballFixedPlasmaDrag(
                                        x
                                    ),
                                    ionosphere_model=iono[2],
                                    eop_data=_eop_data,
                                ),
                            ),
                        ),
                        backend[2],
                        _BC_i,
                    )

                    @test f_fd ≈ f_ad
                    @test df_fd ≈ df_ad rtol = 1e-4
                end
            end
        end
    end
end
