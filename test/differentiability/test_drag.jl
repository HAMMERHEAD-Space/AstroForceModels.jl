@testset "Drag Differentiability State" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Enzyme" && atmo[1] ∈ _ENZYME_RUNTIME_ACTIVITY)
                    backend = (
                        "Enzyme",
                        AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)),
                    )
                end

                f_fd, df_fd = value_and_jacobian(
                    (x) -> acceleration(
                        x,
                        _p,
                        _t,
                        DragAstroModel(;
                            satellite_drag_model=_satellite_drag_model,
                            atmosphere_model=atmo[2],
                            eop_data=_eop_data,
                        ),
                    ),
                    AutoFiniteDiff(),
                    _state,
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(
                        acceleration(
                            x,
                            _p,
                            _t,
                            DragAstroModel(;
                                satellite_drag_model=_satellite_drag_model,
                                atmosphere_model=atmo[2],
                                eop_data=_eop_data,
                            ),
                        ),
                    ),
                    backend[2],
                    _state,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-5
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Drag Differentiability Time" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Enzyme" && atmo[1] ∈ _ENZYME_RUNTIME_ACTIVITY)
                    backend = (
                        "Enzyme",
                        AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)),
                    )
                end
                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(
                        _state,
                        _p,
                        x,
                        DragAstroModel(;
                            satellite_drag_model=_satellite_drag_model,
                            atmosphere_model=atmo[2],
                            eop_data=_eop_data,
                        ),
                    ),
                    AutoFiniteDiff(),
                    _t,
                )

                if !(backend[1] == "Zygote" && atmo[1] == "NoAtmosphere")
                    f_ad, df_ad = value_and_derivative(
                        (x) -> Array(
                            acceleration(
                                _state,
                                _p,
                                x,
                                DragAstroModel(;
                                    satellite_drag_model=_satellite_drag_model,
                                    atmosphere_model=atmo[2],
                                    eop_data=_eop_data,
                                ),
                            ),
                        ),
                        backend[2],
                        _t,
                    )
                    @test f_fd ≈ f_ad
                    @test df_fd ≈ df_ad atol = 2e-1
                else
                    try
                        f_ad, df_ad = value_and_derivative(
                            (x) -> Array(
                                acceleration(
                                    _state,
                                    _p,
                                    x,
                                    DragAstroModel(;
                                        satellite_drag_model=_satellite_drag_model,
                                        atmosphere_model=atmo[2],
                                        eop_data=_eop_data,
                                    ),
                                ),
                            ),
                            backend[2],
                            _t,
                        )
                    catch err
                        @test err isa Exception
                        @test startswith(
                            sprint(showerror, err),
                            "MethodError: no method matching iterate(::Nothing)",
                        )
                    end
                end
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Drag Differentiability Model Parameters" begin
    SpaceIndices.init()

    for backend in _BACKENDS
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Enzyme")
                    backend = (
                        "Enzyme",
                        AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)),
                    )
                end

                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(
                        _state,
                        _p,
                        _t,
                        DragAstroModel(;
                            satellite_drag_model=CannonballFixedDrag(x),
                            atmosphere_model=atmo[2],
                            eop_data=_eop_data,
                        ),
                    ),
                    AutoFiniteDiff(),
                    _BC,
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(
                        acceleration(
                            _state,
                            _p,
                            _t,
                            DragAstroModel(;
                                satellite_drag_model=CannonballFixedDrag(x),
                                atmosphere_model=atmo[2],
                                eop_data=_eop_data,
                            ),
                        ),
                    ),
                    backend[2],
                    _BC,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad rtol = 1e-4
            end
        end
    end
    SpaceIndices.destroy()
end
