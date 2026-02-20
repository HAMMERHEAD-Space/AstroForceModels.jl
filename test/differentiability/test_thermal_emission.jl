@testset "Thermal Emission Differentiability State" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            if backend[1] == "Enzyme"
                backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
            end
            testname = "Thermal Emission State " * shadow[1] * " " * backend[1]
            @testset "$testname" begin
                thermal_model = ThermalEmissionAstroModel(;
                    satellite_thermal_model=_thermal_sat_model,
                    sun_data=_sun_model,
                    eop_data=_eop_data,
                    shadow_model=shadow[2],
                )

                f_fd, df_fd = value_and_jacobian(
                    (x) -> acceleration(x, _p, _t, thermal_model), AutoFiniteDiff(), _state
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(acceleration(x, _p, _t, thermal_model)), backend[2], _state
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-10
            end
        end
    end
end

@testset "Thermal Emission Differentiability Time" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            if backend[1] == "Enzyme"
                backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
            end
            testname = "Thermal Emission Time " * shadow[1] * " " * backend[1]
            @testset "$testname" begin
                thermal_model = ThermalEmissionAstroModel(;
                    satellite_thermal_model=_thermal_sat_model,
                    sun_data=_sun_model,
                    eop_data=_eop_data,
                    shadow_model=shadow[2],
                )

                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(_state, _p, x, thermal_model), AutoFiniteDiff(), _t
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(acceleration(_state, _p, x, thermal_model)), backend[2], _t
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-10
            end
        end
    end
end

@testset "Thermal Emission Differentiability Parameters" begin
    for shadow in _SHADOW_MODELS
        for backend in _BACKENDS
            testname = "Thermal Emission Parameters " * shadow[1] * " " * backend[1]
            @testset "$testname" begin
                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(
                        _state,
                        _p,
                        _t,
                        ThermalEmissionAstroModel(;
                            satellite_thermal_model=FixedThermalEmission(x),
                            sun_data=_sun_model,
                            eop_data=_eop_data,
                            shadow_model=shadow[2],
                        ),
                    ),
                    AutoFiniteDiff(),
                    _C_thm,
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(
                        acceleration(
                            _state,
                            _p,
                            _t,
                            ThermalEmissionAstroModel(;
                                satellite_thermal_model=FixedThermalEmission(x),
                                sun_data=_sun_model,
                                eop_data=_eop_data,
                                shadow_model=shadow[2],
                            ),
                        ),
                    ),
                    backend[2],
                    _C_thm,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-10
            end
        end
    end
end
