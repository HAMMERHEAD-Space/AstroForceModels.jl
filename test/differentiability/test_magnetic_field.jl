@testset "Magnetic Field Differentiability State" begin
    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
        end
        testname = "Magnetic Field State Dipole " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _mag_dipole_model), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _mag_dipole_model)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Magnetic Field Differentiability Time" begin
    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
        end
        testname = "Magnetic Field Time Dipole " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _mag_dipole_model), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _mag_dipole_model)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Magnetic Field Differentiability Parameters" begin
    for backend in _BACKENDS
        testname = "Magnetic Field Parameters Dipole " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(
                    _state,
                    _p,
                    _t,
                    MagneticFieldAstroModel(;
                        spacecraft_charge_model=FixedChargeMassRatio(x),
                        geomagnetic_field_model=DipoleMagneticField(),
                        eop_data=_eop_data,
                    ),
                ),
                AutoFiniteDiff(),
                _q_over_m,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(
                    acceleration(
                        _state,
                        _p,
                        _t,
                        MagneticFieldAstroModel(;
                            spacecraft_charge_model=FixedChargeMassRatio(x),
                            geomagnetic_field_model=DipoleMagneticField(),
                            eop_data=_eop_data,
                        ),
                    ),
                ),
                backend[2],
                _q_over_m,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end
