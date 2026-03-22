@testset "Albedo Differentiability State" begin
    for backend in _BACKENDS
        testname = "Albedo Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, _p, _t, _albedo_model), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, _p, _t, _albedo_model)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 2e-1
        end
    end
end

@testset "Albedo Differentiability Time" begin
    for backend in _BACKENDS
        testname = "Albedo Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, x, _albedo_model), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, x, _albedo_model)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end

@testset "Albedo Differentiability Albedo Parameters" begin
    # Pre-compute the Lebedev quadrature and surface positions once — they are constants
    # and should not be differentiated through (Lebedev.jl uses in-place mutation).
    _albedo_ref = AlbedoAstroModel(;
        satellite_shape_model=CannonballFixedSRP(_RC),
        sun_data=_sun_model,
        body_albedo_model=_uniform_albedo_model,
        body_fixed_frame=:ITRF,
        propagation_frame=:ICRF,
    )

    for backend in _BACKENDS
        testname = "Albedo Differentiability " * backend[1]
        @testset "$testname" begin
            # Rebuild only the satellite_shape_model inside the closure, reusing
            # pre-computed surface positions and weights from _albedo_ref.
            function _albedo_with_rc(x)
                AlbedoAstroModel(
                    CannonballFixedSRP(x),
                    _albedo_ref.sun_data,
                    _albedo_ref.body_albedo_model,
                    _albedo_ref.body_fixed_frame,
                    _albedo_ref.propagation_frame,
                    _albedo_ref.solar_irradiance,
                    _albedo_ref.speed_of_light,
                    _albedo_ref.AU,
                    _albedo_ref.surface_positions_ecef,
                    _albedo_ref.weights,
                    _albedo_ref.compiled_rotation3,
                )
            end

            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(_state, _p, _t, _albedo_with_rc(x)),
                AutoFiniteDiff(),
                _RC,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(_state, _p, _t, _albedo_with_rc(x))),
                backend[2],
                _RC,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end
