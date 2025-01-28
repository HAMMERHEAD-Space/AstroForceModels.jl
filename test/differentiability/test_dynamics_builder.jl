@testset "Dynamics Builder State Differentiability" begin
    SpaceIndices.init()

    for backend in _BACKENDS
        testname = "Dynamics Builder State Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> build_dynamics_model(x, _p, _t, _model_list), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(build_dynamics_model(x, _p, _t, _model_list)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-4
        end
    end
    SpaceIndices.destroy()
end

@testset "Dynamics Builder Time Differentiability" begin
    SpaceIndices.init()

    for backend in _BACKENDS
        testname = "Dynamics Builder Time Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> build_dynamics_model(_state, _p, x, _model_list), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(build_dynamics_model(_state, _p, x, _model_list)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-1
        end
    end
    SpaceIndices.destroy()
end

@testset "Dynamics Builder Parameter Differentiability" begin
    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    JD = date_to_jd(2024, 1, 5, 12, 32, 0.0)
    p = ComponentVector(; JD=JD)

    function dynamics_params(x::AbstractArray)
        grav_model = GravityHarmonicsAstroModel(;
            gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
        )
        sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
        moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

        satellite_srp_model = CannonballFixedSRP(x[1])
        srp_model = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=sun_third_body,
            eop_data=eop_data,
            shadow_model=Conical(),
        )

        satellite_drag_model = CannonballFixedDrag(x[2])
        drag_model = DragAstroModel(;satellite_drag_model=satellite_drag_model, atmosphere_model=JB2008(), eop_data=eop_data)

        state = [
            -1076.225324679696
            -6765.896364327722
            -332.3087833503755
            9.356857417032581
            -3.3123476319597557
            -1.1880157328553503
        ] #km, km/s

        t = 0.0
        model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

        return Array(build_dynamics_model(state, p, t, model_list))
    end

    p_spacecraft = [0.2, 0.2]

    for backend in _BACKENDS
        testname = "Dynamics Builder Parameter Differentiability " * backend[1]
        @testset "$testname" begin
            if backend[1] == "Enzyme"
                backend = ("Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)))
            end
            f_fd, df_fd = value_and_jacobian(
                dynamics_params, AutoFiniteDiff(), p_spacecraft
            )

            f_ad, df_ad = value_and_jacobian(dynamics_params, backend[2], p_spacecraft)

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 1e-2
        end
    end
    SpaceIndices.destroy()
end
