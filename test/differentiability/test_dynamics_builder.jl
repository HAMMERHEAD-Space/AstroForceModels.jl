@testset "Dynamics Builder State Differentiability" begin
    SpaceIndices.init()

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "Dynamics Builder State Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> build_dynamics_model(x, _p, _t, _model_list),
                AutoFiniteDiff(),
                _state,
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(build_dynamics_model(x, _p, _t, _model_list)),
                backend[2],
                _state,
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
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "Dynamics Builder Time Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> build_dynamics_model(_state, _p, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(build_dynamics_model(_state, _p, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-1
        end
    end
    SpaceIndices.destroy()
end

@testset "Dynamics Builder Parameter Differentiability" begin
    SpaceIndices.init()

    function dynamics_params(x::AbstractArray{T}) where {T<:Number}
        satellite_srp_model = CannonballFixedSRP(x[1])
        srp_model = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=_sun_model,
            eop_data=_eop_data,
            shadow_model=Conical(),
        )

        satellite_drag_model = CannonballFixedDrag(x[2])
        drag_model = DragAstroModel(;
            satellite_drag_model=satellite_drag_model,
            atmosphere_model=JB2008(),
            eop_data=_eop_data,
        )

        models = (_sun_model, _moon_model, srp_model, drag_model)
        model_list = CentralBodyDynamicsModel(_grav_model, models)

        return Array(build_dynamics_model(_state, _p, _t, model_list))
    end

    p_spacecraft = [0.2, 0.2]

    for backend in _BACKENDS
        testname = "Dynamics Builder Parameter Differentiability " * backend[1]
        @testset "$testname" begin
            if backend[1] == "Enzyme"
                backend = (
                    "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
                )
            end
            f_fd, df_fd = value_and_jacobian(
                dynamics_params, AutoFiniteDiff(), p_spacecraft
            )

            f_ad, df_ad = value_and_jacobian(dynamics_params, backend[2], p_spacecraft)

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-9
        end
    end
    SpaceIndices.destroy()
end
