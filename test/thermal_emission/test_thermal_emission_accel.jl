@testset "Thermal Emission Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = create_test_params(; JD=JD, eop_data=eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_model = test_sun_model()

    @testset "Consistency with SRP structure" begin
        C_thm = 0.2
        thermal_sat = FixedThermalEmission(C_thm)
        srp_sat = CannonballFixedSRP(C_thm)

        thermal_model = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            shadow_model=NoShadow(),
            R_Occulting=AstroForceModels.R_EARTH,
        )

        srp_model = SRPAstroModel(;
            satellite_srp_model=srp_sat,
            sun_data=sun_model,
            shadow_model=NoShadow(),
            R_Occulting=AstroForceModels.R_EARTH,
        )

        thermal_accel = acceleration(state, p, 0.0, thermal_model)
        srp_accel_val = acceleration(state, p, 0.0, srp_model)

        # With the same coefficient and no shadow, the formulas are identical
        @test thermal_accel ≈ srp_accel_val rtol = 1e-14
    end

    @testset "additional_occulters: empty tuple equals single-body result" begin
        thermal_sat = FixedThermalEmission(0.03)
        m1 = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            shadow_model=NoShadow(),
            R_Occulting=AstroForceModels.R_EARTH,
        )
        m2 = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            shadow_model=NoShadow(),
            R_Occulting=AstroForceModels.R_EARTH,
            additional_occulters=(),
        )
        @test acceleration(state, p, 0.0, m1) ≈ acceleration(state, p, 0.0, m2) rtol = 1e-14
    end

    @testset "additional_occulters: body larger than Sun's disk zeros the acceleration" begin
        # A body at the Sun's position with radius ≫ R_SUN but ≪ spacecraft-Sun
        # distance must fully eclipse the Sun via the Conical model. See the
        # corresponding SRP test for the geometric rationale.
        thermal_sat = FixedThermalEmission(0.03)
        sun_tb = test_sun_model(; frames=p.frames)
        giant = OccultingBody(sun_tb, 1.0e7)

        model = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            shadow_model=Conical(),
            R_Occulting=AstroForceModels.R_EARTH,
            additional_occulters=(giant,),
        )
        a = acceleration(state, p, 0.0, model)
        @test norm(a) ≈ 0.0 atol = 1e-15
    end

    @testset "Direction is along Sun-spacecraft line" begin
        thermal_sat = FixedThermalEmission(0.05)

        thermal_model = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            shadow_model=NoShadow(),
            R_Occulting=AstroForceModels.R_EARTH,
        )

        accel = acceleration(state, p, 0.0, thermal_model)

        t_ft = AstroForceModels.ft_time(p, 0.0)
        sun_pos = get_position(sun_model.ephem_type, sun_model.body, p.frames, t_ft)
        r_sc_sun = SVector{3}(
            state[1] - sun_pos[1], state[2] - sun_pos[2], state[3] - sun_pos[3]
        )
        r_sc_sun_hat = r_sc_sun / norm(r_sc_sun)

        accel_hat = accel / norm(accel)
        @test accel_hat ≈ r_sc_sun_hat rtol = 1e-12
    end
end
