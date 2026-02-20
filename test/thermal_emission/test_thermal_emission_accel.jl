@testset "Thermal Emission Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    @testset "Consistency with SRP structure" begin
        C_thm = 0.2
        thermal_sat = FixedThermalEmission(C_thm)
        srp_sat = CannonballFixedSRP(C_thm)

        thermal_model = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            eop_data=eop_data,
            shadow_model=NoShadow(),
        )

        srp_model = SRPAstroModel(;
            satellite_srp_model=srp_sat,
            sun_data=sun_model,
            eop_data=eop_data,
            shadow_model=NoShadow(),
        )

        thermal_accel = acceleration(state, p, 0.0, thermal_model)
        srp_accel_val = acceleration(state, p, 0.0, srp_model)

        # With the same coefficient and no shadow, the formulas are identical
        @test thermal_accel ≈ srp_accel_val rtol = 1e-14
    end

    @testset "Direction is along Sun-spacecraft line" begin
        thermal_sat = FixedThermalEmission(0.05)

        thermal_model = ThermalEmissionAstroModel(;
            satellite_thermal_model=thermal_sat,
            sun_data=sun_model,
            eop_data=eop_data,
            shadow_model=NoShadow(),
        )

        accel = acceleration(state, p, 0.0, thermal_model)

        sun_pos = sun_model(AstroForceModels.current_jd(p, 0.0), Position()) ./ 1E3
        r_sc_sun = SVector{3}(
            state[1] - sun_pos[1], state[2] - sun_pos[2], state[3] - sun_pos[3]
        )
        r_sc_sun_hat = r_sc_sun / norm(r_sc_sun)

        accel_hat = accel / norm(accel)
        @test accel_hat ≈ r_sc_sun_hat rtol = 1e-12
    end
end
