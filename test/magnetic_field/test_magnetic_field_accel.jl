@testset "Magnetic Field Acceleration" begin
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

    @testset "IGRF Model basic computation" begin
        charge_model = FixedChargeMassRatio(1e-3)

        mag_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=charge_model,
            geomagnetic_field_model=IGRFField(),
            eop_data=eop_data,
        )

        accel = acceleration(state, p, 0.0, mag_model)
        @test accel isa SVector{3}
        @test all(isfinite, accel)
        @test norm(accel) > 0
    end

    @testset "Dipole Model basic computation" begin
        charge_model = FixedChargeMassRatio(1e-3)

        mag_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=charge_model,
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel = acceleration(state, p, 0.0, mag_model)
        @test accel isa SVector{3}
        @test all(isfinite, accel)
        @test norm(accel) > 0
    end

    @testset "IGRF and Dipole agree in order of magnitude" begin
        charge_model = FixedChargeMassRatio(1e-3)

        igrf_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=charge_model,
            geomagnetic_field_model=IGRFField(),
            eop_data=eop_data,
        )

        dipole_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=charge_model,
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel_igrf = acceleration(state, p, 0.0, igrf_model)
        accel_dipole = acceleration(state, p, 0.0, dipole_model)

        ratio = norm(accel_igrf) / norm(accel_dipole)
        @test 0.1 < ratio < 10.0
    end

    @testset "Acceleration perpendicular to B and v_rel" begin
        charge_model = FixedChargeMassRatio(1e-3)

        mag_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=charge_model,
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel = acceleration(state, p, 0.0, mag_model)

        ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)
        r = SVector{3}(state[1], state[2], state[3])
        v = SVector{3}(state[4], state[5], state[6])
        v_rel = v - cross(ω_vec, r)

        dot_v = dot(accel, v_rel) / (norm(accel) * norm(v_rel))
        @test abs(dot_v) < 1e-10
    end

    @testset "Linearity with charge-to-mass ratio" begin
        q_m_1 = 1e-3
        q_m_2 = 2e-3

        model_1 = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(q_m_1),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        model_2 = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(q_m_2),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel_1 = acceleration(state, p, 0.0, model_1)
        accel_2 = acceleration(state, p, 0.0, model_2)

        @test accel_2 ≈ 2.0 * accel_1 rtol = 1e-14
    end

    @testset "Zero charge gives zero acceleration" begin
        model = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(0.0),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel = acceleration(state, p, 0.0, model)
        @test norm(accel) ≈ 0.0 atol = 1e-30
    end

    @testset "Sign reversal with opposite charge" begin
        model_pos = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(1e-3),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        model_neg = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(-1e-3),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel_pos = acceleration(state, p, 0.0, model_pos)
        accel_neg = acceleration(state, p, 0.0, model_neg)

        @test accel_neg ≈ -accel_pos rtol = 1e-14
    end

    @testset "StateChargeModel" begin
        const_q_m = 1e-3
        state_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=StateChargeModel((u, p, t) -> const_q_m),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        fixed_model = MagneticFieldAstroModel(;
            spacecraft_charge_model=FixedChargeMassRatio(const_q_m),
            geomagnetic_field_model=DipoleMagneticField(),
            eop_data=eop_data,
        )

        accel_state = acceleration(state, p, 0.0, state_model)
        accel_fixed = acceleration(state, p, 0.0, fixed_model)

        @test accel_state ≈ accel_fixed rtol = 1e-14
    end
end
