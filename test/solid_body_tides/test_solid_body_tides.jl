@testset "Solid Body Tides Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    @testset "Regression" begin
        model = SolidBodyTidesModel(eop_data)
        accel = acceleration(state, p, t, model)

        # Cross-validated against Orekit 13.1.4 SolidTides (IERS 2010, no pole tide)
        # at the same epoch and Cartesian state. Orekit computes acceleration via
        # spherical harmonic coefficient corrections (ΔC_nm, ΔS_nm) evaluated with
        # Holmes-Featherstone — a completely independent computational path from
        # our direct Legendre addition-theorem gradient.
        #
        # Orekit reference (km/s²):
        #   ax = -8.020132521332073e-11
        #   ay =  1.280231565156116e-10
        #   az = -6.803247555720357e-11
        #
        # Our Step 1 (frequency-independent) implementation agrees in direction
        # and magnitude. The ~5% norm difference is accounted for by:
        #   - Orekit includes Step 2 frequency-dependent corrections (ΔC₂₀,
        #     ΔC₂₁, ΔS₂₁, ΔC₂₂, ΔS₂₂)
        #   - Orekit removes the permanent tide for zero-tide gravity fields
        #   - Orekit uses order-dependent complex Love numbers (k₂₀≠k₂₁≠k₂₂)
        #   - Different ephemeris sources (JPL DE vs SatelliteToolboxCelestialBodies,
        #     confirmed to contribute < 0.01% of the difference)
        orekit_ref = [-8.020132521332073e-11, 1.280231565156116e-10, -6.803247555720357e-11]
        orekit_dir = orekit_ref / norm(orekit_ref)
        accel_dir = accel / norm(accel)
        @test dot(accel_dir, orekit_dir) ≈ 1.0 atol = 0.005
        @test norm(accel) / norm(orekit_ref) ≈ 1.0 atol = 0.10
    end

    @testset "Physical consistency" begin
        jd = AstroForceModels.current_jd(p, t)
        r_sun = sun_model(jd, AstroForceModels.Position()) ./ 1E3
        r_moon = moon_model(jd, AstroForceModels.Position()) ./ 1E3

        # Isolate Moon contribution by passing only the Moon
        accel_moon = solid_body_tides_accel(
            state, ((r_moon, AstroForceModels.μ_MOON),); include_degree_3=false
        )
        # Isolate Sun contribution by passing only the Sun
        accel_sun = solid_body_tides_accel(
            state, ((r_sun, AstroForceModels.μ_SUN),); include_degree_3=false
        )

        ratio = norm(accel_moon) / norm(accel_sun)
        @test ratio > 0.5
        @test ratio < 5.0

        # Acceleration scales linearly with k2
        bodies = ((r_sun, AstroForceModels.μ_SUN), (r_moon, AstroForceModels.μ_MOON))
        accel_k2_1 = solid_body_tides_accel(state, bodies; k2=0.302, include_degree_3=false)
        accel_k2_2 = solid_body_tides_accel(state, bodies; k2=0.604, include_degree_3=false)
        @test accel_k2_2 ≈ 2.0 * accel_k2_1 rtol = 1e-14

        # Zero Love number should give zero acceleration
        accel_zero = solid_body_tides_accel(state, bodies; k2=0.0, k3=0.0)
        @test norm(accel_zero) ≈ 0.0 atol = 1e-30
    end

    @testset "Degree 2 analytical verification" begin
        Re = AstroForceModels.R_EARTH
        k2 = 0.30190
        μ_body = AstroForceModels.μ_MOON

        r_sat_val = 7000.0
        r_body_val = 384400.0

        u_simple = [r_sat_val, 0.0, 0.0, 0.0, 7.5, 0.0]
        r_body_vec = SVector{3}(r_body_val, 0.0, 0.0)

        # For collinear case (γ = 0, ξ = 1):
        #   a_x = (3 k2 μ Re^5) / (r_body^3 r^4) * (-2 + 1) = -(3 k2 μ Re^5) / (r_body^3 r^4)
        expected_ax = -3.0 * k2 * μ_body * Re^5 / (r_body_val^3 * r_sat_val^4)

        accel = solid_body_tides_accel(
            u_simple, ((r_body_vec, μ_body),); k2=k2, include_degree_3=false
        )

        @test accel[1] ≈ expected_ax rtol = 1e-10
        @test abs(accel[2]) < 1e-30
        @test abs(accel[3]) < 1e-30
    end

    @testset "Degree 3 analytical verification" begin
        # For collinear case (γ = 0, ξ = 1):
        #   a_x = (k3 μ Re^7) / (2 r_body^4 r^5) * (-20 + 12) = -4 k3 μ Re^7 / (r_body^4 r^5)

        Re = AstroForceModels.R_EARTH
        k3 = 0.093
        μ_body = AstroForceModels.μ_MOON

        r_sat_val = 7000.0
        r_body_val = 384400.0

        u_simple = [r_sat_val, 0.0, 0.0, 0.0, 7.5, 0.0]
        r_body_vec = SVector{3}(r_body_val, 0.0, 0.0)

        expected_ax = -4.0 * k3 * μ_body * Re^7 / (r_body_val^4 * r_sat_val^5)

        accel = solid_body_tides_accel(
            u_simple, ((r_body_vec, μ_body),); k2=0.0, k3=k3, include_degree_3=true
        )

        @test accel[1] ≈ expected_ax rtol = 1e-10
        @test abs(accel[2]) < 1e-30
        @test abs(accel[3]) < 1e-30
    end
end
