
@testset "SRP Shadow Models" verbose = true begin
    eop_data = fetch_iers_eop()
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    R_MOD2J2000 = r_eci_to_eci(MOD(), J2000(), JD, eop_data)

    sun_pos = R_MOD2J2000 * sun_position_mod(JD)

    sat_pos = 6800 * normalize(sun_pos)

    shadow_model_set = [Cylindrical(), Conical(), NoShadow(), SmoothedConical()]

    # No Occulsion
    for model in shadow_model_set
        shadow_scaling = shadow_model(
            sat_pos, sun_pos, model; R_Sun=R_SUN, R_Occulting=R_EARTH
        )
        @test shadow_scaling ≈ 1.0 atol = 1E-10
    end

    # Full Occulsion
    for model in [Cylindrical(), Conical()]
        shadow_scaling = shadow_model(
            -sat_pos, sun_pos, model; R_Sun=R_SUN, R_Occulting=R_EARTH
        )
        @test shadow_scaling == 0.0
    end
    @test shadow_model(-sat_pos, sun_pos, NoShadow(); R_Sun=R_SUN, R_Occulting=R_EARTH) ==
        1.0
    # SmoothedConical returns ≈ 0 in full shadow (sigmoid, not exactly 0)
    @test shadow_model(
        -sat_pos, sun_pos, SmoothedConical(); R_Sun=R_SUN, R_Occulting=R_EARTH
    ) ≈ 0.0 atol = 1E-6

    # Partial Occulsion
    sun_pos_simple = [norm(sun_pos), 0.0, 0.0]
    sat_pos_simple = [-20000.0, 6378.1, 0.0]
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, Cylindrical(); R_Sun=R_SUN, R_Occulting=R_EARTH
    ) ≈ 0.0 atol = 1E-6
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, Conical(); R_Sun=R_SUN, R_Occulting=R_EARTH
    ) ≈ 0.2573669402416321 atol = 1E-6
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, NoShadow(); R_Sun=R_SUN, R_Occulting=R_EARTH
    ) ≈ 1.0 atol = 1E-6
    # SmoothedConical should give a partial shadow factor between 0 and 1
    smoothed_partial = shadow_model(
        sat_pos_simple, sun_pos_simple, SmoothedConical(); R_Sun=R_SUN, R_Occulting=R_EARTH
    )
    @test 0.0 < smoothed_partial < 1.0

    @testset "body_pos kwarg (multi-body shadow geometry)" begin
        # A non-origin body at `offset` shifts the entire shadow geometry by `offset`:
        # a spacecraft at position `sat_pos + offset` in the translated geometry should
        # see the same shadow factor as `sat_pos` in the original (central-body) geometry.
        #
        # Note: we intentionally use sun distance in kilometers (1 AU in km) here so the
        # offset-subtraction stays unit-consistent. The outer testset's `sun_pos` is in
        # meters (from `sun_position_mod`), which would trigger floating-point
        # cancellation when subtracting the km-scale `offset` from a ~1.5e11 value.
        sun_pos_t = [1.47e8, 0.0, 0.0]          # km
        sat_pos_t = [-20000.0, 6378.1, 0.0]     # km, partial-eclipse geometry
        offset = SVector{3}(1.0e5, -2.5e5, 5.0e4)

        for model in (Cylindrical(), Conical(), SmoothedConical())
            F_origin = shadow_model(
                sat_pos_t, sun_pos_t, model; R_Sun=R_SUN, R_Occulting=R_EARTH
            )
            F_shifted = shadow_model(
                sat_pos_t + offset,
                sun_pos_t + offset,
                model;
                R_Sun=R_SUN,
                R_Occulting=R_EARTH,
                body_pos=offset,
            )
            @test F_shifted ≈ F_origin atol = 1E-10
        end

        # NoShadow is trivially invariant.
        @test shadow_model(
            sat_pos_t,
            sun_pos_t,
            NoShadow();
            R_Sun=R_SUN,
            R_Occulting=R_EARTH,
            body_pos=offset,
        ) == 1.0
    end

    @testset "Second body casts shadow (Moon-eclipses-Sun geometry)" begin
        # Place a fictitious occulter directly between the spacecraft and the Sun.
        # Spacecraft sees full Sun from the central body, but the second body fully blocks it.
        sun_pos_t = [1.496e8, 0.0, 0.0]
        sat_pos_t = [7000.0, 0.0, 0.0]                # safely outside Earth's shadow
        occ_pos = SVector{3}(1e6, 0.0, 0.0)         # between sat and Sun
        occ_R = 5e5                                # large enough for full eclipse at this geometry

        # Central body (Earth) alone: full Sun.
        @test shadow_model(
            sat_pos_t, sun_pos_t, Conical(); R_Sun=R_SUN, R_Occulting=R_EARTH
        ) ≈ 1.0

        # Second body alone (via body_pos) fully eclipses.
        F_occ = shadow_model(
            sat_pos_t,
            sun_pos_t,
            Conical();
            R_Sun=R_SUN,
            R_Occulting=occ_R,
            body_pos=occ_pos,
        )
        @test F_occ ≈ 0.0 atol = 1E-10
    end

    @testset "SmoothedConical Properties" begin
        # Default constructor
        sc = SmoothedConical()
        @test sc.cs == 289.78
        @test sc.ct == 1.0

        # Custom parameters
        sc_custom = SmoothedConical(500.0, 0.8)
        @test sc_custom.cs == 500.0
        @test sc_custom.ct == 0.8

        # Monotonicity: shadow factor should increase as satellite moves from shadow to sunlight
        # Move satellite perpendicular to Sun direction at varying offsets
        sun_dir = normalize(sun_pos)
        perp = normalize(cross(sun_dir, [0.0, 0.0, 1.0]))
        factors = Float64[]
        for offset in range(-2.0, 2.0; length=20)
            test_pos = -6800.0 * sun_dir + offset * R_EARTH * perp
            f = shadow_model(
                test_pos, sun_pos, SmoothedConical(); R_Sun=R_SUN, R_Occulting=R_EARTH
            )
            push!(factors, f)
        end
        # Factors should be non-decreasing then non-increasing (peak in the middle when outside shadow)
        # At minimum, the range should include values near 0 and near 1
        @test minimum(factors) < 0.1
        @test maximum(factors) > 0.9
    end
end
