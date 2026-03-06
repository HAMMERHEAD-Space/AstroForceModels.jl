@testset "Plasma Drag Acceleration" begin
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
    ] # km, km/s

    @testset "Full model with ChapmanIonosphere" begin
        sat_model = CannonballFixedPlasmaDrag(0.025)
        iono = ChapmanIonosphere()

        model = PlasmaDragAstroModel(;
            satellite_plasma_drag_model=sat_model, ionosphere_model=iono, eop_data=eop_data
        )

        accel = acceleration(state, p, 0.0, model)

        # Plasma drag should be non-zero at LEO
        @test norm(accel) > 0.0

        # Plasma drag is typically 1-2 orders of magnitude smaller than atmospheric drag
        # Atmospheric drag ~1e-8 to 1e-10 km/s², so plasma drag ~1e-10 to 1e-12 km/s²
        mag = norm(accel)
        @test mag > 1e-15
        @test mag < 1e-7

        # Plasma drag should oppose the velocity (negative work)
        velocity = SVector{3}(state[4], state[5], state[6])
        @test dot(accel, velocity) < 0.0

        # Direction should be roughly anti-velocity (like atmospheric drag)
        cos_angle = dot(normalize(accel), normalize(velocity))
        @test cos_angle < -0.9
    end

    @testset "Full model with ConstantIonosphere" begin
        sat_model = CannonballFixedPlasmaDrag(0.025)
        iono = ConstantIonosphere(; rho_i=1e-17)

        model = PlasmaDragAstroModel(;
            satellite_plasma_drag_model=sat_model, ionosphere_model=iono, eop_data=eop_data
        )

        accel = acceleration(state, p, 0.0, model)
        @test norm(accel) > 0.0

        # Compare two different densities: higher density => stronger drag
        iono_high = ConstantIonosphere(; rho_i=1e-15)
        model_high = PlasmaDragAstroModel(;
            satellite_plasma_drag_model=sat_model,
            ionosphere_model=iono_high,
            eop_data=eop_data,
        )
        accel_high = acceleration(state, p, 0.0, model_high)
        @test norm(accel_high) > norm(accel)

        # Ratio should match density ratio (linear in density)
        @test norm(accel_high) / norm(accel) ≈ 1e-15 / 1e-17 rtol = 1e-10
    end

    @testset "NoIonosphere gives zero" begin
        sat_model = CannonballFixedPlasmaDrag(0.025)
        iono = NoIonosphere()

        model = PlasmaDragAstroModel(;
            satellite_plasma_drag_model=sat_model, ionosphere_model=iono, eop_data=eop_data
        )

        accel = acceleration(state, p, 0.0, model)
        @test all(accel .== 0.0)
    end

    @testset "Low-level plasma_drag_accel" begin
        ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)
        rho_i = 1e-17  # kg/m³
        BC_i = 0.025   # m²/kg

        accel = plasma_drag_accel(state, rho_i, BC_i, ω_vec)

        # Same physics as drag_accel — verify consistency
        expected = drag_accel(state, rho_i, BC_i, ω_vec)
        @test accel ≈ expected
    end

    @testset "Comparison with atmospheric drag magnitude" begin
        SpaceIndices.init()

        # Atmospheric drag
        sat_atmo = CannonballFixedDrag(0.2)
        drag_model = DragAstroModel(;
            satellite_drag_model=sat_atmo, atmosphere_model=JB2008(), eop_data=eop_data
        )
        atmo_accel = acceleration(state, p, 0.0, drag_model)

        # Plasma drag
        sat_plasma = CannonballFixedPlasmaDrag(0.025)
        iono = ChapmanIonosphere()
        plasma_model = PlasmaDragAstroModel(;
            satellite_plasma_drag_model=sat_plasma, ionosphere_model=iono, eop_data=eop_data
        )
        plasma_accel = acceleration(state, p, 0.0, plasma_model)

        # Plasma drag should be smaller than atmospheric drag at ~400 km
        @test norm(plasma_accel) < norm(atmo_accel)

        # Both should oppose velocity
        velocity = SVector{3}(state[4], state[5], state[6])
        @test dot(atmo_accel, velocity) < 0.0
        @test dot(plasma_accel, velocity) < 0.0
    end
end
