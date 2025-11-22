@testset "Albedo Acceleration Tests" begin
    @testset "Reference Values from Literature" begin
        # Test based on Vielberg & Kusche (2020) "Extended forward and inverse modeling 
        # of radiation pressure accelerations for LEO satellites"
        # Reference: For GRACE during January 2010, average ERP acceleration ≈ 1.4e-8 m/s²

        # GRACE-like satellite configuration (January 1, 2010, 12:00 UTC)
        JD = date_to_jd(2010, 1, 1, 16, 0, 0.0)
        p = ComponentVector(; JD=JD)

        # Fetch EOP data for coordinate transformations
        eop_data = fetch_iers_eop()

        # GRACE-like orbital state (altitude ~460 km, polar orbit)
        # Position in ECI frame [km]
        state = [6834.0, 0.0, 0.0, 0.0, 7.63715140464081, 0.0]

        # GRACE-like satellite properties
        # GRACE specs: approx 432 kg mass, approx 3.1×2.0×0.7m body + solar panels
        # Effective cross-sectional area: approx 8-12 m² (depending on orientation)
        # Area-to-mass ratio: approx 8-12 m²/432 kg ≈ 0.019-0.028 m²/kg
        # With reflectivity factor approx 1.3-1.8: RC ≈ 0.025-0.050 m²/kg
        # Using conservative middle value for GRACE
        satellite_shape_model = CannonballFixedSRP(0.030)  # m²/kg

        # Sun position data
        sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

        # Uniform albedo model with literature-based values
        # visible_albedo=0.3, infrared_emissivity=0.7 (Earth-averaged values)
        uniform_albedo_model = UniformAlbedoModel(
            visible_albedo=0.3, infrared_emissivity=0.7
        )

        # Create albedo force model
        albedo_model = AlbedoAstroModel(;
            satellite_shape_model=satellite_shape_model,
            sun_data=sun_model,
            body_albedo_model=uniform_albedo_model,
            eop_data=eop_data,
        )

        # Compute albedo acceleration
        albedo_accel = acceleration(state, p, 0.0, albedo_model)

        # Expected magnitude based on Vielberg & Kusche (2020)
        # Earth Radiation Pressure: ~1.4e-8 m/s² for GRACE
        expected_magnitude = 1.4e-11  # km/s²

        # Compute actual magnitude
        actual_magnitude = norm(albedo_accel)

        # Test that the magnitude is within reasonable bounds
        # Allow for factor of 2-3 variation due to:
        # - Different orbital geometry
        # - Simplified integration parameters
        # - Earth orientation differences
        # - Seasonal variations
        @test actual_magnitude > expected_magnitude / 3.0
        @test actual_magnitude < expected_magnitude * 3.0

        # Additional sanity checks
        @test length(albedo_accel) == 3
        @test all(isfinite.(albedo_accel))
        @test actual_magnitude > 0.0
    end
end
