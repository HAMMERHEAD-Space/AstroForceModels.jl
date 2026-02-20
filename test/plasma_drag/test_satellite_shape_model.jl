@testset "Plasma Drag Satellite Shape Models" begin
    @testset "CannonballFixedPlasmaDrag" begin
        # From precomputed ballistic coefficient
        model = CannonballFixedPlasmaDrag(0.025)
        @test model.ballistic_coeff ≈ 0.025

        # From physical properties
        r = 1.0   # m
        m = 500.0 # kg
        Cd = 4.0
        model2 = CannonballFixedPlasmaDrag(r, m, Cd)
        @test model2.ballistic_coeff ≈ Cd * π * r^2 / m

        # Ballistic coefficient accessor
        u = zeros(6)
        p = ComponentVector(; JD=0.0)
        @test ion_ballistic_coefficient(u, p, 0.0, model2) ≈ model2.ballistic_coeff

        # Argument validation
        @test_throws ArgumentError CannonballFixedPlasmaDrag(-1.0)
        @test_throws ArgumentError CannonballFixedPlasmaDrag(1.0, 1.0, -1.0)
        @test_throws ArgumentError CannonballFixedPlasmaDrag(-1.0, 1.0, 1.0)
        @test_throws ArgumentError CannonballFixedPlasmaDrag(1.0, -1.0, 1.0)
    end

    @testset "StatePlasmaDragModel" begin
        model = StatePlasmaDragModel()
        u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05]
        p = ComponentVector(; JD=0.0)
        @test ion_ballistic_coefficient(u, p, 0.0, model) ≈ 0.05

        # Custom index
        model8 = StatePlasmaDragModel(8)
        u8 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03]
        @test ion_ballistic_coefficient(u8, p, 0.0, model8) ≈ 0.03
    end
end
