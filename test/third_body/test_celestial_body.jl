@testset "Celestial Body Construction" begin
    sun = SunBody()
    @test sun.μ == μ_SUN
    @test nameof(sun) == :Sun
    @test sun.Req == R_SUN

    moon = MoonBody()
    @test moon.μ == μ_MOON
    @test nameof(moon) == :Moon
    @test moon.Req == R_MOON
end

@testset "Test Low Precision Celestial Body Construction" begin
    sun = SunBody(; T=Float32)
    @test sun.μ == Float32(μ_SUN)

    typeof(sun.μ)

    @test nameof(sun) == :Sun
    @test sun.Req == Float32(R_SUN)

    moon = MoonBody(; T=Float32)
    @test moon.μ == Float32(μ_MOON)
    @test nameof(moon) == :Moon
    @test moon.Req == Float32(R_MOON)
end

@testset "Test Custom Celestial Body Construction" begin
    tatooine = CelestialBody(:Tatooine, :Dual_Sun, 99999, 1.2E8, 3456.0)

    @test nameof(tatooine) == :Tatooine
    @test tatooine.central_body == :Dual_Sun
    @test tatooine.jpl_code == 99999
    @test tatooine.μ == 1.2E8
    @test tatooine.Req == 3456
end
