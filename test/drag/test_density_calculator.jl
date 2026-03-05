# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Satellite Toolbox already tests the density calculations. These tests just check if our 
# wrapper function returns the same value, avoids the invalid regions, throws proper errors.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
@testset "Compute Density" verbose = true begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    pos_eci = [-1076.225324679696; -6765.896364327722; -332.3087833503755] #km
    R_J2000_ITRF = r_eci_to_ecef(J2000(), ITRF(), JD, eop_data)
    pos_ecef = R_J2000_ITRF * pos_eci

    ϕ_gd, λ_gd, h = ecef_to_geodetic(pos_ecef * 1E3)

    # JB2008 Check
    rho = AtmosphericModels.jb2008(JD, ϕ_gd, λ_gd, h).total_density
    ρ = compute_density(JD, pos_eci, eop_data, JB2008())
    @test rho == ρ

    # JR1971 Check
    rho = AtmosphericModels.jr1971(JD, ϕ_gd, λ_gd, h).total_density
    ρ = compute_density(JD, pos_eci, eop_data, JR1971())
    @test rho == ρ

    # MSIS2000 Check
    rho = AtmosphericModels.nrlmsise00(JD, h, ϕ_gd, λ_gd).total_density
    ρ = compute_density(JD, pos_eci, eop_data, MSIS2000())
    @test rho == ρ

    # Exponential Check
    rho = AtmosphericModels.exponential(h)
    ρ = compute_density(JD, pos_eci, eop_data, ExpAtmo())
    @test rho == ρ

    # Harris-Priester Check
    rho = AtmosphericModels.harrispriester(JD, ϕ_gd, λ_gd, h)
    ρ = compute_density(JD, pos_eci, eop_data, HarrisPriester())
    @test rho == ρ

    # Harris-Priester with custom n
    rho = AtmosphericModels.harrispriester(JD, ϕ_gd, λ_gd, h; n=2)
    ρ = compute_density(JD, pos_eci, eop_data, HarrisPriester(; n=2))
    @test rho == ρ

    # Modified Harris-Priester Check
    rho = AtmosphericModels.harrispriester_modified(JD, ϕ_gd, λ_gd, h)
    ρ = compute_density(JD, pos_eci, eop_data, HarrisPriesterModified())
    @test rho == ρ

    # Modified Harris-Priester with custom n
    rho = AtmosphericModels.harrispriester_modified(JD, ϕ_gd, λ_gd, h; n=6)
    ρ = compute_density(JD, pos_eci, eop_data, HarrisPriesterModified(; n=6))
    @test rho == ρ
end
