@testset "Ionosphere Models" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()

    # Compute the ECI→ECEF rotation DCM for the test epoch
    R_J2000_ITRF = r_eci_to_ecef(J2000(), ITRF(), JD, eop_data)

    # LEO state at ~400 km altitude
    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    @testset "ChapmanIonosphere" begin
        iono = ChapmanIonosphere()

        rho = compute_ion_density(JD, state, R_J2000_ITRF, iono)

        # Ion density should be positive at LEO altitudes
        @test rho > 0.0

        # Ion mass density at ~400 km: n_i ≈ 10^10-10^12 m⁻³, m_O+ ≈ 2.66e-26 kg
        # => ρ_i ≈ 10^(-16) to 10^(-14) kg/m³
        @test rho > 1e-20
        @test rho < 1e-12

        # Custom parameters (solar maximum)
        iono_max = ChapmanIonosphere(; Nmax=1e12, hmax=400.0, Hs=70.0)
        rho_max = compute_ion_density(JD, state, R_J2000_ITRF, iono_max)
        @test rho_max > 0.0

        # Higher Nmax should give higher density at similar altitudes
        iono_low = ChapmanIonosphere(; Nmax=1e10)
        rho_low = compute_ion_density(JD, state, R_J2000_ITRF, iono_low)
        @test rho_low < rho

        # Very high altitude (above 1500 km) should give zero
        high_state = [0.0, 0.0, 7878.0 + R_EARTH, 0.0, 0.0, 0.0]
        rho_high = compute_ion_density(JD, high_state, R_J2000_ITRF, iono)
        @test rho_high == 0.0
    end

    @testset "ConstantIonosphere" begin
        rho_fixed = 1e-17
        iono = ConstantIonosphere(; rho_i=rho_fixed)
        rho = compute_ion_density(JD, state, R_J2000_ITRF, iono)
        @test rho == rho_fixed
    end

    @testset "NoIonosphere" begin
        iono = NoIonosphere()
        rho = compute_ion_density(JD, state, R_J2000_ITRF, iono)
        @test rho == 0.0
    end
end
