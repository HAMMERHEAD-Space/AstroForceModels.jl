@testset "Gravitational Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    central_body_accel =
        -(GravityModels.gravity_constant(grav_coeffs) / 1E9) / (norm(state[1:3])^3) *
        state[1:3]
    zonal_accel = acceleration(state, p, 0.0, grav_model) - central_body_accel

    expected_zonal_accel =
        [0.0019242559697225815, 0.011610900519485680, 0.0017986584026060844] ./ 1E3

    @test expected_zonal_accel ≈ zonal_accel rtol = 1E-4
end

@testset "Gravitational Potential" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    U = potential(state, p, 0.0, grav_model)

    #TODO: Get a better test value
    expected_U = -58.14048760520612

    @test expected_U ≈ U rtol = 1E-8
end

@testset "Gravitational Potential Time Derivative" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    U = potential_time_derivative(state, p, 0.0, grav_model)

    #TODO: Get a better test value
    expected_U = 3.6750762128561695e-8

    @test expected_U ≈ U rtol = 1E-8
end
