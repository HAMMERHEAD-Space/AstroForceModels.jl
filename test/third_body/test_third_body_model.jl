@testset "Third Body Model Ephemeris" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    frames = create_test_frames(eop_data)

    sun_model = test_sun_model()

    t_ft = (JD - 2451545.0) * 86400.0  # Convert JD to seconds since J2000

    sun_pos = get_position(sun_model.ephem_type, sun_model.body, frames, t_ft)

    # Old expected values were in meters; new API returns km, so divide by 1e3
    expected_sun_pos = [36327254721.23421, -130787280761.2391, -56694897390.09085] ./ 1e3

    @test expected_sun_pos ≈ sun_pos rtol = 1E-4

    sun_pos_vel, sun_vel = get_velocity(sun_model.ephem_type, sun_model.body, frames, t_ft)

    # Old expected values were in m/s; new API returns km/s, so divide by 1e3
    expected_sun_vel = [29356.8910247322, 6844.1181149359, 2965.5158876916] ./ 1e3

    @test expected_sun_vel ≈ sun_vel rtol = 1E-2

    moon_model = test_moon_model()

    moon_pos = get_position(moon_model.ephem_type, moon_model.body, frames, t_ft)

    # Old expected values were in meters; new API returns km, so divide by 1e3
    expected_moon_pos = [-344989050.2810446, -175734543.4089319, -81943013.49241804] ./ 1e3

    @test expected_moon_pos ≈ moon_pos rtol = 1E-3
end
