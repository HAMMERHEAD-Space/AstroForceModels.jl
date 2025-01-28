@testset "Drag Acceleration" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_drag_model = CannonballFixedDrag(0.2)

    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=JB2008(),
        eop_data=eop_data,
    )

    drag_accel = acceleration(state, p, 0.0, drag_model)

    # Generated with Orekit Harris-Priester
    expected_acceleration = [
        -4.987575495487041e-09, 1.8197638453380163e-09, 6.678597193592705e-10
    ] # km/s^2

    # There is no common atmosphere between the two so we just test direction
    @test dot(normalize(drag_accel), normalize(expected_acceleration)) ≈ 1.0
end
