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
    velocity = SVector{3}(state[4], state[5], state[6])

    # --- JB2008 drag model ---
    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=JB2008(),
        eop_data=eop_data,
    )

    drag_accel = acceleration(state, p, 0.0, drag_model)

    drag_magnitude = norm(drag_accel)
    @test drag_magnitude > 1e-10
    @test drag_magnitude < 1e-7
    @test dot(drag_accel, velocity) < 0

    # --- Harris-Priester drag model ---
    # Reference generated with Python Orekit 13.1 Harris-Priester (test/orekit_hp_reference.py)
    orekit_hp_acceleration = [
        -1.4923697153110742e-08, 5.4450504968106981e-09, 1.9983709103332451e-09
    ] # km/s^2

    hp_drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=HarrisPriester(),
        eop_data=eop_data,
    )
    hp_accel = acceleration(state, p, 0.0, hp_drag_model)

    @test norm(hp_accel) ≈ norm(orekit_hp_acceleration) rtol = 1e-3
    @test dot(normalize(hp_accel), normalize(orekit_hp_acceleration)) ≈ 1.0 atol = 1e-6
end
