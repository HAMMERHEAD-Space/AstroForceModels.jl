# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from a 3rd Body Represented as a Point Mass
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/acceleration_models/third_body_acceleration.html
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, third_body_model::ThirdBodyModel)

Computes the third-body gravitational acceleration acting on a spacecraft.

# Arguments
- `u::AbstractVector`: Current state of the simulation.
- `p::FrameAwareParams`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `third_body_model::ThirdBodyModel`: Third body model struct.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional third-body acceleration [km/s²].
"""
function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::TT, third_body::ThirdBodyModel
) where {TT}
    body_pos = get_position(
        third_body.ephem_type,
        third_body.body,
        p.frames,
        ft_time(p, t),
        third_body.compiled_vector3,
    )

    return third_body_accel(u, third_body.body.μ, body_pos)
end

export third_body_accel

"""
    third_body_accel(u::AbstractVector, μ_body::Number, body_pos::AbstractVector) -> SVector{3}

Compute the Acceleration from a 3rd Body Represented as a Point Mass

Since the central body is also being acted upon by the third body, the acceleration of body 𝐁 acting on
spacecraft 𝐀 in the orbiting body's 𝐂 is part of the force not acting on the central body

                a = ∇UB(rA) - ∇UB(rC)

# Arguments

- `u::AbstractVector`: The current state of the spacecraft in the central body's inertial frame.
- `μ_body`: Gravitation Parameter of the 3rd body.
- `body_pos::AbstractVector`: The current position of the 3rd body in the central body's inertial frame [km].

# Returns

- `SVector{3}`: Inertial acceleration from the 3rd body [km/s²].
"""
@inline function third_body_accel(
    u::AbstractVector{PT}, μ_body::Number, body_pos::AbstractVector{BT}
) where {PT,BT}
    RT = promote_type(PT, BT)

    # Compute Position Vectors for the Spacecraft w.r.t the Central and 3rd Body Respectively
    sat_pos = SVector{3,PT}(u[1], u[2], u[3])
    r_spacecraft_to_body = body_pos - sat_pos

    # Calculate and Return the Acceleration from the Difference in Potential
    r_sc_body = norm(r_spacecraft_to_body)
    r_body = norm(body_pos)

    return SVector{3,RT}(
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[1] -
        (μ_body / (r_body^3)) * body_pos[1],
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[2] -
        (μ_body / (r_body^3)) * body_pos[2],
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[3] -
        (μ_body / (r_body^3)) * body_pos[3],
    )
end
