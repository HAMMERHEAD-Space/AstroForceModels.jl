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
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, third_body_model::ThirdBodyModel)

Computes the drag acceleration acting on a spacecraft given a drag model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `third_body_model::ThirdBodyModel`: Third body model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional srp acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::TT, third_body::ThirdBodyModel
) where {TT}
    body_pos = third_body(p.JD + t / 86400.0, Position()) ./ 1E3

    return third_body_accel(u, third_body.body.μ, body_pos)
end

export third_body_accel

"""
    third_body_accel(u::AbstractVector, μ_body::Number, body_pos::AbstractVector, h::Number) -> SVector{3}{Number}

Compute the Acceleration from a 3rd Body Represented as a Point Mass

Since the central body is also being acted upon by the third body, the acceleration of body 𝐁 acting on 
spacecraft 𝐀 in the orbiting body's 𝐂 is part of the force not acting on the central body

                a = ∇UB(rA) - ∇UB(rC)

# Arguments

- `u::AbstractVector`: The current state of the spacecraft in the central body's inertial frame.
- `μ_body`: Gravitation Parameter of the 3rd body.
- `body_pos::AbstractVector`: The current position of the 3rd body in the central body's inertial frame.

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
@inline function third_body_accel(
    u::AbstractVector{PT}, μ_body::Number, body_pos::AbstractVector{BT}
) where {PT,BT}
    RT = promote_type(PT, BT)

    # Compute Position Vectors for the Spacecraft w.r.t the Central and 3rd Body Respectively
    sat_pos = @view(u[1:3])
    r_spacecraft_to_body = body_pos - sat_pos

    # Calculate and Return the Acceleration from the Difference in Potential
    r_sc_body = norm(r_spacecraft_to_body)
    r_body = norm(body_pos)

    return SVector{3,RT}(
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[1] - (μ_body / (r_body^3)) * body_pos[1],
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[2] - (μ_body / (r_body^3)) * body_pos[2],
        (μ_body / (r_sc_body^3)) * r_spacecraft_to_body[3] - (μ_body / (r_body^3)) * body_pos[3],
    )
end
