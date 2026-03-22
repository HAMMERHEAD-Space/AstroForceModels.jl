# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Drag
#
#   Earth-only. Requires Earth atmosphere model and ITRF frame in FrameSystem.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#
#   [1] Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications (4th ed.). Microcosm Press.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export DragAstroModel, drag_accel

"""
    DragAstroModel

Atmospheric drag force model for spacecraft orbital dynamics.

Earth-only. Requires Earth atmosphere model and an ITRF frame registered in the FrameSystem.

# Constructor
    DragAstroModel(; satellite_drag_model, atmosphere_model, body_fixed_frame=:ITRF,
        propagation_frame=:ICRF, frames=nothing, rts=nothing, P=nothing)

# Arguments
- `satellite_drag_model`: Satellite drag model providing the ballistic coefficient.
- `atmosphere_model::AtmosphericModelType`: Atmospheric density model (JB2008, JR1971, etc.).
- `body_fixed_frame::Symbol`: Body-fixed frame for geodetic conversion (default: `:ITRF`).
- `propagation_frame::Symbol`: Frame in which the state vector is expressed (default: `:ICRF`).
- `frames`: Optional `FrameSystem`. When provided, pre-compiles the rotation between
  `propagation_frame` and `body_fixed_frame` for allocation-free evaluation.
- `rts`: Optional pre-allocated roots container for atmospheric models.
- `P`: Optional pre-allocated matrix for atmospheric models.
"""
struct DragAstroModel{
    ST<:AbstractSatelliteDragModel,
    AT<:AtmosphericModelType,
    RT<:Union{Nothing,AbstractVector},
    PT<:Union{Nothing,AbstractMatrix},
    CR3,
} <: AbstractNonPotentialBasedForce
    satellite_drag_model::ST
    atmosphere_model::AT
    body_fixed_frame::Symbol
    propagation_frame::Symbol
    rts::RT
    P::PT
    compiled_rotation3::CR3
end

function DragAstroModel(;
    satellite_drag_model::ST,
    atmosphere_model::AT,
    body_fixed_frame::Symbol=:ITRF,
    propagation_frame::Symbol=:ICRF,
    frames=nothing,
    rts::RT=nothing,
    P::PT=nothing,
    # Legacy: accept eop_data but ignore it (kept for backward compat during transition)
    eop_data=nothing,
) where {ST<:AbstractSatelliteDragModel,AT<:AtmosphericModelType,RT,PT}
    cr3 = nothing
    if !isnothing(frames) && propagation_frame != body_fixed_frame
        cr3 = compile_rotation3(frames, propagation_frame, body_fixed_frame)
    end
    return DragAstroModel(
        satellite_drag_model,
        atmosphere_model,
        body_fixed_frame,
        propagation_frame,
        rts,
        P,
        cr3,
    )
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, drag_model::DragAstroModel)

Computes the drag acceleration acting on a spacecraft.

Earth-only. Uses the FrameSystem to rotate the state into the body-fixed frame for
geodetic coordinate computation required by atmospheric density models.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `drag_model::DragAstroModel`: Drag model.

# Returns
- `SVector{3}`: Drag acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, drag_model::DragAstroModel
)
    t_ft = ft_time(p, t)

    # Get propagation→body-fixed rotation
    R_rot = if isnothing(drag_model.compiled_rotation3)
        rotation3(p.frames, drag_model.propagation_frame, drag_model.body_fixed_frame, t_ft)
    else
        drag_model.compiled_rotation3(t_ft)
    end
    R_prop2bf = R_rot.m[1]

    # Compute density at the satellite's current position
    rho = compute_density(
        current_jd(p, t),
        u,
        R_prop2bf,
        drag_model.atmosphere_model;
        roots_container=drag_model.rts,
        P=drag_model.P,
    )

    ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)

    # Compute the ballistic coefficient
    BC = ballistic_coefficient(u, p, t, drag_model.satellite_drag_model)

    # Return the 3-Dimensional Drag Force
    return drag_accel(u, rho, BC, ω_vec)
end

"""
    drag_accel(u::AbstractVector, rho::Number, BC::Number, ω_vec::AbstractVector) -> SVector{3}

Compute the Acceleration from Atmospheric Drag.

The atmosphere is treated as a solid revolving with the Earth and the apparent velocity of the satellite is computed
using the transport theorem

                𝐯_app = 𝐯 - 𝛚 x 𝐫

The acceleration from drag is then computed with a cannonball model as

                𝐚 = 1/2 * ρ * BC * |𝐯_app|₂^2 * v̂

# Arguments

- `u::AbstractVector`: The current state of the spacecraft [km, km/s].
- `rho::Number`: Atmospheric density at (t, u) [kg/m^3].
- `BC::Number`: The ballistic coefficient (area/mass) * drag coefficient [m^2/kg].
- `ω_vec::AbstractVector`: The angular velocity vector of Earth [rad/s].

# Returns

- `SVector{3}`: Inertial acceleration from drag [km/s²].
"""
@inline function drag_accel(
    u::AbstractVector{UT}, rho::RT, BC::BT, ω_vec::AbstractVector{WT}
) where {UT,RT,BT,WT}
    AT = promote_type(UT, RT, BT, WT)

    # Compute Apparent Velocity w.r.t the Atmosphere using the Transport Theorem
    apparent_vel = SVector{3}(u[4], u[5], u[6]) - cross(ω_vec, SVector{3}(u[1], u[2], u[3]))

    # Scaled by 1E3 to convert to km/s
    drag_force = -0.5 * BC * rho * norm(apparent_vel) * 1E3
    accel = SVector{3,AT}(
        drag_force * apparent_vel[1],
        drag_force * apparent_vel[2],
        drag_force * apparent_vel[3],
    )

    return accel
end
