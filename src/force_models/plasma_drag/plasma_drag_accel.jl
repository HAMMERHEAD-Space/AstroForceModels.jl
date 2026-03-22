# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from ionospheric plasma (ion) drag.
#
#   Earth-only. Requires Earth ionosphere model and ITRF frame in FrameSystem.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Lafleur, T. (2023). "Charged aerodynamics: Ionospheric plasma drag on objects in
#       low-Earth orbit." Acta Astronautica, 212, 370-386.
#
#   [2] Li, L.-S. (2011). "Perturbation effect of the Coulomb drag on the orbital
#       elements of the Earth satellite in the ionosphere." Acta Astronautica, 68, 717-721.
#
#   [3] Chapra, K.P. (1961). "Interaction of rapidly moving bodies on terrestrial
#       atmosphere." Rev. Mod. Phys., 33, 152-198.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export PlasmaDragAstroModel, plasma_drag_accel

"""
    PlasmaDragAstroModel

Ionospheric plasma drag force model for spacecraft in low-Earth orbit.

Earth-only. Requires Earth ionosphere model and ITRF frame in the FrameSystem.

# Constructor
    PlasmaDragAstroModel(; satellite_plasma_drag_model, ionosphere_model,
        body_fixed_frame=:ITRF, propagation_frame=:ICRF, frames=nothing)

# Arguments
- `satellite_plasma_drag_model`: Model for ion ballistic coefficient.
- `ionosphere_model`: Model for ion mass density.
- `body_fixed_frame::Symbol`: Body-fixed frame for geodetic conversion (default: `:ITRF`).
- `propagation_frame::Symbol`: Frame in which the state vector is expressed (default: `:ICRF`).
- `frames`: Optional `FrameSystem`. When provided, pre-compiles the rotation between
  `propagation_frame` and `body_fixed_frame` for allocation-free evaluation.
"""
struct PlasmaDragAstroModel{
    ST<:AbstractSatellitePlasmaDragModel,IT<:AbstractIonosphereModel,CR3
} <: AbstractNonPotentialBasedForce
    satellite_plasma_drag_model::ST
    ionosphere_model::IT
    body_fixed_frame::Symbol
    propagation_frame::Symbol
    compiled_rotation3::CR3
end

function PlasmaDragAstroModel(;
    satellite_plasma_drag_model::ST,
    ionosphere_model::IT,
    body_fixed_frame::Symbol=:ITRF,
    propagation_frame::Symbol=:ICRF,
    frames=nothing,
    # Legacy: accept eop_data but ignore it
    eop_data=nothing,
) where {ST<:AbstractSatellitePlasmaDragModel,IT<:AbstractIonosphereModel}
    cr3 = nothing
    if !isnothing(frames) && propagation_frame != body_fixed_frame
        cr3 = compile_rotation3(frames, propagation_frame, body_fixed_frame)
    end
    return PlasmaDragAstroModel(
        satellite_plasma_drag_model,
        ionosphere_model,
        body_fixed_frame,
        propagation_frame,
        cr3,
    )
end

"""
    acceleration(u, p, t, model::PlasmaDragAstroModel) -> SVector{3}

Compute the plasma drag acceleration on a spacecraft.

Earth-only.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `model::PlasmaDragAstroModel`: Plasma drag model.

# Returns
- `SVector{3}`: Plasma drag acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, model::PlasmaDragAstroModel
)
    t_ft = ft_time(p, t)

    R_rot = if isnothing(model.compiled_rotation3)
        rotation3(p.frames, model.propagation_frame, model.body_fixed_frame, t_ft)
    else
        model.compiled_rotation3(t_ft)
    end
    R_prop2bf = R_rot.m[1]

    rho_i = compute_ion_density(current_jd(p, t), u, R_prop2bf, model.ionosphere_model)

    ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)

    BC_i = ion_ballistic_coefficient(u, p, t, model.satellite_plasma_drag_model)

    return plasma_drag_accel(u, rho_i, BC_i, ω_vec)
end

"""
    plasma_drag_accel(u, rho_i, BC_i, ω_vec) -> SVector{3}

Low-level computation of plasma drag acceleration.

# Returns
- `SVector{3}`: Plasma drag acceleration [km/s²].
"""
@inline function plasma_drag_accel(
    u::AbstractVector{UT}, rho_i::RT, BC_i::BT, ω_vec::AbstractVector{WT}
) where {UT,RT,BT,WT}
    AT = promote_type(UT, RT, BT, WT)

    apparent_vel = SVector{3}(u[4], u[5], u[6]) - cross(ω_vec, SVector{3}(u[1], u[2], u[3]))

    force_mag = -0.5 * BC_i * rho_i * norm(apparent_vel) * 1E3
    accel = SVector{3,AT}(
        force_mag * apparent_vel[1],
        force_mag * apparent_vel[2],
        force_mag * apparent_vel[3],
    )

    return accel
end
