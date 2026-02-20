# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from ionospheric plasma (ion) drag.
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
    PlasmaDragAstroModel{ST,IT,EoT} <: AbstractNonPotentialBasedForce

Ionospheric plasma drag force model for spacecraft in low-Earth orbit.

Plasma drag arises from direct collection of ionospheric ions (primarily O⁺ in the 
F2 region) and their momentum transfer to the spacecraft. The acceleration follows 
the same cannonball drag formulation as atmospheric drag, but uses the ion mass 
density from the ionosphere rather than the neutral atmospheric density:

    𝐚 = -½ ⋅ BC_i ⋅ ρᵢ ⋅ |𝐯_app| ⋅ 𝐯_app

where BC_i = C_{D,i} A/m is the ion ballistic coefficient, ρᵢ is the ion mass density, 
and 𝐯_app is the spacecraft velocity relative to the co-rotating ionosphere [1, 3].

At LEO altitudes (300-600 km), plasma drag can contribute 5-35% of total aerodynamic 
force, with larger contributions at higher altitudes where neutral density decreases 
faster than ion density [1].

# Type Parameters
- `ST <: AbstractSatellitePlasmaDragModel`: Satellite plasma drag model type
- `IT <: AbstractIonosphereModel`: Ionospheric density model type
- `EoT <: Union{EopIau1980,EopIau2000A}`: Earth Orientation Parameters type

# Fields
- `satellite_plasma_drag_model::ST`: Model for computing ion ballistic coefficient
- `ionosphere_model::IT`: Model for computing ion mass density
- `eop_data::EoT`: Earth Orientation Parameters for coordinate transformations
"""
Base.@kwdef struct PlasmaDragAstroModel{
    ST<:AbstractSatellitePlasmaDragModel,
    IT<:AbstractIonosphereModel,
    EoT<:Union{EopIau1980,EopIau2000A},
} <: AbstractNonPotentialBasedForce
    satellite_plasma_drag_model::ST
    ionosphere_model::IT
    eop_data::EoT
end

"""
    acceleration(u, p, t, model::PlasmaDragAstroModel) -> SVector{3}

Compute the plasma drag acceleration on a spacecraft.

# Arguments
- `u::AbstractVector`: State vector [r; v] in J2000 ECI [km; km/s].
- `p::ComponentVector`: Parameters; must include `p.JD` (Julian Date at epoch).
- `t::Number`: Elapsed time since epoch [s].
- `model::PlasmaDragAstroModel`: Plasma drag force model.

# Returns
- `SVector{3}`: Plasma drag acceleration in J2000 ECI [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, model::PlasmaDragAstroModel
)
    rho_i = compute_ion_density(current_jd(p, t), u, model.eop_data, model.ionosphere_model)

    ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)

    BC_i = ion_ballistic_coefficient(u, p, t, model.satellite_plasma_drag_model)

    return plasma_drag_accel(u, rho_i, BC_i, ω_vec)
end

"""
    plasma_drag_accel(u, rho_i, BC_i, ω_vec) -> SVector{3}

Low-level computation of plasma drag acceleration.

The ionospheric plasma is treated as co-rotating with the Earth. The apparent 
velocity of the spacecraft relative to the plasma is:

    𝐯_app = 𝐯 - 𝛚 × 𝐫

and the plasma drag acceleration is:

    𝐚 = -½ ⋅ BC_i ⋅ ρᵢ ⋅ |𝐯_app| ⋅ 𝐯_app

The factor of 1E3 converts from m/s² to km/s² (density is in kg/m³, velocity in km/s).

# Arguments
- `u::AbstractVector`: State [r; v] in J2000 ECI [km; km/s].
- `rho_i::Number`: Ion mass density [kg/m³].
- `BC_i::Number`: Ion ballistic coefficient C_{D,i} A/m [m²/kg].
- `ω_vec::AbstractVector`: Earth angular velocity vector [rad/s].

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
