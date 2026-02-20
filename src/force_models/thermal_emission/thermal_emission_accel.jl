# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Spacecraft Thermal Emission (Thermal Re-Radiation)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Duan, B. & Hugentobler, U. (2022). "Estimating surface optical properties and
#       thermal thrust for Galileo satellite body and solar panels." GPS Solutions, 26, 135.
#       https://doi.org/10.1007/s10291-022-01324-1
#
#   [2] Vigue, Y., Schutz, B. E. & Abusali, P. (1994). "Thermal force modeling for global
#       positioning system using the finite element method." Journal of Spacecraft and
#       Rockets, 31(5), 855-859.
#
#   [3] Milani, A., Nobili, A. M. & Farinella, P. (1987). "Non-gravitational Perturbations
#       and Satellite Geodesy." Adam Hilger, Bristol.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ThermalEmissionAstroModel, thermal_emission_accel

"""
    ThermalEmissionAstroModel{TT,SDT,EoT,SMT,RST,ROT,PT,AUT} <: AbstractNonPotentialBasedForce

Spacecraft thermal emission (thermal re-radiation) force model.

When a spacecraft absorbs solar radiation, the surface heats up and re-emits thermal
infrared radiation. If the emission is anisotropic (e.g., due to different emissivities on
front and back surfaces of solar panels), the net photon recoil produces a small but
persistent perturbative acceleration. This effect is most significant during eclipse season
transitions for GNSS satellites, where it can cause orbit errors of 1-3 cm [1].

In the steady-state (equilibrium) approximation, the thermal emission acceleration has
the same functional form as solar radiation pressure:

    𝐚 = ν * C_thm * Ψ * (AU / d_sun)² * d̂_sun

where ν is the shadow factor, C_thm is the thermal emission coefficient [m²/kg],
Ψ is the solar flux at 1 AU [N/m²], d_sun is the Sun-spacecraft distance, and d̂_sun
is the unit vector from the Sun to the spacecraft.

For a flat plate with front emissivity ε_f and back emissivity ε_b, the thermal emission
coefficient in equilibrium is [1,2]:

    C_thm = (2/3) * (A/M) * [(ε_f - ξ·ε_b) / (ε_f + ξ·ε_b)] * (1 - η) * α

# Type Parameters
- `TT <: AbstractSatelliteThermalModel`: Type of the satellite thermal model
- `SDT <: ThirdBodyModel`: Type of the Sun data model
- `EoT <: Union{EopIau1980,EopIau2000A}`: Type of Earth Orientation Parameters
- `SMT <: ShadowModelType`: Type of the shadow model

# Fields
- `satellite_thermal_model::TT`: Satellite thermal model providing the emission coefficient
- `sun_data::SDT`: Model to compute the Sun's position and velocity
- `eop_data::EoT`: Earth Orientation Parameters for coordinate transformations
- `shadow_model::SMT`: Shadow model type — defaults to `Conical()`
- `R_Sun::RST`: Radius of the Sun [km] — defaults to `R_SUN`
- `R_Occulting::ROT`: Radius of the occulting body [km] — defaults to `R_EARTH`
- `Ψ::PT`: Solar flux constant at 1 AU [N/m²] — defaults to `SOLAR_FLUX`
- `AU::AUT`: Astronomical Unit [km] — defaults to `ASTRONOMICAL_UNIT / 1E3`
"""
Base.@kwdef struct ThermalEmissionAstroModel{
    TT<:AbstractSatelliteThermalModel,
    SDT<:ThirdBodyModel,
    EoT<:Union{EopIau1980,EopIau2000A},
    SMT<:ShadowModelType,
    RST<:Number,
    ROT<:Number,
    PT<:Number,
    AUT<:Number,
} <: AbstractNonPotentialBasedForce
    satellite_thermal_model::TT
    sun_data::SDT
    eop_data::EoT
    shadow_model::SMT = Conical()

    R_Sun::RST = R_SUN
    R_Occulting::ROT = R_EARTH
    Ψ::PT = SOLAR_FLUX
    AU::AUT = ASTRONOMICAL_UNIT / 1E3
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, model::ThermalEmissionAstroModel)

Compute the acceleration from spacecraft thermal emission given a thermal model and the
current state and parameters of the spacecraft.

# Arguments
- `u::AbstractVector`: Current state of the simulation [km, km/s].
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation [s].
- `model::ThermalEmissionAstroModel`: Thermal emission model.

# Returns
- `SVector{3}`: The 3-dimensional thermal emission acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, model::ThermalEmissionAstroModel
)
    sun_pos = model.sun_data(current_jd(p, t), Position()) ./ 1E3

    C_thm = thermal_emission_coefficient(u, p, t, model.satellite_thermal_model)

    return thermal_emission_accel(
        u,
        sun_pos,
        C_thm;
        ShadowModel=model.shadow_model,
        R_Sun=model.R_Sun,
        R_Occulting=model.R_Occulting,
        Ψ=model.Ψ,
        AU=model.AU,
    )
end

"""
    thermal_emission_accel(u, sun_pos, C_thm; ShadowModel, R_Sun, R_Occulting, Ψ, AU)

Compute the acceleration from spacecraft thermal emission (thermal re-radiation).

Spacecraft surfaces absorb solar radiation, heat up, and re-emit thermal photons. When
emission is anisotropic (different emissivities on opposite sides), the net photon recoil
produces a perturbative force. In the steady-state equilibrium approximation, the
acceleration has the form [1]:

    𝐚 = ν * C_thm * Ψ * (AU / d_sun)² * d̂_sun

This is structurally similar to solar radiation pressure but captures the distinct physical
mechanism of absorbed-then-re-emitted thermal radiation rather than reflected/scattered
photons. The thermal emission coefficient C_thm encodes the surface emissivity asymmetry
and absorptivity properties of the spacecraft.

# Arguments
- `u::AbstractVector`: Current state of the spacecraft in the central body inertial frame [km, km/s].
- `sun_pos::AbstractVector`: Current position of the Sun [km].
- `C_thm::Number`: Thermal emission coefficient [m²/kg].

# Keyword Arguments
- `ShadowModel::ShadowModelType`: Shadow model. Default: `Conical()`.
- `R_Sun::Number`: Radius of the Sun [km]. Default: `R_SUN`.
- `R_Occulting::Number`: Radius of the occulting body [km]. Default: `R_EARTH`.
- `Ψ::Number`: Solar radiation pressure at 1 AU [N/m²]. Default: `SOLAR_FLUX`.
- `AU::Number`: Astronomical Unit [km]. Default: `ASTRONOMICAL_UNIT / 1E3`.

# Returns
- `SVector{3}`: Inertial acceleration from thermal emission [km/s²].
"""
@inline function thermal_emission_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector,
    C_thm::Number;
    ShadowModel::ShadowModelType=Conical(),
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    Ψ::Number=SOLAR_FLUX,
    AU::Number=ASTRONOMICAL_UNIT / 1E3,
) where {UT}
    sat_pos = SVector{3,UT}(u[1], u[2], u[3])

    F = shadow_model(sat_pos, sun_pos, ShadowModel; R_Sun=R_Sun, R_Occulting=R_Occulting)

    R_spacecraft_Sun = sat_pos - sun_pos
    R_sc_sun = norm(R_spacecraft_Sun)

    F_thm = F * C_thm * Ψ * (AU / R_sc_sun)^2 / 1E3

    return SVector{3}(
        F_thm * R_spacecraft_Sun[1] / R_sc_sun,
        F_thm * R_spacecraft_Sun[2] / R_sc_sun,
        F_thm * R_spacecraft_Sun[3] / R_sc_sun,
    )
end
