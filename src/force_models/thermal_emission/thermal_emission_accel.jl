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
    ThermalEmissionAstroModel

Spacecraft thermal emission (thermal re-radiation) force model.

# Fields
- `satellite_thermal_model`: Satellite thermal model providing the emission coefficient.
- `sun_data::ThirdBodyModel`: Model to compute the Sun's position via FrameEphemeris.
- `shadow_model::ShadowModelType`: Shadow model type — defaults to `Conical()`.
- `R_Sun::Number`: Radius of the Sun [km].
- `R_Occulting::Number`: Radius of the occulting body [km]. No default — must be specified.
- `Ψ::Number`: Solar flux constant at 1 AU [N/m²].
- `AU::Number`: Astronomical Unit [km].
"""
Base.@kwdef struct ThermalEmissionAstroModel{
    TT<:AbstractSatelliteThermalModel,
    SDT<:ThirdBodyModel,
    SMT<:ShadowModelType,
    RST<:Number,
    ROT<:Number,
    PT<:Number,
    AUT<:Number,
} <: AbstractNonPotentialBasedForce
    satellite_thermal_model::TT
    sun_data::SDT
    shadow_model::SMT = Conical()

    R_Sun::RST = R_SUN
    R_Occulting::ROT
    Ψ::PT = SOLAR_FLUX
    AU::AUT = ASTRONOMICAL_UNIT / 1E3
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, model::ThermalEmissionAstroModel)

Compute the acceleration from spacecraft thermal emission.

# Returns
- `SVector{3}`: Thermal emission acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, model::ThermalEmissionAstroModel
)
    sun_pos = get_position(
        model.sun_data.ephem_type,
        model.sun_data.body,
        p.frames,
        ft_time(p, t),
        model.sun_data.compiled_vector3,
    )

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

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `sun_pos::AbstractVector`: Current position of the Sun [km].
- `C_thm::Number`: Thermal emission coefficient [m²/kg].

# Keyword Arguments
- `ShadowModel::ShadowModelType`: Shadow model. Default: `Conical()`.
- `R_Sun::Number`: Radius of the Sun [km].
- `R_Occulting::Number`: Radius of the occulting body [km].
- `Ψ::Number`: Solar radiation pressure at 1 AU [N/m²].
- `AU::Number`: Astronomical Unit [km].

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
