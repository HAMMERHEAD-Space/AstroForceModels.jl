# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Satellite Thermal Emission Models to compute the thermal emission coefficient
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AbstractSatelliteThermalModel, FixedThermalEmission, FlatPlateThermalModel
export thermal_emission_coefficient

"""
Abstract satellite thermal emission model.
"""
abstract type AbstractSatelliteThermalModel end

"""
    FixedThermalEmission{CT<:Number} <: AbstractSatelliteThermalModel

Fixed thermal emission coefficient model. The user directly specifies the
effective thermal emission coefficient C_thm [m²/kg].

The thermal emission acceleration magnitude is computed as:

    |a| = ν * C_thm * Ψ * (AU / d_sun)²

where ν is the shadow factor, Ψ is the solar radiation pressure at 1 AU,
and d_sun is the Sun-spacecraft distance.

# Fields
- `thermal_emission_coeff::Number`: The thermal emission coefficient [m²/kg].
"""
struct FixedThermalEmission{CT<:Number} <: AbstractSatelliteThermalModel
    thermal_emission_coeff::CT
end

"""
    FlatPlateThermalModel{AT,MT,EFT,EBT,ABT,XT,ET,CT} <: AbstractSatelliteThermalModel

Flat plate thermal emission model based on the solar panel thermal imbalance
formulation of Duan & Hugentobler (2022) and Vigue et al. (1994).

Computes the thermal emission coefficient from physical surface properties using
the equilibrium temperature assumption:

    C_thm = (2/3) * (A/M) * [(ε_f - ξ·ε_b) / (ε_f + ξ·ε_b)] * (1 - η) * α

The resulting acceleration is directed along the Sun-spacecraft line.

# Fields
- `area::Number`: Emitting surface area [m²]
- `mass::Number`: Spacecraft mass [kg]
- `ε_front::Number`: Front-side (Sun-facing) thermal emissivity
- `ε_back::Number`: Back-side thermal emissivity
- `absorptivity::Number`: Surface absorptivity for solar radiation
- `ξ::Number`: Temperature ratio factor between front and back sides (default: 1.0)
- `η::Number`: Electric conversion efficiency of solar panels (default: 0.0)
- `thermal_emission_coeff::Number`: Computed thermal emission coefficient [m²/kg]
"""
struct FlatPlateThermalModel{
    AT<:Number,
    MT<:Number,
    EFT<:Number,
    EBT<:Number,
    ABT<:Number,
    XT<:Number,
    ET<:Number,
    CT<:Number,
} <: AbstractSatelliteThermalModel
    area::AT
    mass::MT
    ε_front::EFT
    ε_back::EBT
    absorptivity::ABT
    ξ::XT
    η::ET
    thermal_emission_coeff::CT
end

"""
    FlatPlateThermalModel(; area, mass, ε_front, ε_back, absorptivity, ξ=1.0, η=0.0)

Construct a flat plate thermal model from physical surface properties.

The thermal emission coefficient is computed as:

    C_thm = (2/3) * (A/M) * [(ε_f - ξ·ε_b) / (ε_f + ξ·ε_b)] * (1 - η) * α

# Keyword Arguments
- `area::Number`: Emitting surface area [m²]
- `mass::Number`: Spacecraft mass [kg]
- `ε_front::Number`: Front-side (Sun-facing) thermal emissivity (0-1)
- `ε_back::Number`: Back-side thermal emissivity (0-1)
- `absorptivity::Number`: Surface absorptivity for solar radiation (0-1)
- `ξ::Number`: Temperature ratio factor between front and back sides (default: 1.0)
- `η::Number`: Electric conversion efficiency of solar panels (default: 0.0)
"""
function FlatPlateThermalModel(;
    area::Number,
    mass::Number,
    ε_front::Number,
    ε_back::Number,
    absorptivity::Number,
    ξ::Number=1.0,
    η::Number=0.0,
)
    (area < 0) && throw(ArgumentError("area must be ≥ 0"))
    (mass ≤ 0) && throw(ArgumentError("mass must be > 0"))
    (ε_front < 0 || ε_front > 1) && throw(ArgumentError("ε_front must be in [0, 1]"))
    (ε_back < 0 || ε_back > 1) && throw(ArgumentError("ε_back must be in [0, 1]"))
    (absorptivity < 0 || absorptivity > 1) &&
        throw(ArgumentError("absorptivity must be in [0, 1]"))
    (η < 0 || η > 1) && throw(ArgumentError("η must be in [0, 1]"))

    C_thm =
        (2 / 3) *
        (area / mass) *
        ((ε_front - ξ * ε_back) / (ε_front + ξ * ε_back)) *
        (1 - η) *
        absorptivity

    return FlatPlateThermalModel(area, mass, ε_front, ε_back, absorptivity, ξ, η, C_thm)
end

"""
    thermal_emission_coefficient(u, p, t, model::FixedThermalEmission)

Return the thermal emission coefficient for a fixed thermal emission model.
"""
@inline function thermal_emission_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::FixedThermalEmission
)
    return model.thermal_emission_coeff
end

"""
    thermal_emission_coefficient(u, p, t, model::FlatPlateThermalModel)

Return the pre-computed thermal emission coefficient for a flat plate model.
"""
@inline function thermal_emission_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::FlatPlateThermalModel
)
    return model.thermal_emission_coeff
end
