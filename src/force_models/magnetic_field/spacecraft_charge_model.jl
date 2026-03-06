# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Spacecraft charge-to-mass ratio models for the Lorentz force computation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AbstractSpacecraftChargeModel, FixedChargeMassRatio, StateChargeModel
export charge_mass_ratio

"""
    AbstractSpacecraftChargeModel

Abstract type for spacecraft charge-to-mass ratio models.

Subtypes must implement `charge_mass_ratio(u, p, t, model)` returning q/m in C/kg.
"""
abstract type AbstractSpacecraftChargeModel end

"""
    FixedChargeMassRatio{QT<:Number} <: AbstractSpacecraftChargeModel

Fixed charge-to-mass ratio model. The user directly specifies q/m [C/kg].

Typical values range from 1e-6 to 1e-3 C/kg for natural spacecraft charging in
Low Earth Orbit, and up to ~0.03 C/kg for proposed active Lorentz-augmented orbits.

# Fields
- `q_over_m::QT`: Charge-to-mass ratio [C/kg].
"""
struct FixedChargeMassRatio{QT<:Number} <: AbstractSpacecraftChargeModel
    q_over_m::QT
end

"""
    StateChargeModel{FT<:Function} <: AbstractSpacecraftChargeModel

State-dependent charge-to-mass ratio model. The user provides a function
`f(u, p, t) -> q/m` that returns the charge-to-mass ratio [C/kg] as a function
of the current state, parameters, and time.

# Fields
- `charge_function::FT`: Function `(u, p, t) -> q/m` [C/kg].
"""
struct StateChargeModel{FT<:Function} <: AbstractSpacecraftChargeModel
    charge_function::FT
end

"""
    charge_mass_ratio(u, p, t, model::FixedChargeMassRatio)

Return the fixed charge-to-mass ratio [C/kg].
"""
@inline function charge_mass_ratio(
    u::AbstractVector, p::AbstractVector, t::Number, model::FixedChargeMassRatio
)
    return model.q_over_m
end

"""
    charge_mass_ratio(u, p, t, model::StateChargeModel)

Evaluate the state-dependent charge-to-mass ratio function [C/kg].
"""
@inline function charge_mass_ratio(
    u::AbstractVector, p::AbstractVector, t::Number, model::StateChargeModel
)
    return model.charge_function(u, p, t)
end
