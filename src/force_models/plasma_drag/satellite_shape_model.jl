# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Satellite plasma drag models for computing the ion ballistic coefficient.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractSatellitePlasmaDragModel, CannonballFixedPlasmaDrag, StatePlasmaDragModel
export ion_ballistic_coefficient

"""
    AbstractSatellitePlasmaDragModel

Abstract type for satellite models used to compute ion ballistic coefficients for 
plasma drag.

Ion drag coefficients typically range from 2.0 to 4.0 depending on surface 
accommodation and electrostatic sheath effects, compared to 2.0-2.5 for neutral 
atmospheric drag.
"""
abstract type AbstractSatellitePlasmaDragModel end

"""
    CannonballFixedPlasmaDrag{RT,MT,DT,BT} <: AbstractSatellitePlasmaDragModel

Cannonball plasma drag model with a fixed ion ballistic coefficient.

The ion ballistic coefficient is computed as BC_i = C_D,i × A / m, where C_D,i 
is the ion drag coefficient. For a sphere in the free molecular flow regime with 
full surface accommodation, C_D,i ≈ 4.0 (Chapra 1961, Lafleur 2023).

# Fields
- `radius::RT`: Spacecraft radius [m].
- `mass::MT`: Spacecraft mass [kg].
- `drag_coeff::DT`: Ion drag coefficient (dimensionless). Typical range 2.0-4.0.
- `ballistic_coeff::BT`: Precomputed C_D,i × A / m [m²/kg].

# Example
```julia
# From physical properties (sphere r=1m, m=500kg, Cd_i=4.0)
sat = CannonballFixedPlasmaDrag(1.0, 500.0, 4.0)

# From precomputed ballistic coefficient
sat = CannonballFixedPlasmaDrag(0.025)
```
"""
struct CannonballFixedPlasmaDrag{RT<:Number,MT<:Number,DT<:Number,BT<:Number} <:
       AbstractSatellitePlasmaDragModel
    radius::RT
    mass::MT
    drag_coeff::DT
    ballistic_coeff::BT
end

"""
    CannonballFixedPlasmaDrag(ballistic_coeff::Number)

Construct a fixed ion ballistic coefficient plasma drag model.

# Arguments
- `ballistic_coeff::Number`: Precomputed C_D,i × A / m [m²/kg].
"""
function CannonballFixedPlasmaDrag(ballistic_coeff::Number)
    (ballistic_coeff < 0.0) &&
        throw(ArgumentError("Ion ballistic coefficient should be ≥ 0"))
    return CannonballFixedPlasmaDrag(1.0, 1.0, ballistic_coeff, ballistic_coeff)
end

"""
    CannonballFixedPlasmaDrag(radius::Number, mass::Number, drag_coeff::Number)

Construct a cannonball plasma drag model from physical properties.

The ion ballistic coefficient is BC_i = C_D,i × π r² / m.

# Arguments
- `radius::Number`: Spacecraft radius [m].
- `mass::Number`: Spacecraft mass [kg].
- `drag_coeff::Number`: Ion drag coefficient (dimensionless).
"""
function CannonballFixedPlasmaDrag(radius::Number, mass::Number, drag_coeff::Number)
    (drag_coeff < 0.0) && throw(ArgumentError("Ion drag coefficient should be ≥ 0"))
    (radius < 0.0) && throw(ArgumentError("Radius should be ≥ 0"))
    (mass < 0.0) && throw(ArgumentError("Mass should be ≥ 0"))

    area = π * radius^2.0
    return CannonballFixedPlasmaDrag(radius, mass, drag_coeff, drag_coeff * area / mass)
end

"""
    ion_ballistic_coefficient(u, p, t, model::CannonballFixedPlasmaDrag)

Returns the fixed ion ballistic coefficient C_D,i × A / m [m²/kg].
"""
@inline function ion_ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::CannonballFixedPlasmaDrag
)
    return model.ballistic_coeff
end

"""
    StatePlasmaDragModel <: AbstractSatellitePlasmaDragModel

Plasma drag model that reads the ion ballistic coefficient from the state vector.
Intended for estimation problems where the ballistic coefficient is a solve-for parameter.
"""
struct StatePlasmaDragModel{IT<:Integer} <: AbstractSatellitePlasmaDragModel
    state_index::IT
end

StatePlasmaDragModel() = StatePlasmaDragModel(7)

"""
    ion_ballistic_coefficient(u, p, t, model::StatePlasmaDragModel)

Returns the ion ballistic coefficient from the state vector at `model.state_index`.
"""
@inline function ion_ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::StatePlasmaDragModel
)
    return u[model.state_index]
end
