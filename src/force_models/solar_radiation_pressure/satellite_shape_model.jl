# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Satellite SRP Models to easily compute the reflectivity ballistic coefficient
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractSatelliteSRPModel, CannonballFixedSRP, StateSRPModel
export reflectivity_ballistic_coefficient

"""
Abstract satellite SRP model to help compute solar radiation pressure forces.
"""
abstract type AbstractSatelliteSRPModel end

"""
Cannonball Fixed SRP Model struct
Contains information to compute the reflectivity coefficient of a cannonball drag model with a 
fixed reflectivity coefficient.

# Fields
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `reflectivity_coeff::Number`: The reflectivity coefficient of the spacecraft.
- `reflectivity_ballistic_coeff::Number`: The fixed ballistic coefficient to use.
"""
struct CannonballFixedSRP{RT<:Number,MT<:Number,RcT<:Number,RbT<:Number} <:
       AbstractSatelliteSRPModel
    radius::RT
    mass::MT
    reflectivity_coeff::RcT
    reflectivity_ballistic_coeff::RbT
end

"""
    CannonballFixedSRP(ballistic_coeff::Number)

Constructor for a fixed ballistic coefficient SRP model.

# Arguments
- `ballistic_coeff::Number`: The fixed ballistic coefficient to use.

# Returns 
- `srp_model::CannonballFixedDrag`: A fixed ballistic coefficient SRP model.

"""
function CannonballFixedSRP(reflectivity_ballistic_coeff::Number)
    (reflectivity_ballistic_coeff < 0.0) &&
        throw(ArgumentError("Reflectivity ballistic coefficient should be ≥ 0"))

    return CannonballFixedSRP(
        1.0, 1.0, reflectivity_ballistic_coeff, reflectivity_ballistic_coeff
    )
end

"""
    CannonballFixedSRP(radius::Number, mass::Number, drag_coeff::Number)

Constructor for a fixed cannonball ballistic coefficient SRP model.

The ballistic coefficient is computed with the following equation:

                RC = CD * area/mass

where area is the 2D projection of a sphere

                area = π * r^2

# Arguments
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `reflectivity_coeff::Number`: The reflectivity coefficient of the spacecraft.

# Returns 
- srp_model::CannonballFixedSRP`: A fixed ballistic coefficient SRP model.

"""
function CannonballFixedSRP(radius::Number, mass::Number, reflectivity_coeff::Number)
    (reflectivity_coeff < 0.0) &&
        throw(ArgumentError("Reflectivity coefficient should be ≥ 0"))
    (radius < 0.0) && throw(ArgumentError("Radius should be ≥ 0"))
    (mass < 0.0) && throw(ArgumentError("Mass should be ≥ 0"))

    area = π * radius^2.0

    return CannonballFixedSRP(
        radius, mass, reflectivity_coeff, reflectivity_coeff * area / mass
    )
end

"""
reflectivity_ballistic_coefficient(
    u::AbstractVector, 
    p::ComponentVector, 
    t::Number, 
    model::CannonballFixedSRP)

Returns the ballistic coefficient for a SRP model given the model and current state of the simulation.

# Arguments
    - `u::AbstractVector`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::CannonballFixedSRP`: SRP model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.

"""
@inline function reflectivity_ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::CannonballFixedSRP
)
    return model.reflectivity_ballistic_coeff
end

"""
State SRP Model struct
Empty struct used for when the simulation state includes the reflectivity ballistic coefficient.
"""
struct StateSRPModel <: AbstractSatelliteSRPModel end

"""
ballistic_coefficient(
    u::AbstractVector, 
    p::ComponentVector, 
    t::Number, 
    model::StateSRPModel)

Returns the ballistic coefficient for a SRP model given the model and current state 
of the simulation.

# Arguments
    - `u::AbstractVector`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::StateSRPModel`: The SRP model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.
"""
@inline function reflectivity_ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::StateSRPModel
)
    #TODO: GENERALIZE THE INDEX, COMPONENT VECTOR?
    return u[8]
end
