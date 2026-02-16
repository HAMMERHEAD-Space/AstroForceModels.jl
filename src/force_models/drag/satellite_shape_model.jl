# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Satellite Drag Models to easily compute the ballistic coefficient
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractSatelliteDragModel, CannonballFixedDrag, StateDragModel
export ballistic_coefficient

"""
Abstract satellite drag model to help compute drag forces.
"""
abstract type AbstractSatelliteDragModel end

"""
Cannonball Fixed Drag Model struct
Contains information to compute the ballistic coefficient of a cannonball drag model with a 
fixed ballistic coefficient.

# Fields
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `drag_coeff::Number`: The drag coefficient of the spacecraft.
- `ballistic_coeff::Number`: The fixed ballistic coefficient to use.
"""
struct CannonballFixedDrag{RT<:Number,MT<:Number,DT<:Number,BT<:Number} <:
       AbstractSatelliteDragModel
    radius::RT
    mass::MT
    drag_coeff::DT
    ballistic_coeff::BT
end

"""
    CannonballFixedDrag(ballistic_coeff::Number)

Constructor for a fixed ballistic coefficient drag model.

# Arguments
- `ballistic_coeff::Number`: The fixed ballistic coefficient to use.

# Returns 
- `drag_model::CannonballFixedDrag`: A fixed ballistic coefficient drag model.

"""
function CannonballFixedDrag(ballistic_coeff::Number)
    (ballistic_coeff < 0.0) && throw(ArgumentError("Ballistic coefficient should be ≥ 0"))
    return CannonballFixedDrag(1.0, 1.0, ballistic_coeff, ballistic_coeff)
end

"""
    CannonballFixedDrag(radius::Number, mass::Number, drag_coeff::Number)

Constructor for a fixed cannonball ballistic coefficient drag model.

The ballistic coefficient is computed with the following equation:

                BC = CD * area/mass

where area is the 2D projection of a sphere

                area = π * r^2

# Arguments
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `drag_coeff::Number`: The drag coefficient of the spacecraft.

# Returns 
- `drag_model::CannonballFixedDrag`: A fixed ballistic coefficient drag model.

"""
function CannonballFixedDrag(radius::Number, mass::Number, drag_coeff::Number)
    (drag_coeff < 0.0) && throw(ArgumentError("Drag coefficient should be ≥ 0"))
    (radius < 0.0) && throw(ArgumentError("Radius should be ≥ 0"))
    (mass < 0.0) && throw(ArgumentError("Mass should be ≥ 0"))

    area = π * radius^2.0

    return CannonballFixedDrag(radius, mass, drag_coeff, drag_coeff * area / mass)
end

"""
ballistic_coefficient(
    u::AbstractVector, 
    p::ComponentVector, 
    t::Number, 
    model::CannonballFixedDrag)

Returns the ballistic coefficient for a drag model given the model and current state of the simulation.

# Arguments
    - `u::AbstractVector`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::CannonballFixedDrag`: Drag model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.

"""
@inline function ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::CannonballFixedDrag
)
    return model.ballistic_coeff
end

"""
State Drag Model struct
Empty struct used for when the simulation state includes the ballistic coefficient.
"""
struct StateDragModel <: AbstractSatelliteDragModel end

"""
ballistic_coefficient(
    u::AbstractVector, 
    p::ComponentVector, 
    t::Number, 
    model::StateDragModel)

Returns the ballistic coefficient for a drag model given the model and current state 
of the simulation.

# Arguments
    - `u::AbstractVector`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::StateDragModel`: The drag model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.
"""
@inline function ballistic_coefficient(
    u::AbstractVector, p::AbstractVector, t::Number, model::StateDragModel
)
    #TODO: GENERALIZE THE INDEX, COMPONENT VECTOR?
    return u[7]
end
