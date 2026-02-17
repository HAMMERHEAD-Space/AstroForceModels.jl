# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Low-Thrust Propulsion
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications (4th ed.). 
#       Microcosm Press.
#   [2] Conway, B. A. (Ed.) (2010). Spacecraft Trajectory Optimization. Cambridge 
#       University Press.
#   [3] Betts, J. T. (2010). Practical Methods for Optimal Control and Estimation Using 
#       Nonlinear Programming (2nd ed.). SIAM.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export LowThrustAstroModel

"""
    LowThrustAstroModel{TM<:AbstractThrustModel,FR<:AbstractThrustFrame} <: AbstractNonPotentialBasedForce

Low-thrust propulsion force model for spacecraft orbital dynamics.

This model computes the acceleration due to low-thrust propulsion (e.g., electric/ion 
engines, solar sails, or any continuous thrust source). The thrust acceleration is 
determined by the underlying [`AbstractThrustModel`](@ref), which defines the magnitude 
and direction of the thrust vector, and the [`AbstractThrustFrame`](@ref), which specifies
how the thrust components are oriented.

# Type Parameters
- `TM <: AbstractThrustModel`: Type of the thrust model
- `FR <: AbstractThrustFrame`: Type of the reference frame

# Fields
- `thrust_model::TM`: The thrust model defining the acceleration vector (e.g., 
  [`ConstantCartesianThrust`](@ref), [`ConstantTangentialThrust`](@ref), 
  [`StateThrustModel`](@ref), [`PiecewiseConstantThrust`](@ref))
- `frame::FR`: The reference frame in which the thrust model's output is expressed.
  Defaults to [`InertialFrame`](@ref). Use [`RTNFrame`](@ref) or [`VNBFrame`](@ref) to
  specify thrust in orbital coordinates.

# Constructors

    LowThrustAstroModel(; thrust_model, frame=InertialFrame())
"""
Base.@kwdef struct LowThrustAstroModel{TM<:AbstractThrustModel,FR<:AbstractThrustFrame} <:
                   AbstractNonPotentialBasedForce
    thrust_model::TM
    frame::FR = InertialFrame()
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, lt_model::LowThrustAstroModel)

Computes the low-thrust acceleration acting on a spacecraft in the inertial frame.

The thrust model first produces a 3-component acceleration in the model's reference
frame, which is then transformed to the inertial frame using
[`transform_to_inertial`](@ref).

# Arguments
- `u::AbstractVector`: Current state of the simulation [x, y, z, vx, vy, vz] in km and km/s.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation [s].
- `lt_model::LowThrustAstroModel`: Low-thrust model containing the thrust configuration.

# Returns
- `SVector{3}`: The 3-dimensional thrust acceleration in the inertial frame [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, lt_model::LowThrustAstroModel
)
    a_local = thrust_acceleration(u, p, t, lt_model.thrust_model)
    return transform_to_inertial(a_local, u, lt_model.frame)
end
