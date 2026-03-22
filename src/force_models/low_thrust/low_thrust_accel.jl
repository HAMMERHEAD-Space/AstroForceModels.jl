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
    LowThrustAstroModel

Low-thrust propulsion force model. Frame-independent (RTN/VNB are orbital frames).

# Fields
- `thrust_model::AbstractThrustModel`: Thrust model defining the acceleration vector.
- `frame::AbstractThrustFrame`: Reference frame for thrust output. Default: `InertialFrame()`.
"""
Base.@kwdef struct LowThrustAstroModel{TM<:AbstractThrustModel,FR<:AbstractThrustFrame} <:
                   AbstractNonPotentialBasedForce
    thrust_model::TM
    frame::FR = InertialFrame()
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, lt_model::LowThrustAstroModel)

Computes the low-thrust acceleration in the inertial frame.

# Returns
- `SVector{3}`: Thrust acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, lt_model::LowThrustAstroModel
)
    a_local = thrust_acceleration(u, p, t, lt_model.thrust_model)
    return transform_thrust_to_state_frame(a_local, u, p, t, lt_model.frame)
end
