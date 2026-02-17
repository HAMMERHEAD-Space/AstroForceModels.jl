# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Thrust Models for Low-Thrust Propulsion
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
#   [4] Sims, J. A. and Flanagan, S. N. (1999). "Preliminary Design of Low-Thrust
#       Interplanetary Missions." AAS/AIAA Astrodynamics Specialist Conference, AAS 99-338.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AbstractThrustModel,
    ConstantCartesianThrust,
    ConstantTangentialThrust,
    StateThrustModel,
    PiecewiseConstantThrust
export thrust_acceleration

"""
    AbstractThrustModel

Abstract base type for all low-thrust propulsion models.

Concrete subtypes define how the thrust acceleration vector is computed from the spacecraft
state, simulation parameters, and time. Each subtype must implement:

    thrust_acceleration(u, p, t, model) → SVector{3}

returning a 3-component thrust acceleration [km/s²]. The components are interpreted in the
reference frame specified by the parent [`LowThrustAstroModel`](@ref)'s `frame` field (see
[`AbstractThrustFrame`](@ref)).

!!! note
    [`ConstantTangentialThrust`](@ref) is an exception — it always returns the acceleration 
    in the inertial frame and should be used with [`InertialFrame`](@ref) (the default).
"""
abstract type AbstractThrustModel end

# ======================================================================================== #
#                            ConstantCartesianThrust
# ======================================================================================== #

"""
    ConstantCartesianThrust{T1,T2,T3} <: AbstractThrustModel

Constant thrust acceleration expressed as a fixed 3-component vector.

The components are interpreted in the reference frame set on the parent
[`LowThrustAstroModel`](@ref):
- With [`InertialFrame`](@ref): `(ax, ay, az)` are ECI components
- With [`RTNFrame`](@ref): `(ax, ay, az)` are `(a_R, a_T, a_N)` components
- With [`VNBFrame`](@ref): `(ax, ay, az)` are `(a_V, a_N, a_B)` components

# Type Parameters
- `T1 <: Number`: Type of the first-component acceleration
- `T2 <: Number`: Type of the second-component acceleration
- `T3 <: Number`: Type of the third-component acceleration

# Fields
- `ax::T1`: First-component thrust acceleration [km/s²]
- `ay::T2`: Second-component thrust acceleration [km/s²]
- `az::T3`: Third-component thrust acceleration [km/s²]

# Constructors

    ConstantCartesianThrust(ax, ay, az)

Create from acceleration components directly [km/s²].

    ConstantCartesianThrust(direction, thrust, mass)

Create from a thrust direction vector (will be normalized), thrust magnitude [N], and 
spacecraft mass [kg]. The acceleration is computed as:

        a = (T / m) / 1000 * d̂   [km/s²]
"""
struct ConstantCartesianThrust{T1<:Number,T2<:Number,T3<:Number} <: AbstractThrustModel
    ax::T1
    ay::T2
    az::T3
end

"""
    ConstantCartesianThrust(direction::AbstractVector, thrust::Number, mass::Number)

Construct a constant Cartesian thrust model from a direction vector, thrust [N], and mass [kg].
"""
function ConstantCartesianThrust(direction::AbstractVector, thrust::Number, mass::Number)
    (thrust < 0.0) && throw(ArgumentError("Thrust magnitude should be ≥ 0"))
    (mass ≤ 0.0) && throw(ArgumentError("Mass should be > 0"))
    d = normalize(direction)
    accel_mag = thrust / (mass * 1.0E3)
    return ConstantCartesianThrust(accel_mag * d[1], accel_mag * d[2], accel_mag * d[3])
end

"""
    thrust_acceleration(u::AbstractVector, p::AbstractVector, t::Number, model::ConstantCartesianThrust)

Returns the constant Cartesian thrust acceleration vector [km/s²].

# Arguments
- `u::AbstractVector`: Current state of the simulation.
- `p::AbstractVector`: Parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `model::ConstantCartesianThrust`: Constant Cartesian thrust model.

# Returns
- `SVector{3}`: Thrust acceleration [km/s²] in the frame specified by the parent model.
"""
@inline function thrust_acceleration(
    u::AbstractVector, p::AbstractVector, t::Number, model::ConstantCartesianThrust
)
    return SVector{3}(model.ax, model.ay, model.az)
end

# ======================================================================================== #
#                           ConstantTangentialThrust
# ======================================================================================== #

"""
    ConstantTangentialThrust{MT<:Number} <: AbstractThrustModel

Constant-magnitude thrust directed along (or opposite to) the velocity vector.

When the magnitude is positive, thrust is aligned with the velocity direction (orbit raising).
When negative, thrust opposes the velocity direction (orbit lowering).

!!! note
    This model always computes the velocity-aligned direction internally and returns the
    acceleration in the **inertial** frame. It should be used with the default
    [`InertialFrame`](@ref). For velocity-aligned thrust in a different frame, use
    `ConstantCartesianThrust(magnitude, 0, 0)` with [`VNBFrame`](@ref).

# Type Parameters
- `MT <: Number`: Type of the acceleration magnitude

# Fields
- `magnitude::MT`: Thrust acceleration magnitude [km/s²]. Positive = along velocity, negative = anti-velocity.

# Constructors

    ConstantTangentialThrust(magnitude)

Create from acceleration magnitude directly [km/s²].

    ConstantTangentialThrust(thrust, mass)

Create from thrust [N] and spacecraft mass [kg]. The acceleration magnitude is computed as:

        |a| = T / (m × 1000)   [km/s²]
"""
struct ConstantTangentialThrust{MT<:Number} <: AbstractThrustModel
    magnitude::MT
end

"""
    ConstantTangentialThrust(thrust::Number, mass::Number)

Construct a constant tangential thrust model from thrust [N] and mass [kg].
"""
function ConstantTangentialThrust(thrust::Number, mass::Number)
    (mass ≤ 0.0) && throw(ArgumentError("Mass should be > 0"))
    return ConstantTangentialThrust(thrust / (mass * 1.0E3))
end

"""
    thrust_acceleration(u::AbstractVector, p::AbstractVector, t::Number, model::ConstantTangentialThrust)

Returns the tangential thrust acceleration vector directed along the velocity [km/s²].

The result is always in the inertial frame regardless of the parent model's frame setting.

# Arguments
- `u::AbstractVector`: Current state of the simulation.
- `p::AbstractVector`: Parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `model::ConstantTangentialThrust`: Constant tangential thrust model.

# Returns
- `SVector{3}`: Thrust acceleration in the inertial frame [km/s²].
"""
@inline function thrust_acceleration(
    u::AbstractVector{UT}, p::AbstractVector, t::Number, model::ConstantTangentialThrust{MT}
) where {UT,MT}
    RT = promote_type(UT, MT)

    v = SVector{3,UT}(u[4], u[5], u[6])
    v_norm = norm(v)

    if v_norm < eps(RT)
        z = zero(RT)
        return SVector{3}(z, z, z)
    end

    scale = model.magnitude / v_norm
    return SVector{3,RT}(scale * v[1], scale * v[2], scale * v[3])
end

# ======================================================================================== #
#                              StateThrustModel
# ======================================================================================== #

"""
    StateThrustModel{F} <: AbstractThrustModel

A user-defined thrust model where the thrust acceleration vector is computed by a 
provided function.

This is the most flexible thrust model, allowing arbitrary state- and time-dependent 
thrust profiles such as feedback control laws, scheduled thrust arcs, or optimization-based 
thrust histories.

The 3-component return value is interpreted in the frame specified by the parent
[`LowThrustAstroModel`](@ref).

# Type Parameters
- `F`: Type of the callable (function, functor, or closure)

# Fields
- `f::F`: Callable with signature `f(u, p, t) → SVector{3}` returning thrust acceleration [km/s²]
"""
struct StateThrustModel{F} <: AbstractThrustModel
    f::F
end

"""
    thrust_acceleration(u::AbstractVector, p::AbstractVector, t::Number, model::StateThrustModel)

Returns the thrust acceleration vector computed by the user-provided function [km/s²].

# Arguments
- `u::AbstractVector`: Current state of the simulation.
- `p::AbstractVector`: Parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `model::StateThrustModel`: User-defined thrust model.

# Returns
- `SVector{3}`: Thrust acceleration [km/s²] in the frame specified by the parent model.
"""
@inline function thrust_acceleration(
    u::AbstractVector, p::AbstractVector, t::Number, model::StateThrustModel
)
    return SVector{3}(model.f(u, p, t))
end

# ======================================================================================== #
#                          PiecewiseConstantThrust
# ======================================================================================== #

"""
    PiecewiseConstantThrust{N,TT<:Number,AT<:Number} <: AbstractThrustModel

A piecewise-constant thrust schedule defined by a sequence of time-tagged arcs.

Each arc specifies a constant 3-component acceleration vector that is active from its
start time until the next arc begins. This is the natural representation for
Sims-Flanagan trajectory segments, finite-burn maneuver schedules, and any thrust
profile that changes discretely between arcs.

The 3-component vectors are interpreted in the frame specified by the parent
[`LowThrustAstroModel`](@ref).

# Type Parameters
- `N`: Number of arcs (compile-time constant for type stability)
- `TT <: Number`: Element type of the time breakpoints
- `AT <: Number`: Element type of the acceleration components

# Fields
- `times::SVector{N,TT}`: Sorted arc start times [s]. Must be strictly increasing.
- `accelerations::SVector{N,SVector{3,AT}}`: Acceleration vector for each arc [km/s²].

Before `times[1]` and after the last arc, the last arc's acceleration is held.

# Constructors

    PiecewiseConstantThrust(times, accelerations)

`times` is an `AbstractVector` of arc start times and `accelerations` is an 
`AbstractVector` of 3-component acceleration vectors (or `SVector{3}`).
"""
struct PiecewiseConstantThrust{N,TT<:Number,AT<:Number} <: AbstractThrustModel
    times::SVector{N,TT}
    accelerations::SVector{N,SVector{3,AT}}
end

function PiecewiseConstantThrust(
    times::AbstractVector{TT}, accelerations::AbstractVector{<:AbstractVector{AT}}
) where {TT<:Number,AT<:Number}
    N = length(times)
    N == length(accelerations) || throw(
        ArgumentError(
            "times and accelerations must have the same length (got $N and $(length(accelerations)))",
        ),
    )
    N ≥ 1 || throw(ArgumentError("At least one arc is required"))

    for i in 2:N
        times[i] > times[i - 1] || throw(
            ArgumentError(
                "Arc start times must be strictly increasing (times[$i] ≤ times[$(i-1)])",
            ),
        )
    end

    sv_times = SVector{N,TT}(times...)
    sv_accels = SVector{N,SVector{3,AT}}(SVector{3,AT}(a...) for a in accelerations)
    return PiecewiseConstantThrust{N,TT,AT}(sv_times, sv_accels)
end

"""
    thrust_acceleration(u::AbstractVector, p::AbstractVector, t::Number, model::PiecewiseConstantThrust)

Returns the piecewise-constant thrust acceleration for the arc containing time `t` [km/s²].

Uses a reverse linear scan to find the active arc (last arc whose start time ≤ `t`).
If `t` is before the first arc, the first arc's acceleration is returned.

# Arguments
- `u::AbstractVector`: Current state of the simulation.
- `p::AbstractVector`: Parameters of the simulation.
- `t::Number`: Current time of the simulation [s].
- `model::PiecewiseConstantThrust`: Piecewise-constant thrust schedule.

# Returns
- `SVector{3}`: Thrust acceleration [km/s²] in the frame specified by the parent model.
"""
@inline function thrust_acceleration(
    u::AbstractVector, p::AbstractVector, t::Number, model::PiecewiseConstantThrust{N}
) where {N}
    idx = 1
    for i in 2:N
        if t ≥ model.times[i]
            idx = i
        end
    end
    return model.accelerations[idx]
end
