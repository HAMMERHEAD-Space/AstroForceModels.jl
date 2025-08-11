# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Full Dynamics Model and acceleration interface
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export CentralBodyDynamicsModel, build_dynamics_model

"""
    CentralBodyDynamicsModel{N,GT,PT} <: AbstractDynamicsModel

A dynamics model structure that efficiently combines a central body gravity model with 
multiple perturbing force models to compute the total acceleration acting on a spacecraft.

The model treats gravity as the dominant central body force and adds perturbations from 
other effects such as atmospheric drag, solar radiation pressure, third-body gravity, 
and relativistic effects.

# Type Parameters
- `N::Int`: Number of perturbing force models
- `GT <: AbstractGravityAstroModel`: Type of the gravity model
- `PT <: Tuple`: Type of the tuple containing perturbing models

# Fields
- `gravity_model::GT`: The central body gravitational force model (e.g., KeplerianGravityAstroModel, GravityHarmonicsAstroModel)
- `perturbing_models::PT`: Tuple of perturbing force models (e.g., drag, SRP, third-body, relativistic effects)

# Constructors

```julia
# With explicit gravity model and perturbations
CentralBodyDynamicsModel(gravity_model, (drag_model, srp_model, ...))

# With only perturbations (uses Keplerian gravity by default)
CentralBodyDynamicsModel((drag_model, srp_model, ...))

# With only gravity model (no perturbations)
CentralBodyDynamicsModel(gravity_model)
```

# Example

```julia
using AstroForceModels

# Create force models
gravity = GravityHarmonicsAstroModel(...)
drag = DragAstroModel(...)
srp = SRPAstroModel(...)
third_body = ThirdBodyAstroModel(...)

# Combine into dynamics model
dynamics = CentralBodyDynamicsModel(gravity, (drag, srp, third_body))

# Use in ODE integration
function satellite_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics)
end
```

# See Also
- [`build_dynamics_model`](@ref): Function to compute total acceleration
- [`AbstractDynamicsModel`](@ref): Parent abstract type
- [`acceleration`](@ref): Individual force model acceleration interface
"""
struct CentralBodyDynamicsModel{N,GT<:AbstractGravityAstroModel,PT<:Tuple} <:
       AbstractDynamicsModel where {N<:Int}
    gravity_model::GT
    perturbing_models::PT
end

"""
    CentralBodyDynamicsModel(gravity_model, perturbing_models)

Create a dynamics model with explicit gravity and perturbing force models.

# Arguments
- `gravity_model::AbstractGravityAstroModel`: Central body gravity model
- `perturbing_models::NTuple{N,AbstractAstroForceModel}`: Tuple of perturbing force models
"""
function CentralBodyDynamicsModel(
    gravity_model::AbstractGravityAstroModel, models::NTuple{N,AbstractAstroForceModel}
) where {N}
    return CentralBodyDynamicsModel{N,typeof(gravity_model),typeof(models)}(
        gravity_model, models
    )
end

"""
    CentralBodyDynamicsModel(perturbing_models)

Create a dynamics model with perturbing forces and default Keplerian gravity.

# Arguments
- `perturbing_models::NTuple{N,AbstractAstroForceModel}`: Tuple of perturbing force models
"""
function CentralBodyDynamicsModel(models::NTuple{N,AbstractAstroForceModel}) where {N}
    gravity_model = KeplerianGravityAstroModel()
    return CentralBodyDynamicsModel(gravity_model, models)
end

"""
    CentralBodyDynamicsModel(gravity_model)

Create a dynamics model with only a gravity model (no perturbations).

# Arguments
- `gravity_model::AbstractGravityAstroModel`: Central body gravity model
"""
function CentralBodyDynamicsModel(gravity_model::AbstractGravityAstroModel)
    return CentralBodyDynamicsModel(gravity_model, ())
end

"""
    build_dynamics_model(u::AbstractVector, p::ComponentVector, t::Number, models::CentralBodyDynamicsModel)

Compute the total acceleration acting on a spacecraft using a [`CentralBodyDynamicsModel`](@ref).

This function efficiently combines the central body gravity acceleration with all perturbing 
accelerations to produce the total acceleration vector. It is optimized for use in ODE 
integration routines and provides better performance than manually summing individual 
force model accelerations.

# Arguments
- `u::AbstractVector`: Current spacecraft state vector [rx, ry, rz, vx, vy, vz] in km and km/s
- `p::ComponentVector`: Simulation parameters including time references (JD), spacecraft properties, etc.
- `t::Number`: Current simulation time (typically seconds since epoch)
- `models::CentralBodyDynamicsModel`: Combined dynamics model containing gravity and perturbing forces

# Returns
- `SVector{3}`: The 3-dimensional total acceleration vector [ax, ay, az] in km/s²

# Performance Notes
The function is marked `@inline` for performance and uses compile-time optimizations 
through the tuple-based storage of perturbing models. This approach provides:
- Zero-cost abstraction over manual force summation
- Type stability for all force combinations
- Efficient memory access patterns

# Example
```julia
# Create dynamics model
dynamics = CentralBodyDynamicsModel(gravity_model, (drag_model, srp_model))

# Compute acceleration at current state
state = [6678.137, 0.0, 0.0, 0.0, 7.66, 0.0]  # km, km/s
params = ComponentVector(JD = date_to_jd(2024, 1, 1, 12, 0, 0))
time = 0.0

total_accel = build_dynamics_model(state, params, time, dynamics)
```

# See Also
- [`CentralBodyDynamicsModel`](@ref): Dynamics model structure
- [`acceleration`](@ref): Individual force model interface
- [`sum_accelerations`](@ref): Internal acceleration summation function
"""
@inline function build_dynamics_model(
    u::AbstractVector, p::ComponentVector, t::Number, models::CentralBodyDynamicsModel
)
    perturbing_accel = sum_accelerations(u, p, t, models.perturbing_models)
    central_body_accel = acceleration(u, p, t, models.gravity_model)
    return SVector{3}(
        central_body_accel[1] + perturbing_accel[1],
        central_body_accel[2] + perturbing_accel[2],
        central_body_accel[3] + perturbing_accel[3],
    )
end

"""
    sum_accelerations(u::AbstractVector, p::ComponentVector, t::Number, models::Tuple)

Internal recursive function to efficiently sum accelerations from multiple force models.

This function uses compile-time recursion over the tuple of force models to compute
the total perturbation acceleration. The recursive approach allows the compiler to
unroll the loop and inline all acceleration computations, providing optimal performance.

# Arguments
- `u::AbstractVector`: Current spacecraft state vector
- `p::ComponentVector`: Simulation parameters  
- `t::Number`: Current simulation time
- `models::Tuple`: Tuple of force models to sum over

# Returns
- `SVector{3}`: Sum of all acceleration vectors from the force models

# Implementation Notes
- Uses tail recursion for compile-time unrolling
- Base case returns zero acceleration for empty tuple
- Each recursive call processes one force model and continues with the tail
"""
@inline function sum_accelerations(
    u::AbstractVector, p::ComponentVector, t::Number, models::Tuple
)
    sum_accel = sum_accelerations(u, p, t, Base.tail(models))
    current_accel = acceleration(u, p, t, first(models))
    return SVector{3}(
        sum_accel[1] + current_accel[1],
        sum_accel[2] + current_accel[2],
        sum_accel[3] + current_accel[3],
    )
end

# Base case: empty tuple returns zero acceleration
@inline sum_accelerations(u::AbstractVector, p::ComponentVector, t::Number, models::Tuple{}) = SVector{
    3
}(
    0.0, 0.0, 0.0
)
