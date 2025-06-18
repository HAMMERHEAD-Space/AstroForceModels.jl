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
Central Body Dynamics Model Struct
Contains information to compute the full acceleration acting on a spacecraft in a central body frame.

# Fields
- `gravity_model::AbstractGravityModel`: The gravitational AstroForceModel.
- `peturbing_models::NTuple{N,AbstractAstroForceModel}`: The perturbing AstroForceModels.
"""
struct CentralBodyDynamicsModel{N,GT<:AbstractGravityAstroModel,PT<:Tuple} <:
       AbstractDynamicsModel where {N<:Int}
    gravity_model::GT
    perturbing_models::PT
end

function CentralBodyDynamicsModel(
    gravity_model::AbstractGravityAstroModel, models::NTuple{N,AbstractAstroForceModel}
) where {N}
    return CentralBodyDynamicsModel{N,typeof(gravity_model),typeof(models)}(
        gravity_model, models
    )
end

function CentralBodyDynamicsModel(models::NTuple{N,AbstractAstroForceModel}) where {N}
    gravity_model = KeplerianGravityAstroModel()
    return CentralBodyDynamicsModel(gravity_model, models)
end

function CentralBodyDynamicsModel(gravity_model::AbstractGravityAstroModel)
    return CentralBodyDynamicsModel(gravity_model, ())
end

"""
    build_dynamics_model(u::AbstractVector, p::ComponentVector, t::Number, models::AbstractVector{AbstractAstroForceModel})

Convenience function to compute the overall acceleration acting on the spacecraft.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `models::AbstractVector{AbstractAstroForceModel}`: Array of acceleration models acting on the spacecraft

# Returns
- `acceleration: SVector{3}`: The 3-dimensional drag acceleration acting on the spacecraft.

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

@inline sum_accelerations(u::AbstractVector, p::ComponentVector, t::Number, models::Tuple{}) = SVector{
    3
}(
    0.0, 0.0, 0.0
)
