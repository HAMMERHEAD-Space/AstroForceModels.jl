export build_dynamics_model

"""
    build_dynamics_model(u::AbstractVector, p::ComponentVector, t::Number, models::AbstractVector{AbstractAstroForceModel})

Convenience funcction to compute the overall acceleration acting on the spacecraft.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `models::AbstractVector{AbstractAstroForceModel}`: Array of acceleration models acting on the spacecraft

# Returns
- `acceleration: SVector{3}`: The 3-dimensional drag acceleration acting on the spacecraft.

"""
@generated function build_dynamics_model(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AbstractAstroForceModel},
) where {N}
    exprs = [:(acceleration(u, p, t, models[$i])) for i in 1:N]
    return :(SVector{3}($(foldl((a, b) -> :($a[1] + $b[1], $a[2] + $b[2], $a[3] + $b[3]), exprs))))
end
