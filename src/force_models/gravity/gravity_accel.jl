# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Gravitational Acceleration Models (Keplerian and Spherical Harmonics)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Montenbruck, O., & Gill, E. (2000). Satellite Orbits: Models, Methods, and Applications. Springer.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractGravityAstroModel, GravityHarmonicsAstroModel, KeplerianGravityAstroModel

abstract type AbstractGravityAstroModel <: AbstractPotentialBasedForce end

"""
    GravityHarmonicsAstroModel

Gravitational harmonics force model using frame-system rotations.

# Constructor
    GravityHarmonicsAstroModel(; gravity_model, body_fixed_frame, propagation_frame=:ICRF,
        order=-1, degree=-1, P=nothing, dP=nothing, frames=nothing)

# Arguments
- `gravity_model::AbstractGravityModel`: The gravitational potential model and coefficient data.
- `body_fixed_frame::Symbol`: The body-fixed frame (e.g., `:ITRF`, `:ErosBodyFixed`).
- `propagation_frame::Symbol`: The frame in which the state vector is expressed (default: `:ICRF`).
- `order::Int`: Maximum order of the spherical harmonic expansion (-1 for full model).
- `degree::Int`: Maximum degree of the spherical harmonic expansion (-1 for full model).
- `P`: Pre-allocated Legendre polynomial buffer (e.g., `MMatrix{N,N,Float64}`).
- `dP`: Pre-allocated Legendre polynomial derivative buffer.
- `frames`: Optional `FrameSystem`. When provided, pre-compiles the rotation between
  `propagation_frame` and `body_fixed_frame` for allocation-free evaluation. Falls back to
  runtime `rotation3`/`rotation6` when `nothing` or when the frame pair is not directly connected.

# Example
```julia
grav_model = GravityHarmonicsAstroModel(;
    gravity_model = grav_coeffs,
    body_fixed_frame = :ITRF,
    propagation_frame = :ICRF,
    order = 36, degree = 36,
    P  = MMatrix{37,37,Float64}(zeros(37, 37)),
    dP = MMatrix{37,37,Float64}(zeros(37, 37)),
    frames = my_frame_system,   # optional: enables allocation-free rotation
)
```
"""
struct GravityHarmonicsAstroModel{
    GT<:AbstractGravityModel,
    V<:Int,
    PT<:Union{AbstractArray,Nothing},
    DPT<:Union{AbstractArray,Nothing},
    CR3,
    CR6,
} <: AbstractGravityAstroModel
    gravity_model::GT
    body_fixed_frame::Symbol
    propagation_frame::Symbol
    order::V
    degree::V
    P::PT
    dP::DPT
    compiled_rotation3::CR3
    compiled_rotation6::CR6
end

function GravityHarmonicsAstroModel(;
    gravity_model::GT,
    body_fixed_frame::Symbol,
    propagation_frame::Symbol=:ICRF,
    order::V=-1,
    degree::V=-1,
    P::PT=nothing,
    dP::DPT=nothing,
    frames=nothing,
) where {GT<:AbstractGravityModel,V<:Int,PT,DPT}
    cr3 = nothing
    cr6 = nothing
    if !isnothing(frames) && propagation_frame != body_fixed_frame
        try
            cr3 = compile_rotation3(frames, propagation_frame, body_fixed_frame)
            cr6 = compile_rotation6(frames, propagation_frame, body_fixed_frame)
        catch
            cr3 = nothing
            cr6 = nothing
        end
    end
    return GravityHarmonicsAstroModel(
        gravity_model, body_fixed_frame, propagation_frame, order, degree, P, dP, cr3, cr6
    )
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the gravitational harmonics acceleration using frame-system rotations.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::GravityHarmonicsAstroModel`: Gravity harmonics model.

# Returns
- `SVector{3}`: Gravitational acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector,
    p::FrameAwareParams,
    t::Number,
    grav_model::GravityHarmonicsAstroModel,
)
    t_ft = ft_time(p, t)

    # Get the rotation from propagation frame to body-fixed frame
    R_rot = if isnothing(grav_model.compiled_rotation3)
        rotation3(p.frames, grav_model.propagation_frame, grav_model.body_fixed_frame, t_ft)
    else
        grav_model.compiled_rotation3(t_ft)
    end
    R_prop2bf = R_rot.m[1]  # Extract DCM from Rotation

    # Compute the body-fixed position (convert km to m for GravityModels)
    bf_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3]) .* 1E3

    # Compute the body-fixed acceleration (returns m/s²)
    time = ft_time(p, t)
    accel_bf =
        GravityModels.gravitational_acceleration(
            grav_model.gravity_model,
            bf_pos,
            time;
            max_degree=grav_model.degree,
            max_order=grav_model.order,
            P=grav_model.P,
            dP=grav_model.dP,
        ) ./ 1E3  # m/s² → km/s²

    # Rotate back to propagation frame
    return R_prop2bf' * accel_bf
end

"""
    potential(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the gravitational potential using frame-system rotations.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::GravityHarmonicsAstroModel`: Gravity harmonics model.

# Returns
- `Number`: Gravitational potential [km²/s²].
"""
function potential(
    u::AbstractVector,
    p::FrameAwareParams,
    t::Number,
    grav_model::GravityHarmonicsAstroModel,
)
    t_ft = ft_time(p, t)

    R_rot = if isnothing(grav_model.compiled_rotation3)
        rotation3(p.frames, grav_model.propagation_frame, grav_model.body_fixed_frame, t_ft)
    else
        grav_model.compiled_rotation3(t_ft)
    end
    R_prop2bf = R_rot.m[1]

    bf_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3]) .* 1E3

    time = ft_time(p, t)

    U =
        -GravityModels.gravitational_potential(
            grav_model.gravity_model,
            bf_pos,
            time;
            max_degree=grav_model.degree,
            max_order=grav_model.order,
            P=grav_model.P,
        ) / 1E6  # m²/s² → km²/s²

    return U
end

"""
    potential_time_derivative(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the time derivative of the gravitational potential using frame-system rotations.
Uses rotation6 to get the rotation rate from the frame system instead of the IAU 2006 model.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::GravityHarmonicsAstroModel`: Gravity harmonics model.

# Returns
- `Number`: Time derivative of gravitational potential [km²/s³].
"""
function potential_time_derivative(
    u::AbstractVector,
    p::FrameAwareParams,
    t::Number,
    grav_model::GravityHarmonicsAstroModel,
)
    t_ft = ft_time(p, t)

    # Get DCM and its time derivative via rotation6
    R6 = if isnothing(grav_model.compiled_rotation6)
        rotation6(p.frames, grav_model.propagation_frame, grav_model.body_fixed_frame, t_ft)
    else
        grav_model.compiled_rotation6(t_ft)
    end

    # Extract the DCM and its time derivative from Rotation{2}
    R_prop2bf = R6.m[1]  # DCM
    dR = R6.m[2]         # dDCM/dt

    # Compute the body-fixed position (convert km to m)
    bf_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3]) .* 1E3

    time = ft_time(p, t)

    # ω_skew = dR * R' gives the skew-symmetric angular velocity matrix
    ω_skew = dR * R_prop2bf'
    # The Earth rotation rate equivalent is extracted from the z-component of the angular velocity
    # For a general body, we need the magnitude of the angular velocity
    # The skew-symmetric matrix has ω components at off-diagonal positions:
    # [0, -ωz, ωy; ωz, 0, -ωx; -ωy, ωx, 0]
    # The field derivative uses the tangential velocity component perpendicular to the body-fixed z-axis

    # Get the gravitational field derivative in body-fixed frame
    field_deriv = GravityModels.gravitational_field_derivative(
        grav_model.gravity_model,
        bf_pos,
        time;
        max_degree=grav_model.degree,
        max_order=grav_model.order,
        P=grav_model.P,
        dP=grav_model.dP,
    )

    # The potential time derivative is ω · (r × ∇U) in the body-fixed frame
    # For a body rotating about z-axis: ∂U/∂t = ω_z * ∂U/∂λ
    # More generally, use the angular velocity vector dotted with the cross product
    ω_vec = SVector{3}(-ω_skew[2, 3], ω_skew[1, 3], -ω_skew[1, 2])
    ω_mag = norm(ω_vec)

    # field_deriv[3] is the tangential component (∂U/∂λ / r) for a body rotating about z-axis
    # For the general case, we use the rotation rate magnitude times the tangential derivative
    ∇Uₜ = ω_mag * field_deriv[3] / 1E6  # m²/s³ → km²/s³

    return ∇Uₜ
end

"""
    KeplerianGravityAstroModel{MT<:Number} <: AbstractGravityAstroModel

Keplerian (point-mass) gravity model.

# Fields
- `μ::Number`: Gravitational parameter of the central body [km³/s²].
"""
Base.@kwdef struct KeplerianGravityAstroModel{MT<:Number} <: AbstractGravityAstroModel
    μ::MT
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the Keplerian gravitational acceleration.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::KeplerianGravityAstroModel`: Keplerian gravity model.

# Returns
- `SVector{3}`: Gravitational acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector,
    p::FrameAwareParams,
    t::Number,
    grav_model::KeplerianGravityAstroModel,
)
    r = SVector{3}(u[1], u[2], u[3])
    r_norm = norm(r)

    grav_force = -grav_model.μ / (r_norm^3)

    return SVector{3}(grav_force * r[1], grav_force * r[2], grav_force * r[3])
end

"""
    potential(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the Keplerian gravitational potential.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::KeplerianGravityAstroModel`: Keplerian gravity model.

# Returns
- `Number`: Gravitational potential [km²/s²].
"""
function potential(
    u::AbstractVector,
    p::FrameAwareParams,
    t::Number,
    grav_model::KeplerianGravityAstroModel,
)
    U = -grav_model.μ / norm(SVector{3}(u[1], u[2], u[3]))

    return U
end

"""
    potential_time_derivative(u::AbstractVector, p::FrameAwareParams, t::Number, grav_model::KeplerianGravityAstroModel)

Time derivative of the Keplerian gravitational potential (always zero for point mass).

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `grav_model::KeplerianGravityAstroModel`: Keplerian gravity model.

# Returns
- `Number`: Always zero for point mass (potential is time-independent).
"""
function potential_time_derivative(
    u::AbstractVector{UT},
    p::FrameAwareParams,
    t::TT,
    grav_model::KeplerianGravityAstroModel,
) where {UT<:Number,TT<:Number}
    return zero(promote_type(UT, TT))
end
