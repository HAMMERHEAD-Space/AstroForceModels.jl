# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Relativity Force Models for Schwarzschild, Lense-Thirring, and De Sitter Effects
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#   [1] https://link.springer.com/article/10.1007/s10569-021-10014-y
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export RelativityModel,
    relativity_accel,
    lense_thirring_acceleration,
    schwarzschild_acceleration,
    de_Sitter_acceleration

"""
    RelativityModel

Relativistic perturbation model (Schwarzschild, Lense-Thirring, de Sitter effects).

# Fields
- `central_body::ThirdBodyModel`: Central body (provides μ via `body.μ`).
- `sun_body::ThirdBodyModel`: Sun body model (for de Sitter effect, uses FrameEphemeris).
- `J::Union{SVector{3,<:Number}, Nothing}`: Angular momentum per unit mass of the central body
  in the propagation frame [km²/s]. If `nothing`, Lense-Thirring is disabled.
- `c::Number`: Speed of light [km/s].
- `γ::Number`: Post-Newtonian parameter (1.0 in GR).
- `β::Number`: Post-Newtonian parameter (1.0 in GR).
- `schwarzschild_effect::Bool`: Include Schwarzschild effect.
- `lense_thirring_effect::Bool`: Include Lense-Thirring effect.
- `de_Sitter_effect::Bool`: Include de Sitter effect.
"""
Base.@kwdef struct RelativityModel{
    CBT<:ThirdBodyModel,
    SBT<:ThirdBodyModel,
    JT<:Union{SVector{3,<:Number},Nothing},
    CT<:Number,
    GT<:Number,
    BT<:Number,
    ET<:Bool,
} <: AbstractNonPotentialBasedForce
    central_body::CBT
    sun_body::SBT
    J::JT = nothing
    c::CT = SPEED_OF_LIGHT
    γ::GT = 1.0
    β::BT = 1.0
    schwarzschild_effect::ET = true
    lense_thirring_effect::ET = true
    de_Sitter_effect::ET = true
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, relativity_model::RelativityModel)

Computes the relativistic acceleration acting on a spacecraft.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `relativity_model::RelativityModel`: Relativity model.

# Returns
- `SVector{3}`: Relativistic acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector{UT}, p::FrameAwareParams, t::TT, relativity_model::RelativityModel
) where {UT,TT}
    RT = promote_type(UT, TT)
    z = zero(RT)
    accel = SVector{3}(z, z, z)

    t_ft = ft_time(p, t)
    μ_body = relativity_model.central_body.body.μ
    c = relativity_model.c
    γ = relativity_model.γ

    if relativity_model.schwarzschild_effect
        accel =
            accel + schwarzschild_acceleration(u, μ_body; c=c, γ=γ, β=relativity_model.β)
    end

    if relativity_model.lense_thirring_effect && relativity_model.J !== nothing
        accel = accel + lense_thirring_acceleration(u, μ_body, relativity_model.J; c=c, γ=γ)
    end

    if relativity_model.de_Sitter_effect
        sun_pos, sun_vel = get_velocity(
            relativity_model.sun_body.ephem_type,
            relativity_model.sun_body.body,
            p.frames,
            t_ft,
            relativity_model.sun_body.compiled_vector6,
        )
        μ_Sun = relativity_model.sun_body.body.μ
        accel = accel + de_Sitter_acceleration(u, sun_pos, sun_vel, μ_Sun; c=c, γ=γ)
    end

    return accel
end

"""
    relativity_accel(u, r_sun, v_sun, μ_body, μ_Sun, J; c, γ, β, schwarzschild_effect, lense_thirring_effect, de_Sitter_effect)

Computes the combined relativistic acceleration (low-level function).

# Arguments
- `u::AbstractVector`: Current state.
- `r_sun::AbstractVector`: Sun position [km].
- `v_sun::AbstractVector`: Sun velocity [km/s].
- `μ_body::Number`: Central body gravitational parameter [km³/s²].
- `μ_Sun::Number`: Sun gravitational parameter [km³/s²].
- `J::AbstractVector`: Angular momentum per unit mass of central body [km²/s].

# Returns
- `SVector{3}`: Relativistic acceleration [km/s²].
"""
function relativity_accel(
    u::AbstractVector{UT},
    r_sun::AbstractVector{RT},
    v_sun::AbstractVector{VT},
    μ_body::MT,
    μ_Sun::MT2,
    J::AbstractVector{JT};
    c::CT=SPEED_OF_LIGHT,
    γ::GT=1.0,
    β::BT=1.0,
    schwarzschild_effect::Bool=true,
    lense_thirring_effect::Bool=true,
    de_Sitter_effect::Bool=true,
) where {UT,RT,VT,MT,MT2,JT,CT,GT,BT}
    AT = promote_type(UT, RT, VT, MT, MT2, JT, CT, GT, BT)
    z = zero(AT)
    accel = SVector{3}(z, z, z)

    if schwarzschild_effect
        accel = accel + schwarzschild_acceleration(u, μ_body; c=c, γ=γ, β=β)
    end
    if lense_thirring_effect
        accel = accel + lense_thirring_acceleration(u, μ_body, J; c=c, γ=γ)
    end
    if de_Sitter_effect
        accel = accel + de_Sitter_acceleration(u, r_sun, v_sun, μ_Sun; c=c, γ=γ)
    end

    return accel
end

"""
    schwarzschild_acceleration(u, μ_body; c, γ, β)

Computes the Schwarzschild relativistic acceleration.

# Returns
- `SVector{3}`: Schwarzschild acceleration [km/s²].
"""
@inline function schwarzschild_acceleration(
    u::AbstractVector{UT}, μ_body::MT; c::CT=SPEED_OF_LIGHT, γ::GT=1.0, β::BT=1.0
) where {UT,MT,CT,GT,BT}
    RT = promote_type(UT, MT, CT, GT, BT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    r_norm = norm(r)
    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    schwarzschild_pos_force = μ_body / ((c^2.0) * (r_norm^3.0))
    schwarzschild_dir =
        ((2.0 * (β + γ)) * (μ_body / r_norm) - γ * dot(ṙ, ṙ)) * r +
        2.0 * (1.0 + γ) * dot(r, ṙ) * ṙ

    schwarzschild = SVector{3,RT}(
        schwarzschild_pos_force * schwarzschild_dir[1],
        schwarzschild_pos_force * schwarzschild_dir[2],
        schwarzschild_pos_force * schwarzschild_dir[3],
    )

    return schwarzschild
end

"""
    lense_thirring_acceleration(u, μ_body, J; c, γ)

Computes the Lense-Thirring relativistic acceleration.

# Arguments
- `u::AbstractVector`: Current state.
- `μ_body::Number`: Central body gravitational parameter.
- `J::AbstractVector`: Angular momentum per unit mass of the central body [km²/s].
- `c::Number`: Speed of light [km/s].
- `γ::Number`: Post-Newtonian parameter.

# Returns
- `SVector{3}`: Lense-Thirring acceleration [km/s²].
"""
@inline function lense_thirring_acceleration(
    u::AbstractVector{UT},
    μ_body::MT,
    J::AbstractVector{JT};
    c::CT=SPEED_OF_LIGHT,
    γ::GT=1.0,
) where {UT,MT,JT,CT,GT}
    RT = promote_type(UT, MT, JT, CT, GT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    r_norm = norm(r)
    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    lense_thirring_force = (1.0 + γ) * (μ_body / ((c^2.0) * (r_norm^3.0)))
    lense_thirring_dir = ((3.0 / r_norm^2) * cross(r, ṙ) * dot(r, J) + cross(ṙ, J))

    lense_thirring = SVector{3,RT}(
        lense_thirring_force * lense_thirring_dir[1],
        lense_thirring_force * lense_thirring_dir[2],
        lense_thirring_force * lense_thirring_dir[3],
    )

    return lense_thirring
end

"""
    de_Sitter_acceleration(u, r_sun, v_sun, μ_Sun; c, γ)

Computes the de Sitter relativistic acceleration.

# Arguments
- `u::AbstractVector`: Current state.
- `r_sun::AbstractVector`: Sun position relative to central body [km].
- `v_sun::AbstractVector`: Sun velocity relative to central body [km/s].
- `μ_Sun::Number`: Sun gravitational parameter [km³/s²].

# Returns
- `SVector{3}`: de Sitter acceleration [km/s²].
"""
@inline function de_Sitter_acceleration(
    u::AbstractVector{UT},
    r_sun::AbstractVector{ST},
    v_sun::AbstractVector{VT},
    μ_Sun::MT;
    c::CT=SPEED_OF_LIGHT,
    γ::GT=1.0,
) where {UT,ST,VT,MT,CT,GT}
    RT = promote_type(UT, ST, VT, MT, CT, GT)

    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    de_sitter_force = (1.0 + 2.0 * γ) * (-μ_Sun / ((c^2.0) * (norm(-r_sun)^3.0)))
    de_sitter_dir = cross(cross(-v_sun, -r_sun), ṙ)

    de_sitter = SVector{3,RT}(
        de_sitter_force * de_sitter_dir[1],
        de_sitter_force * de_sitter_dir[2],
        de_sitter_force * de_sitter_dir[3],
    )

    return de_sitter
end
