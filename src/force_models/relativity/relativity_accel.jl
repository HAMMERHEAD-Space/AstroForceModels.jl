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
#TODO: CAN THIS BE INCORPORATED INTO THE POTENTIAL FORCES IN DROMO
"""
Relativity Astro Model struct
Contains information to compute the acceleration of relativity acting on a spacecraft.

# Fields
- `central_body::ThirdBodyModel`: The data to compute the central body's gravitational parameter.
- `sun_body::ThirdBodyModel`: The data to compute the Sun's position, velocity, and gravitational parameter.
- `eop_data::Union{EopIau1980,EopIau2000A}`: Earth Orientation Parameter data.
- `c::Number`: The speed of light [km/s].
- `γ::Number`: Post-Newtonian Parameterization parameter. γ=1 in General Relativity.
- `β::Number`: Post-Newtonian Parameterization parameter. β=1 in General Relativity.
- `schwarzschild_effect::Bool`: Include the Schwarzschild relativity effect.
- `lense_thirring_effect::Bool`: Include the Lense Thirring relativity effect.
- `de_Sitter_effect::Bool`: Include the De Sitter relativity effect.
"""
Base.@kwdef struct RelativityModel{
    CBT<:ThirdBodyModel,
    SBT<:ThirdBodyModel,
    EoT<:Union{EopIau1980,EopIau2000A},
    CT<:Number,
    GT<:Number,
    BT<:Number,
    ET<:Bool,
} <: AbstractNonPotentialBasedForce
    central_body::CBT = ThirdBodyModel(; body=EarthBody(), eop_data=fetch_iers_eop())
    sun_body::SBT = ThirdBodyModel(; body=SunBody(), eop_data=fetch_iers_eop())
    eop_data::EoT = fetch_iers_eop()
    c::CT = SPEED_OF_LIGHT
    γ::GT = 1.0
    β::BT = 1.0
    schwarzschild_effect::ET = true
    lense_thirring_effect::ET = true
    de_Sitter_effect::ET = true
end

"""
    RelativityModel(eop_data; kwargs...)

Convenience constructor that shares a single EOP dataset across all sub-models,
avoiding redundant `fetch_iers_eop()` calls.
"""
function RelativityModel(eop_data::Union{EopIau1980,EopIau2000A}; kwargs...)
    return RelativityModel(;
        central_body=ThirdBodyModel(; body=EarthBody(), eop_data=eop_data),
        sun_body=ThirdBodyModel(; body=SunBody(), eop_data=eop_data),
        eop_data=eop_data,
        kwargs...,
    )
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, relativity_model::RelativityModel)

Computes the relativistic acceleration acting on a spacecraft given a relativity model and current 
state and parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `relativity_model::RelativityModel`: Relativity model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional relativity acceleration acting on the spacecraft.

"""
@inline function acceleration(
    u::AbstractVector{UT}, p::ComponentVector{PT}, t::TT, relativity_model::RelativityModel
) where {UT,PT,TT}
    RT = promote_type(UT, PT, TT)
    z = zero(RT)
    accel = SVector{3}(z, z, z)

    current_time = current_jd(p, t)
    μ_body = relativity_model.central_body.body.μ
    c = relativity_model.c
    γ = relativity_model.γ

    if relativity_model.schwarzschild_effect
        accel =
            accel + schwarzschild_acceleration(u, μ_body; c=c, γ=γ, β=relativity_model.β)
    end

    if relativity_model.lense_thirring_effect
        R_ITRF2J2000::SatelliteToolboxTransformations.DCM{RT} = r_ecef_to_eci(
            ITRF(), J2000(), current_time, relativity_model.eop_data
        )
        J =
            SVector{3}(R_ITRF2J2000[1, 3], R_ITRF2J2000[2, 3], R_ITRF2J2000[3, 3]) *
            EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS
        accel = accel + lense_thirring_acceleration(u, μ_body, J; c=c, γ=γ)
    end

    if relativity_model.de_Sitter_effect
        sun_pos = relativity_model.sun_body(current_time, Position()) ./ 1E3
        sun_vel = relativity_model.sun_body(current_time, Velocity()) ./ 1E3
        μ_Sun = relativity_model.sun_body.body.μ
        accel = accel + de_Sitter_acceleration(u, sun_pos, sun_vel, μ_Sun; c=c, γ=γ)
    end

    return accel
end

"""
    relativity_accel(
        u::AbstractVector,
        r_sun::AbstractVector,
        v_sun::AbstractVector,
        μ_body::Number,
        μ_Sun::Number,
        J::AbstractVector;
        c::Number=SPEED_OF_LIGHT,
        γ::Number=1.0,
        β::Number=1.0,
        schwarzschild_effect::Bool=true,
        lense_thirring_effect::Bool=true,
        de_Sitter_effect::Bool=true,
    )

Computes the relativity acceleration acting on a spacecraft given a relativity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `r_sun::AbstractVector`: The position of the sun in the Earth inertial frame.
- `v_sun::AbstractVector`: The velocity of the sun in the Earth inertial frame.
- `μ_body::Number`: Gravitation Parameter of the central body.
- `μ_Sun::Number`: Gravitation Parameter of the Sun. [km^3/s^2]
- `J::AbstractVector`: Angular momentum vector per unit mass of the central body. [km^3/s^2]
- `c::Number`: Speed of Light [km/s]
- `γ::Number`: Post-Newtonian Parameterization parameter. γ=1 in General Relativity.
- `β::Number`: Post-Newtonian Parameterization parameter. β=1 in General Relativity.
- `schwarzschild_effect::Bool`: Include the Schwarzschild relativity effect.
- `lense_thirring_effect::Bool`: Include the Lense Thirring relativity effect.
- `de_Sitter_effect::Bool`: Include the De Sitter relativity effect.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional relativity acceleration acting on the spacecraft.

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
    schwarzschild_acceleration(
        u::AbstractVector, μ_body::Number; c::Number=SPEED_OF_LIGHT, γ::Number=1.0, β::Number=1.0
    )

Computes the relativity acceleration acting on a spacecraft given a relativity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `μ_body::Number`: Gravitation Parameter of the central body.
- `c::Number`: Speed of Light [km/s]
- `γ::Number`: Post-Newtonian Parameterization parameter. γ=1 in General Relativity.
- `β::Number`: Post-Newtonian Parameterization parameter. β=1 in General Relativity.

# Returns
- `schwarzschild_acceleration: SVector{3}`: The 3-dimensional schwarzschild acceleration acting on the spacecraft.

"""
@inline function schwarzschild_acceleration(
    u::AbstractVector{UT}, μ_body::MT; c::CT=SPEED_OF_LIGHT, γ::GT=1.0, β::BT=1.0
) where {UT,MT,CT,GT,BT}
    RT = promote_type(UT, MT, CT, GT, BT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    r_norm = norm(r)
    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    schwarzschild_pos_force = μ_body / ((c^2.0) * (r_norm^3.0))
    schwarzschild_dir =
        ((2.0 * (β + γ)) * (μ_body / r_norm) - γ * dot(ṙ, ṙ)) * r +
        2.0 * (1.0 + γ) * dot(r, ṙ) * ṙ

    schwarzschild = SVector{3,RT}(
        schwarzschild_pos_force * schwarzschild_dir[1],
        schwarzschild_pos_force * schwarzschild_dir[2],
        schwarzschild_pos_force * schwarzschild_dir[3],
    )

    return schwarzschild
end

"""
    lense_thirring_acceleration(
        u::AbstractVector,
        μ_body::Number,
        J::AbstractVector;
        c::Number=SPEED_OF_LIGHT,
        γ::Number=1.0,
    )

Computes the lense thirring relativity acceleration acting on a spacecraft given a relativity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `μ_body::Number`: Gravitation Parameter of the central body.
- `J::AbstractVector`: Angular momentum vector per unit mass of the central body. [km^3/s^2]
- `c::Number`: Speed of Light [km/s]
- `γ::Number`: Post-Newtonian Parameterization parameter. γ=1 in General Relativity.

# Returns
- `lense_thirring_acceleration: SVector{3}`: The 3-dimensional lense thirring acceleration acting on the spacecraft.

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
    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    lense_thirring_force = (1.0 + γ) * (μ_body / ((c^2.0) * (r_norm^3.0)))
    lense_thirring_dir = ((3.0 / r_norm^2) * cross(r, ṙ) * dot(r, J) + cross(ṙ, J))

    lense_thirring = SVector{3,RT}(
        lense_thirring_force * lense_thirring_dir[1],
        lense_thirring_force * lense_thirring_dir[2],
        lense_thirring_force * lense_thirring_dir[3],
    )

    return lense_thirring
end

"""
    de_Sitter_acceleration(
        u::AbstractVector,
        r_sun::AbstractVector,
        v_sun::AbstractVector,
        μ_Sun::Number;
        c::Number=SPEED_OF_LIGHT,
        γ::Number=1.0,
    )

Computes the relativity acceleration acting on a spacecraft given a relativity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `r_sun::AbstractVector`: The position of the sun in the Earth inertial frame.
- `v_sun::AbstractVector`: The velocity of the sun in the Earth inertial frame.
- `μ_Sun::Number`: Gravitation Parameter of the Sun. [km^3/s^2]
- `c::Number`: Speed of Light [km/s]
- `γ::Number`: Post-Newtonian Parameterization parameter. γ=1 in General Relativity.

# Returns
- `de_Sitter_acceleration: SVector{3}`: The 3-dimensional de Sitter acceleration acting on the spacecraft.

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

    ṙ = SVector{3,UT}(u[4], u[5], u[6])

    de_sitter_force = (1.0 + 2.0 * γ) * (-μ_Sun / ((c^2.0) * (norm(-r_sun)^3.0)))
    de_sitter_dir = cross(cross(-v_sun, -r_sun), ṙ)

    de_sitter = SVector{3,RT}(
        de_sitter_force * de_sitter_dir[1],
        de_sitter_force * de_sitter_dir[2],
        de_sitter_force * de_sitter_dir[3],
    )

    return de_sitter
end
