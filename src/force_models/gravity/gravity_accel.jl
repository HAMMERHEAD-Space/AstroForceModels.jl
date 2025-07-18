# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Zonal Harmonics
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
include("./utils.jl")

export GravityHarmonicsAstroModel, KeplerianGravityAstroModel

abstract type AbstractGravityAstroModel <: AbstractPotentialBasedForce end

"""
Gravitational Harmonics Astro Model struct
Contains information to compute the acceleration of a Gravitational Harmonics Model acting on a spacecraft.

# Fields
- `gravity_model::AbstractGravityModel`: The gravitational potential model and coefficient data.
- `eop_data::Union{EopIau1980,EopIau2000A}`: The data compute the Earth's orientation.
- `order::Int`: The maximum order to compute the graviational potential to, a value of -1 compute the maximum order of the supplied model. (Default=-1)
- `degree::Int`: The maximum degree to compute the graviational potential to, a value of -1 compute the maximum degree of the supplied model. (Default=-1)
"""
@with_kw struct GravityHarmonicsAstroModel{GT,EoT,V,PT,DPT} <:
                AbstractGravityAstroModel where {
    GT<:AbstractGravityModel{<:Number,NT} where {NT},
    EoT<:Union{EopIau1980,EopIau2000A},
    V<:Int,
    PT<:Union{AbstractVector,Nothing},
    DPT<:Union{AbstractVector,Nothing},
}
    gravity_model::GT
    eop_data::EoT
    order::V = -1
    degree::V = -1

    #TODO: HANDLE THIS AUTONOMOUSLY, TYPING TRICKY WITH AUTODIFF
    P::PT = nothing
    dP::DPT = nothing
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the gravitational acceleration acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::GravityHarmonicsAstroModel`: Gravity model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional gravity acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel
)
    # Compute the J2000 to ITRF rotation matrix
    R_J2002ITRF = r_eci_to_ecef(J2000(), ITRF(), p.JD + t / 86400.0, grav_model.eop_data)

    # Compute the ITRF position
    itrf_pos = R_J2002ITRF * SVector{3}(u[1], u[2], u[3]) .* 1E3

    # Compute the ITRF acceleration of the spacecraft
    time = (p.JD - JD_J2000) * 86400 + t
    accel_itrf =
        GravityModels.gravitational_acceleration(
            grav_model.gravity_model,
            itrf_pos,
            time;
            max_degree=grav_model.degree,
            max_order=grav_model.order,
            P=grav_model.P,
            dP=grav_model.dP,
        ) ./ 1E3

    # Rotate into the J200 frame
    return R_J2002ITRF' * accel_itrf
end

"""
    potential(u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the gravitational potential acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::GravityHarmonicsAstroModel`: Gravity model struct containing the relevant information to compute the potential.

# Returns
- `potential: Number`: The gravitational potential acting on the spacecraft.

"""
function potential(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel
)

    # Compute the J2000 to ITRF rotation matrix
    R_J2002ITRF = r_eci_to_ecef(J2000(), ITRF(), p.JD + t / 86400.0, grav_model.eop_data)

    # Compute the ITRF position
    itrf_pos = R_J2002ITRF * SVector{3}(u[1], u[2], u[3]) .* 1E3

    time = (p.JD - JD_J2000) * 86400 + t

    U = GravityModels.gravitational_potential(
        grav_model.gravity_model,
        itrf_pos,
        time;
        max_degree=grav_model.degree,
        max_order=grav_model.order,
        P=grav_model.P,
    ) / 1E6

    return U
end

"""
    potential_time_derivative(u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel)

Computes the time derivative of the gravitational potential acting on a spacecraft given a gravity model and current state and 
parameters of an object. Based on the IAU 2006 precession-nutation model, implemention from [1].

[1] Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::GravityHarmonicsAstroModel`: Gravity model struct containing the relevant information to compute the potential time derivative.

# Returns
- `potential_time_derivative: Number`: The time derivative of the gravitational potential acting on the spacecraft.

"""
function potential_time_derivative(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel
)
    curr_jd = p.JD + t / 86400.0
    jd_tt = curr_jd + (get_Δat(curr_jd) + 32.184) / 86400.0

    # Compute the J2000 to ITRF rotation matrix
    R_J2002ITRF = r_eci_to_ecef(J2000(), ITRF(), curr_jd, grav_model.eop_data)

    # Compute the ITRF position
    itrf_pos = R_J2002ITRF * SVector{3}(u[1], u[2], u[3]) .* 1E3

    time = (p.JD - JD_J2000) * 86400 + t

    #TODO: Should This Be Done with AutoDiff?
    ∇Uₜ =
        err_iau2006(jd_tt) * GravityModels.gravitational_field_derivative(
            grav_model.gravity_model,
            itrf_pos,
            time;
            max_degree=grav_model.degree,
            max_order=grav_model.order,
            P=grav_model.P,
            dP=grav_model.dP,
        )[3] / 1E6

    return ∇Uₜ
end

"""
Gravitational Keplerian Astro Model struct
Contains information to compute the acceleration of a Gravitational Harmonics Model acting on a spacecraft.

# Fields
- `μ::Number`: The gravitational potential constant of the central body.
"""
@with_kw struct KeplerianGravityAstroModel{MT} <:
                AbstractGravityAstroModel where {MT<:Number}
    μ::MT = μ_EARTH
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the gravitational acceleration acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::KeplerianGravityAstroModel`: Gravity model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional gravity acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel
)
    r = SVector{3}(u[1], u[2], u[3])
    r_norm = norm(r)

    grav_force = -grav_model.μ / (r_norm^3)

    return SVector{3}(grav_force * r[1], grav_force * r[2], grav_force * r[3])
end

"""
    potential(u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the gravitational potential acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::KeplerianGravityAstroModel`: Gravity model struct containing the relevant information to compute the potential.

# Returns
- `potential: Number`: The gravitational potential acting on the spacecraft.

"""
function potential(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel
)
    U = -grav_model.μ / norm(SVector{3}(u[1], u[2], u[3]))

    return U
end

"""
    potential_time_derivative(u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the time derivative of the gravitational potential acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::KeplerianGravityAstroModel`: Gravity model struct containing the relevant information to compute the potential time derivative.

# Returns
- `potential_time_derivative: Number`: The time derivative of the gravitational potential acting on the spacecraft.

"""
function potential_time_derivative(
    u::AbstractVector, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel
)
    return 0.0
end

