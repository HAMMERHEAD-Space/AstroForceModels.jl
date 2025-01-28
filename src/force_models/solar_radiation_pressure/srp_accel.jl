# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Solar Radiation Pressure
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] https://ai-solutions.com/_freeflyeruniversityguide/solar_radiation_pressure.htm
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export SRPAstroModel, srp_accel

"""
SRP Astro Model struct
Contains information to compute the acceleration of a SRP a spacecraft.

# Fields
- `satellite_srp_model::AbstractSatelliteDragModel`: The satellite srp model for computing the ballistic coefficient.
- `sun_data::ThirdBodyModel`: The data to compute the Sun's position.
"""
@with_kw struct SRPAstroModel{ST,SDT,EoT,SMT,RST,ROT,PT,AUT} <:
                AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EoT<:Union{EopIau1980,EopIau2000A},
    SMT<:ShadowModelType,
}
    satellite_srp_model::ST
    sun_data::SDT
    eop_data::EoT
    shadow_model::SMT = Conical()

    R_Sun::RST = R_SUN
    R_Occulting::ROT = R_EARTH
    Ψ::PT = SOLAR_FLUX
    AU::AUT = ASTRONOMICAL_UNIT
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, srp_model::SRPAstroModel)

Computes the srp acceleration acting on a spacecraft given a srp model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `srp_model::SRPAstroModel`: SRP model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional srp acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, srp_model::SRPAstroModel
)
    # Compute the Sun's Position
    sun_pos = srp_model.sun_data(p.JD + t / 86400.0, Position())

    # Compute the reflectivity ballistic coefficient
    RC = reflectivity_ballistic_coefficient(u, p, t, srp_model.satellite_srp_model)

    # Return the 3-Dimensional SRP Force
    return srp_accel(
        u,
        sun_pos,
        RC;
        ShadowModel=srp_model.shadow_model,
        R_Sun=srp_model.R_Sun,
        R_Occulting=srp_model.R_Occulting,
        Ψ=srp_model.Ψ,
        AU=srp_model.AU,
    )
end

"""
    srp_accel(u::AbstractVector, sun_pos::AbstractVector, R_Sun::Number, R_Occulting::Number, Ψ::Number, RC::Number, t::Number; ShadowModel::ShadowModelType)

Compute the Acceleration from Solar Radiaiton Pressure

Radiation from the Sun reflects off the satellite's surface and transfers momentum perturbing the satellite's trajectory. This
force can be computed using the a Cannonball model with the following equation

                𝐚 = F * RC * Ψ * (AU/(R_sc_Sun))^2 * R̂_sc_Sun


!!! note
    Currently only Cannonball SRP is supported, to use a higher fidelity drag either use a state varying function or compute
    the ballistic coefficient further upstream

# Arguments

- `u::AbstractVector`: The current state of the spacecraft in the central body's inertial frame.
- `sun_pos::AbstractVector`: The current position of the Sun.
- `R_Sun::Number`: The radius of the Sun.
- `R_Occulting::Number`: The radius of the Earth.
- `Ψ::Number`: Solar Constant at 1 Astronomical Unit.
- `RC::Number`: The solar ballistic coefficient of the satellite -- (Area/mass) * Reflectivity Coefficient [kg/m^2].
- `t::Number`: The current time of the Simulation

# Optional Arguments

- `ShadowModel::ShadowModelType`: SRP Earth Shadow Model to use. Current Options -- :Conical, :Conical_Simplified, :Cylinderical

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
function srp_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector,
    RC::Number;
    ShadowModel::ShadowModelType=Conical(),
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    Ψ::Number=SOLAR_FLUX,
    AU::Number=ASTRONOMICAL_UNIT,
) where {UT}
    sat_pos = SVector{3,UT}(u[1], u[2], u[3])

    # Compute the lighting factor
    F = shadow_model(sat_pos, sun_pos, ShadowModel; R_Sun=R_Sun, R_Occulting=R_Occulting)

    # Compute the Vector Between the Satellite and Sun
    R_spacecraft_Sun = sat_pos - sun_pos
    R_sc_sun_norm = norm(R_spacecraft_Sun)

    F_srp = F * RC * Ψ * (AU / R_sc_sun_norm)^2 / 1E3

    #Compute the SRP Force
    return SVector{3}(
        F_srp * R_spacecraft_Sun[1] / R_sc_sun_norm,
        F_srp * R_spacecraft_Sun[2] / R_sc_sun_norm,
        F_srp * R_spacecraft_Sun[3] / R_sc_sun_norm,
    )
end
