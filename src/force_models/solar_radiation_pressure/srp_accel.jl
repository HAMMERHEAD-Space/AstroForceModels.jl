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
@with_kw struct SRPAstroModel{ST,SDT,EoT,SMT} <: AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EoT<:Union{EopIau1980,EopIau2000A},
    SMT<:ShadowModelType,
}
    satellite_srp_model::ST
    sun_data::SDT
    eop_data::EoT
    shadow_model::SMT = Conical()
end

"""
    acceleration(u::AbstractArray, p::ComponentVector, t::Number, srp_model::SRPAstroModel)

Computes the srp acceleration acting on a spacecraft given a srp model and current state and 
parameters of an object.

# Arguments
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `srp_model::SRPAstroModel`: SRP model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional srp acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractArray, p::ComponentVector, t::Number, srp_model::SRPAstroModel
)
    # Compute the Sun's Position
    sun_pos = srp_model.sun_data(p.JD + t / 86400.0, Position())

    # Compute the reflectivity ballistic coefficient
    RC = reflectivity_ballistic_coefficient(u, p, t, srp_model.satellite_srp_model)

    # Return the 3-Dimensional SRP Force
    return srp_accel(u, sun_pos, RC; ShadowModel=srp_model.shadow_model)
end

"""
    srp_accel(u::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, Ψ::Number, RC::Number, t::Number; ShadowModel::ShadowModelType)

Compute the Acceleration from Solar Radiaiton Pressure

Radiation from the Sun reflects off the satellite's surface and transfers momentum perturbing the satellite's trajectory. This
force can be computed using the a Cannonball model with the following equation

                𝐚 = F * RC * Ψ * (AU/(R_sc_Sun))^2 * R̂_sc_Sun


!!! note
    Currently only Cannonball SRP is supported, to use a higher fidelity drag either use a state varying function or compute
    the ballistic coefficient further upstream

# Arguments

- `u::AbstractArray`: The current state of the spacecraft in the central body's inertial frame.
- `sun_pos::AbstractArray`: The current position of the Sun.
- `R_Sun::Number`: The radius of the Sun.
- `R_Earth::Number`: The radius of the Earth.
- `Ψ::Number`: Solar Constant at 1 Astronomical Unit.
- `RC::Number`: The solar ballistic coefficient of the satellite -- (Area/mass) * Reflectivity Coefficient [kg/m^2].
- `t::Number`: The current time of the Simulation

# Optional Arguments

- `ShadowModel::ShadowModelType`: SRP Earth Shadow Model to use. Current Options -- :Conical, :Conical_Simplified, :Cylinderical

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
function srp_accel(
    u::AbstractArray{UT},
    sun_pos::AbstractArray,
    RC::Number;
    ShadowModel::ShadowModelType=Conical(),
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
    Ψ::Number=SOLAR_FLUX,
) where {UT}
    sat_pos = SVector{3,UT}(u[1], u[2], u[3])

    # Compute the lighting factor
    F = shadow_model(sat_pos, sun_pos, ShadowModel; R_Sun=R_Sun, R_Earth=R_Earth)

    # Compute the Vector Between the Satellite and Sun
    R_spacecraft_Sun = sat_pos - sun_pos

    #Compute the SRP Force
    return SVector{3}(
        (
            F * RC * Ψ * (ASTRONOMICAL_UNIT / norm(R_spacecraft_Sun))^2 * R_spacecraft_Sun /
            norm(R_spacecraft_Sun)
        ) / 1E3,
    )
end
