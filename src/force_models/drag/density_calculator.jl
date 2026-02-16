# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function set to compute atmospheric density from atmospheric models provided by 
#   the SatelliteToolbox ecosystem
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AtmosphericModelType, JB2008, JR1971, MSIS2000, ExpAtmo, NoAtmosphere
abstract type AtmosphericModelType end
struct JB2008 <: AtmosphericModelType end
struct JR1971 <: AtmosphericModelType end
struct MSIS2000 <: AtmosphericModelType end
struct ExpAtmo <: AtmosphericModelType end
struct NoAtmosphere <: AtmosphericModelType end

export compute_density
"""
    compute_density(JD::Number, u::AbstractVector, eop_data::EopIau1980, AtmosphereType::AtmosphericModelType)

Computes the atmospheric density at a point given the date, position, eop_data, and atmosphere type

# Arguments
- `JD::Number`: The current time of the simulation in Julian days.
- `u::AbstractVector`: The current state of the simulation.
- `eop_data::EopIau1980`: The earth orientation parameters.
- `AtmosphereType::AtmosphericModelType`: The type of atmospheric model used to compute the density. Available 
    options are Jacchia-Bowman 2008 (JB2008), Jacchia-Roberts 1971 (JR1971), NRL MSIS 2000 (MSIS2000),
    Exponential (ExpAtmo), and NoAtmosphere (NoAtmosphere)

# Returns
- `rho::Number`: The density of the atmosphere at the provided time and point [kg/m^3].

"""
function compute_density end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    eop_data::EopIau1980,
    AtmosphereType::JB2008;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    # Compute the geodetic position of the provided point
    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), JD, eop_data)
    ecef_pos = R_J20002ITRF * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    # Compute the JB2008 density if the point is less than 1000km altitude otherwise it's 0.0
    return (geodetic_pos[3] < 1000E3) *
           AtmosphericModels.jb2008(JD, geodetic_pos...; verbose=Val(false)).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    eop_data::EopIau1980,
    AtmosphereType::JR1971;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    # Compute the geodetic position of the provided point
    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), JD, eop_data)
    ecef_pos = R_J20002ITRF * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    # Compute the JR1971 density if the point is less than 2500km altitude otherwise it's 0.0
    return (geodetic_pos[3] < 2500E3) * AtmosphericModels.jr1971(
        JD, geodetic_pos...; verbose=Val(false), roots_container=roots_container
    ).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    eop_data::EopIau1980,
    AtmosphereType::MSIS2000;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    # Compute the geodetic position of the provided point
    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), JD, eop_data)
    ecef_pos = R_J20002ITRF * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    # Compute the NRLMSISE-00 density if the point is less than 1000km altitude otherwise it's 0.0
    # NRLMSISE-00 is valid from 0 to 1000 km altitude
    return (geodetic_pos[3] < 1000E3) * AtmosphericModels.nrlmsise00(
        JD, geodetic_pos[3], geodetic_pos[1], geodetic_pos[2]; verbose=Val(false), P=P
    ).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    eop_data::EopIau1980,
    AtmosphereType::ExpAtmo;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    # Compute the JR1971 density if the point is less than 2500km altitude otherwise it's 0.0
    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), JD, eop_data)
    ecef_pos = R_J20002ITRF * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    # Compute the ExpAtmo density
    return AtmosphericModels.exponential(geodetic_pos[3])
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    eop_data::EopIau1980,
    AtmosphereType::NoAtmosphere;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    # If NoAtmosphere provided return a type-stable zero
    return zero(eltype(u))
end
