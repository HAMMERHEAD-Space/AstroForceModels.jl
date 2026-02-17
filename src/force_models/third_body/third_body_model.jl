# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Third Body Model and Ephemeris Functions
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractEphemerisType, Vallado
abstract type AbstractEphemerisType end
struct Vallado <: AbstractEphemerisType end

export Position, Velocity
abstract type EphemerisReturn end
struct Position <: EphemerisReturn end
struct Velocity <: EphemerisReturn end

export ThirdBodyModel
"""
Third Body Model Astro Model struct
Contains information to compute the acceleration of a third body force acting on a spacecraft.

# Fields
- `body::CelestialBody`: Celestial body acting on the craft.
- `ephem_type::AbstractEphemerisType`: Ephemeris type used to compute body's position. Options are currently Vallado().
"""
Base.@kwdef struct ThirdBodyModel{
    BT<:CelestialBody,EoT<:Union{EopIau1980,EopIau2000A,Nothing},EpT<:AbstractEphemerisType
} <: AbstractNonPotentialBasedForce
    body::BT = SunBody()
    eop_data::EoT = nothing
    ephem_type::EpT = Vallado()
end

#TODO: EXPAND TO SPICE WITH EXTENSIONS
"""
Computes the position of the celestial body using Vallado's ephemeris

# Arguments
- `ephem_type::Vallado`: Ephemeris type used to compute body's position.
- `body::CelestialBody`: Celestial body acting on the craft.
- `time::Number`: Current time of the simulation in Julian days.

# Returns
- `body_position::SVector{3}`: The 3-dimensional third body position in the J2000 frame [m].
"""
function get_position(
    ephem_type::Vallado, body::CelestialBody, eop_data::T, time::TT
) where {T<:Union{EopIau1980,EopIau2000A,Nothing},TT}

    # Compute the MOD frame in the J2000 frame to rotate the body's position vector
    R_MOD2J2000::SatelliteToolboxTransformations.DCM{TT} = r_eci_to_eci(
        MOD(), J2000(), time, eop_data
    )

    pos_mod = _position_mod(body, time)

    return R_MOD2J2000 * pos_mod
end

_position_mod(::CelestialBody{:Sun}, time) = sun_position_mod(time)
_position_mod(::CelestialBody{:Moon}, time) = moon_position_mod(time)
function _position_mod(body::CelestialBody{Name}, time) where {Name}
    throw(ArgumentError("Vallado position ephemeris is not supported for $Name"))
end

"""
Computes the velocity of the celestial body using Vallado's ephemeris

# Arguments
- `ephem_type::Vallado`: Ephemeris type used to compute body's velocity.
- `body::CelestialBody`: Celestial body acting on the craft.
- `time::Number`: Current time of the simulation in Julian days.

# Returns
- `body_velocity::SVector{3}`: The 3-dimensional third body velocity in the J2000 frame.
"""
function get_velocity(
    ephem_type::Vallado, body::CelestialBody, eop_data::T, time::TT
) where {T<:Union{EopIau1980,EopIau2000A,Nothing},TT}
    vel_mod = _velocity_mod(body, time)

    # Compute the MOD frame in the J2000 frame to rotate the body's velocity vector
    R_MOD2J2000::SatelliteToolboxTransformations.DCM{TT} = r_eci_to_eci(
        MOD(), J2000(), time, eop_data
    )

    return R_MOD2J2000 * vel_mod
end

_velocity_mod(::CelestialBody{:Sun}, time) = sun_velocity_mod(time)
function _velocity_mod(body::CelestialBody{Name}, time) where {Name}
    throw(ArgumentError("Vallado velocity ephemeris is not supported for $Name"))
end

#TODO: ADD FULL STATE WITH SPICE SUPPORT
"""
Convenience to compute the ephemeris position of a CelestialBody in a ThirdBodyModel
Wraps get_position().

# Arguments
- `time::Number`: Current time of the simulation in seconds.

# Returns
- `body_position: SVector{3}`: The 3-dimensional third body position in the J2000 frame.

"""
function (model::ThirdBodyModel)(t::Number, return_type::Position)
    return get_position(model.ephem_type, model.body, model.eop_data, t)
end

function (model::ThirdBodyModel)(t::Number, return_type::Velocity)
    return get_velocity(model.ephem_type, model.body, model.eop_data, t)
end
