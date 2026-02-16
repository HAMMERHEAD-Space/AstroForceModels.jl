#TODO: THIS PROBABLY BELONGS IN SATELLITE TOOLKIT CELESTIAL BODY
export CelestialBody
"""
    CelestialBody{Name, T<:Number}

Celestial body representation with compile-time identity for type-stable dispatch.

# Type Parameters
- `Name::Symbol`: Compile-time body identifier (e.g., `:Sun`, `:Moon`, `:Earth`)
- `T<:Number`: Numeric type for physical parameters

# Fields
- `central_body::Symbol`: Name of the central body being orbited
- `jpl_code::Int`: NAIF ID Code
- `μ::T`: Gravitational parameter [km^3/s^2]
- `Req::T`: Equatorial radius [km]
"""
struct CelestialBody{Name,T<:Number}
    central_body::Symbol
    jpl_code::Int
    μ::T
    Req::T
end

function CelestialBody(
    name::Symbol, central_body::Symbol, jpl_code::Int, μ::T, Req::T
) where {T<:Number}
    return CelestialBody{name,T}(central_body, jpl_code, μ, Req)
end

Base.nameof(::CelestialBody{Name}) where {Name} = Name

export SunBody, MoonBody, EarthBody
function SunBody(; T::DataType=Float64)
    return CelestialBody{:Sun,T}(
        :None,                          # Central Body
        1,                              # NAIF ID Code
        T(μ_SUN),                       # μ [km^3/s^2]
        T(R_SUN),                       # Equatorial Radius [km]
    )
end

function EarthBody(; T::DataType=Float64)
    return CelestialBody{:Earth,T}(
        :Sun,                           # Central Body
        399,                            # NAIF ID Code
        T(μ_EARTH),                     # μ [km^3/s^2]
        T(R_EARTH),                     # Equatorial Radius [km]
    )
end

function MoonBody(; T::DataType=Float64)
    return CelestialBody{:Moon,T}(
        :Earth,                         # Central Body
        301,                            # NAIF ID Code
        T(μ_MOON),                      # μ [km^3/s^2]
        T(R_MOON),                      # Equatorial Radius [km]
    )
end
