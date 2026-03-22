# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function set to compute atmospheric density from atmospheric models provided by
#   the SatelliteToolbox ecosystem
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AtmosphericModelType,
    JB2008, JR1971, MSIS2000, ExpAtmo, HarrisPriester, HarrisPriesterModified, NoAtmosphere
"""
    AtmosphericModelType

Abstract base type for atmospheric density model selectors used by [`compute_density`](@ref).

Concrete subtypes: [`JB2008`](@ref), [`JR1971`](@ref), [`MSIS2000`](@ref), [`ExpAtmo`](@ref),
[`HarrisPriester`](@ref), [`HarrisPriesterModified`](@ref), [`NoAtmosphere`](@ref).
"""
abstract type AtmosphericModelType end

"""
    JB2008 <: AtmosphericModelType

Jacchia-Bowman 2008 atmospheric density model. Valid from 100 to 1000 km altitude.
Requires `SpaceIndices.init()` for automatic solar/geomagnetic index retrieval.
"""
struct JB2008 <: AtmosphericModelType end

"""
    JR1971 <: AtmosphericModelType

Jacchia-Roberts 1971 atmospheric density model. Valid from 100 to 2500 km altitude.
Requires `SpaceIndices.init()` for automatic solar/geomagnetic index retrieval.
"""
struct JR1971 <: AtmosphericModelType end

"""
    MSIS2000 <: AtmosphericModelType

NRLMSISE-00 atmospheric density model. Valid from 0 to 1000 km altitude.
Requires `SpaceIndices.init()` for automatic solar/geomagnetic index retrieval.
"""
struct MSIS2000 <: AtmosphericModelType end

"""
    ExpAtmo <: AtmosphericModelType

Simple exponential atmospheric density model. No space weather indices required.
"""
struct ExpAtmo <: AtmosphericModelType end

"""
    HarrisPriester <: AtmosphericModelType

Harris-Priester atmospheric density model. Valid from 100 to 1000 km altitude.
Uses a static density table for mean solar activity conditions.

# Fields
- `n::Int`: Cosine exponent for diurnal bulge modeling (2 ≤ n ≤ 6). Default: 4.
"""
Base.@kwdef struct HarrisPriester <: AtmosphericModelType
    n::Int = 4
end

"""
    HarrisPriesterModified <: AtmosphericModelType

Modified Harris-Priester atmospheric density model (Hatten & Russell 2017).
Provides C³ continuity, eliminates singularities, and uses cubic dependency on
81-day centered average F10.7 solar flux. Valid from 100 to 1000 km altitude.

# Fields
- `n::Number`: Cosine exponent for diurnal bulge modeling. Default: 4.
"""
Base.@kwdef struct HarrisPriesterModified{NT<:Number} <: AtmosphericModelType
    n::NT = 4
end

"""
    NoAtmosphere <: AtmosphericModelType

Sentinel type that returns zero density. Useful for disabling drag in a dynamics model
without removing the `DragAstroModel` from the perturbation tuple.
"""
struct NoAtmosphere <: AtmosphericModelType end

export compute_density
"""
    compute_density(JD, u, R_prop2bf, AtmosphereType; kwargs...)

Computes the atmospheric density at a spacecraft position.

# Arguments
- `JD::Number`: Current Julian Date.
- `u::AbstractVector`: Spacecraft state vector [km, km/s] in the propagation frame.
- `R_prop2bf`: DCM rotating from propagation frame to body-fixed (ECEF) frame.
- `AtmosphereType::AtmosphericModelType`: Atmospheric density model selector.

# Returns
- `rho::Number`: Atmospheric density [kg/m³].
"""
function compute_density end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::JB2008;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) *
           AtmosphericModels.jb2008(JD, geodetic_pos...; verbose=Val(false)).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::JR1971;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 2500E3) * AtmosphericModels.jr1971(
        JD, geodetic_pos...; verbose=Val(false), roots_container=roots_container
    ).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::MSIS2000;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) * AtmosphericModels.nrlmsise00(
        JD, geodetic_pos[3], geodetic_pos[1], geodetic_pos[2]; verbose=Val(false), P=P
    ).total_density
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::ExpAtmo;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return AtmosphericModels.exponential(geodetic_pos[3])
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::HarrisPriester;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) *
           AtmosphericModels.harrispriester(JD, geodetic_pos...; n=AtmosphereType.n)
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::HarrisPriesterModified;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    ecef_pos = R_prop2bf * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) * AtmosphericModels.harrispriester_modified(
        JD, geodetic_pos...; n=AtmosphereType.n
    )
end

@inline function compute_density(
    JD::Number,
    u::AbstractVector,
    R_prop2bf,
    AtmosphereType::NoAtmosphere;
    roots_container::Union{Nothing,AbstractVector}=nothing,
    P::Union{Nothing,AbstractMatrix}=nothing,
)
    return zero(eltype(u))
end
