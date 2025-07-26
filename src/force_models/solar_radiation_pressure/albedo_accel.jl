# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Earth Albedo Radiation Pressure
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Knocke, P. C., Ries, J. C., and Tapley, B. D. (1988). "Earth radiation pressure 
#       effects on satellites." Proceedings of AIAA/AAS Astrodynamics Conference, pp. 577-587.
#   [2] Borderies, N., & Longaretti, P. Y. (1990). "A new treatment of the albedo radiation 
#       pressure in the case of a uniform albedo and of a spherical satellite." 
#       Celestial Mechanics and Dynamical Astronomy, 49(1), 69-98.
#   [3] Rubincam, D. P., & Weiss, N. R. (1986). "Earth albedo and the orbit of Lageos." 
#       Celestial Mechanics, 38(3), 233-296.
#   [4] Vielberg, K., & Kusche, J. (2020). "Extended forward and inverse modeling of 
#       radiation pressure accelerations for LEO satellites." Journal of Geodesy, 94(4), 1-29.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AlbedoAstroModel, albedo_accel, initialize_albedo_cache!
export AbstractAlbedoModel, UniformAlbedoModel, VariableAlbedoModel

# Set up integration domain in radians
# Latitude: -π/2 to π/2, Longitude: -π to π
const ALBEDO_INTEGRATION_DOMAIN = (SVector{2}(-π/2, -π), SVector{2}(π/2, π))

"""
Abstract type for albedo radiation models used in albedo force calculations.
"""
abstract type AbstractAlbedoModel{AT<:Number, ET<:Number} end

"""
Uniform Albedo Model struct
Simple uniform albedo model with constant reflection and emission coefficients.

# Fields
- `visible_albedo::Number`: Visible light albedo coefficient (0-1)
- `infrared_emissivity::Number`: Infrared emissivity coefficient (0-1)

# Default Values
The default values (visible_albedo=0.3, infrared_emissivity=0.7) represent commonly 
accepted Earth-averaged values used in astrodynamics literature. These are based on 
Earth radiation budget studies and are referenced in works such as Knocke et al. (1988), 
Borderies & Longaretti (1990), and other albedo radiation pressure studies. The 0.3 
value represents approximately 30% of incoming solar radiation being reflected, while 
0.7 represents Earth's average thermal emission characteristics.
"""
@with_kw struct UniformAlbedoModel{AT,ET} <: AbstractAlbedoModel{AT,ET}
    visible_albedo::AT = 0.3      # Earth's average visible albedo (literature consensus)
    infrared_emissivity::ET = 0.7 # Earth's average infrared emissivity (literature consensus)
end

"""
Albedo Astro Model struct
Contains information to compute the acceleration from Earth albedo radiation pressure.

# Fields
- `satellite_shape_model::AbstractSatelliteSRPModel`: The satellite shape model for computing the ballistic coefficient.
- `sun_data::ThirdBodyModel`: The data to compute the Sun's position.
- `body_albedo_model::AbstractAlbedoModel{<:Number, <:Number}`: The Earth albedo radiation model.
- `eop_data::Union{EopIau1980,EopIau2000A}`: Earth orientation parameters.
- `solar_flux::Number`: Solar flux at 1 AU [W/m²].
- `speed_of_light::Number`: Speed of light [km/s].
- `integral_reltol::Number`: Relative tolerance for numerical integration.
- `integral_abstol::Number`: Absolute tolerance for numerical integration.
"""
@with_kw mutable struct AlbedoAstroModel{ST,SDT,EAT,EoT,SFT,CT,AUT,IRT,IAT,ALG,CAT} <:
                AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel{AT, ET} where {AT<:Number, ET<:Number},
    EoT<:Union{EopIau1980,EopIau2000A},
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    IRT<:AbstractFloat,
    IAT<:AbstractFloat,
    ALG<:SciMLBase.AbstractIntegralAlgorithm,
    CAT<:Union{Nothing, Integrals.IntegralCache},
}
    satellite_shape_model::ST
    sun_data::SDT
    body_albedo_model::EAT
    eop_data::EoT
    
    solar_flux::SFT = SOLAR_FLUX
    speed_of_light::CT = SPEED_OF_LIGHT
    AU::AUT = ASTRONOMICAL_UNIT
    integral_algorithm::ALG = HCubatureJL()
    integral_reltol::IRT = 1e-4
    integral_abstol::IAT = 1e-8
    integral_cache::Ref{CAT} = Ref{Union{Nothing, Integrals.IntegralCache}}(nothing)  # Pre-allocated integration cache in a Ref
end

struct AlbedoFunctor{
    ST<:Number,
    RT<:Number,
    SUT<:Number,
    RCT<:Number,
    TT<:Number,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    AM<:AbstractAlbedoModel{<:Number, <:Number},
}
    sat_pos::SVector{3,ST}
    R_ECEF2ECI::DCM{RT}
    sun_pos::SVector{3,SUT}
    RC::RCT
    body_albedo_model::AM
    current_time::TT
    solar_flux::SFT
    speed_of_light::CT
    AU::AUT
end

# call‐overload invokes your existing albedo_integrand
function (f::AlbedoFunctor{ST,RT,SUT,RCT,TT,SFT,CT,AUT,AM})(x, p) where {
    ST<:Number,
    RT<:Number,
    SUT<:Number,
    RCT<:Number,
    TT<:Number,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    AM<:AbstractAlbedoModel{AT,ET} where {AT<:Number, ET<:Number},
}
    return albedo_integrand(x,
                            f.sat_pos,
                            f.R_ECEF2ECI,
                            f.sun_pos,
                            f.RC,
                            f.body_albedo_model,
                            f.current_time,
                            f.solar_flux,
                            f.speed_of_light,
                            f.AU)
end

SciMLBase.isinplace(::AlbedoFunctor, args...; kwargs...) = false
SciMLBase.numargs(::AlbedoFunctor) = 2

"""
Initialize the integration cache for the albedo model to avoid allocations during runtime.
This should be called once after creating the AlbedoAstroModel and before using it for integration.
"""
function initialize_albedo_cache!(albedo_model::AlbedoAstroModel)
    # Create a dummy functor to initialize the cache
    dummy_sat_pos = SVector{3,Float64}(7000.0, 0.0, 0.0)  # km
    dummy_R_ECEF2ECI = DCM(1.0I)  # Identity matrix
    dummy_sun_pos = SVector{3,Float64}(1.496e8, 0.0, 0.0)  # km (1 AU)
    dummy_RC = 0.1  # m²/kg
    dummy_time = 0.0  # Julian days
    
    dummy_functor = AlbedoFunctor(
        dummy_sat_pos, 
        dummy_R_ECEF2ECI, 
        dummy_sun_pos, 
        dummy_RC, 
        albedo_model.body_albedo_model, 
        dummy_time, 
        albedo_model.solar_flux, 
        albedo_model.speed_of_light, 
        albedo_model.AU
    )
    
    # Create integration problem
    prob = IntegralProblem(dummy_functor, ALBEDO_INTEGRATION_DOMAIN)
    
    # Initialize cache using init function
    cache = init(prob, albedo_model.integral_algorithm)
    
    # Store the cache for reuse
    albedo_model.integral_cache[] = cache
    
    return nothing
end

function albedo_integrand(
    x::AbstractVector{XT}, 
    sat_pos::AbstractVector{ST}, 
    R_ECEF2ECI::DCM{RT}, 
    sun_pos::AbstractVector{SUT}, 
    RC::RCT, 
    body_albedo_model::AM, 
    current_time::TT, 
    solar_flux::SFT, 
    speed_of_light::CT, 
    AU::AUT
) where {
    XT<:Number,
    ST<:Number,
    RT<:Number,
    SUT<:Number,
    RCT<:Number,
    TT<:Number,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    AT<:Number, 
    ET<:Number,
    AM<:AbstractAlbedoModel{AT, ET},
}
    RET = promote_type(XT, ST, RT, SUT, RCT, TT, SFT, CT, AUT, AT, ET)

    lat_rad, lon_rad = x

    # Surface element position in ECEF coordinates using SatelliteToolboxTransformations
    # Note: geodetic_to_ecef expects altitude above ellipsoid, so we use h=0 for surface points
    surface_pos_ecef = geodetic_to_ecef(lat_rad, lon_rad, 0.0) ./ 1E3 # Convert meters to km

    # Transform surface position to ECI
    surface_pos = R_ECEF2ECI * surface_pos_ecef

    # Vector from surface element to satellite
    surface_to_sat = sat_pos - surface_pos
    distance = norm(surface_to_sat)

    # Angle between surface normal and satellite direction (numerically stable)
    # Note: angle_between_vectors normalizes internally, so pass raw vectors
    angle_from_nadir = angle_between_vectors(surface_pos, surface_to_sat)
    cos_alpha = cos(angle_from_nadir)

    # Skip if surface element is not visible (behind horizon)
    if cos_alpha <= 0 || angle_from_nadir > π / 2
        return SVector{3,RET}(0.0, 0.0, 0.0)
    end

    # Compute surface element area differential (Jacobian for spherical coordinates)
    # dΩ = R² * cos(lat) * d_lat * d_lon (integrating directly over radians)
    R_earth = norm(surface_pos)
    jacobian = R_earth^2 * cos(lat_rad)  # Jacobian in km²

    # Compute solar illumination of surface element
    # Vector from surface element to Sun
    surface_to_sun = sun_pos - surface_pos
    surface_to_sun_norm = norm(surface_to_sun)

    # Solar flux at Earth's distance (W/m²)
    solar_flux_at_earth = solar_flux * (AU / surface_to_sun_norm)^2

    # Angle between surface normal and sun direction (numerically stable)
    # Note: angle_between_vectors normalizes internally, so pass raw vectors
    solar_zenith_angle = angle_between_vectors(surface_pos, surface_to_sun)
    cos_solar_zenith = cos(solar_zenith_angle)

    # Compute shortwave and longwave fluxes
    total_radiance = compute_earth_radiation_fluxes(
        body_albedo_model,
        lat_rad,
        lon_rad,
        cos_solar_zenith,
        solar_flux_at_earth,
        current_time,
    )

    # Radiation pressure from this surface element
    # Follow SRP pattern: F = RC * pressure * geometric_factor / 1E3
    # For albedo: pressure = total_radiance / speed_of_light
    # geometric_factor = jacobian * cos_alpha / (π * distance^2)
    # Note: jacobian already in km², distance in km, so no unit conversion needed here
    pressure = total_radiance / speed_of_light  # N/m² (pressure units like SOLAR_FLUX)
    geometric_factor = jacobian * cos_alpha / (π * distance^2)  # km² / km² = dimensionless

    # Following SRP formula exactly: F = RC * pressure * geometric_factor / 1E3
    F_albedo = RC * pressure * geometric_factor / 1E3

    # Convert to acceleration vector
    # Direction unit vector from surface to satellite  
    surface_to_sat_unit = surface_to_sat / distance
    element_acceleration = SVector{3,RET}(
        F_albedo * surface_to_sat_unit[1],
        F_albedo * surface_to_sat_unit[2],
        F_albedo * surface_to_sat_unit[3],
    )  # Already in km/s²

    return element_acceleration
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, albedo_model::AlbedoAstroModel)

Computes the albedo acceleration acting on a spacecraft given an albedo model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `albedo_model::AlbedoAstroModel`: Albedo model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional albedo acceleration acting on the spacecraft.
"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, albedo_model::AlbedoAstroModel{ST,SDT,EAT,EoT,SFT,CT,AUT,IRT,IAT,ALG,CAT}
) where {ST,SDT,EAT,EoT,SFT,CT,AUT,IRT,IAT,ALG,CAT}
    # Compute the Sun's position
    sun_pos = albedo_model.sun_data(p.JD + t / 86400.0, Position())

    # Compute the reflectivity ballistic coefficient
    # TODO: Probably not the same as the SRP ballistic coefficient since it would be a function of the satellite's orientation
    RC = reflectivity_ballistic_coefficient(u, p, t, albedo_model.satellite_shape_model)

    # Return the 3-dimensional Albedo Force
    return albedo_accel(
        u,
        sun_pos,
        RC,
        p.JD + t / 86400.0,
        albedo_model.body_albedo_model,
        albedo_model.eop_data;
        solar_flux=albedo_model.solar_flux,
        AU=albedo_model.AU,
        speed_of_light=albedo_model.speed_of_light,
        integral_algorithm=albedo_model.integral_algorithm,
        reltol=albedo_model.integral_reltol,
        abstol=albedo_model.integral_abstol,
        integral_cache=albedo_model.integral_cache,
    )
end

"""
    albedo_accel(
        u::AbstractVector, 
        sun_pos::AbstractVector, 
        RC::Number, 
        current_time::Number,
        body_albedo_model::AbstractAlbedoModel,
        eop_data::Union{EopIau1980,EopIau2000A};
        solar_flux::Number=1360.8,  # Solar flux at 1 AU [W/m²]
        speed_of_light::Number=SPEED_OF_LIGHT,  # Speed of light [m/s]
        AU::Number=149597870.7,  # Astronomical unit [km]
        reltol::Number=1e-4,
        abstol::Number=1e-8
    )

Compute the acceleration from Earth albedo radiation pressure using numerical integration.

Earth albedo radiation pressure arises from two sources:
1. Solar radiation reflected by Earth's surface (shortwave, albedo component)
2. Thermal radiation emitted by Earth (longwave, infrared component)

The total acceleration is computed by integrating over the Earth's surface visible 
to the satellite (field of view), considering both reflected and emitted radiation.

Mathematical formulation based on Knocke et al. (1988):
    a_albedo = (A/m) * RC * ∫∫ [F_SW + F_LW] * cos(γ) * cos(α) * dΩ / (π * c * r²)

Where:
- A/m: Area-to-mass ratio (embedded in RC)
- RC: Reflectivity coefficient  
- F_SW: Shortwave flux (reflected solar radiation)
- F_LW: Longwave flux (thermal emission)
- γ: Angle between incident radiation and satellite surface normal
- α: Angle between Earth surface element normal and satellite direction
- dΩ: Surface element area
- c: Speed of light
- r: Distance from surface element to satellite

# Arguments
- `u::AbstractVector`: The current state of the spacecraft in the central body inertial frame.
- `sun_pos::AbstractVector`: The current position of the Sun [km].
- `RC::Number`: The reflectivity ballistic coefficient of the satellite [m²/kg].
- `body_albedo_model::AbstractAlbedoModel{<:Number, <:Number}`: Earth albedo radiation model.

# Optional Arguments
- `solar_flux::Number`: Solar flux at 1 AU [W/m²].
- `speed_of_light::Number`: Speed of light [km/s].
- `eop_data::Union{EopIau1980,EopIau2000A,Nothing}`: Earth orientation parameters.
- `current_time::Number`: Current time [Julian days].
- `reltol::Number`: Relative tolerance for numerical integration.
- `abstol::Number`: Absolute tolerance for numerical integration.

# Returns
- `SVector{3}{Number}`: Inertial acceleration from albedo radiation pressure [km/s²].
"""
@inline function albedo_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector{ST},
    RC::RCT,
    current_time::TT,
    body_albedo_model::BAM,
    eop_data::EoT;
    solar_flux::SFT=SOLAR_FLUX,
    AU::AUT=ASTRONOMICAL_UNIT,
    speed_of_light::CT=SPEED_OF_LIGHT,
    integral_algorithm::ALG=HCubatureJL(),
    reltol::IRT=1e-4,
    abstol::IAT=1e-8,
    integral_cache::CAT=Ref(nothing),
) where {
    UT<:Number,
    ST<:Number,
    RCT<:Number,
    TT<:Number,
    AT<:Number,
    ET<:Number,
    BAM<:AbstractAlbedoModel{AT,ET},
    EoT<:Union{EopIau1980,EopIau2000A},
    SFT<:Number,
    AUT<:Number,
    CT<:Number,
    ALG<:SciMLBase.AbstractIntegralAlgorithm,
    IRT<:AbstractFloat,
    IAT<:AbstractFloat,
    CAT<:Ref{Union{Nothing, Integrals.IntegralCache}},
}

    RT = promote_type(UT, ST, RCT, TT, AT, ET, SFT, AUT, CT)
    
    sat_pos = SVector{3,RT}(u[1], u[2], u[3])
    
    R_ECEF2ECI = r_ecef_to_eci(ITRF(), J2000(), current_time, eop_data)

    albedo_functor = AlbedoFunctor(sat_pos, R_ECEF2ECI, sun_pos, RC, body_albedo_model, current_time, solar_flux, speed_of_light, AU)
    
    # Create integration problem using optimized function (no struct allocation)
    prob = IntegralProblem(albedo_functor, ALBEDO_INTEGRATION_DOMAIN)
    
    # Use cache if provided, otherwise solve normally
    if integral_cache !== nothing && integral_cache[] !== nothing
        # Update the cache with the new problem and solve
        cache = integral_cache[]
        cache.f = prob.f
        cache.domain = prob.domain
        sol = solve!(cache)
    else
        sol = solve(prob, integral_algorithm; reltol=reltol, abstol=abstol)
    end
    
    return sol.u
end

"""
    compute_earth_radiation_fluxes(
        body_albedo_model::AbstractAlbedoModel, 
        lat::Number, 
        lon::Number, 
        cos_solar_zenith::Number, 
        solar_flux_at_earth::Number,
        current_time::Number
    )

Compute the shortwave (reflected) and longwave (thermal) radiation fluxes from a surface element.

# Arguments
- `body_albedo_model::AbstractAlbedoModel{<:Number, <:Number}`: Earth albedo model containing reflection/emission parameters.
- `lat::Number`: Latitude of surface element [radians].
- `lon::Number`: Longitude of surface element [radians].
- `cos_solar_zenith::Number`: Cosine of solar zenith angle at surface element.
- `solar_flux_at_earth::Number`: Solar flux incident at Earth's surface [W/m²].
- `current_time::Number`: Current time [Julian days].

# Returns
- `(shortwave_flux, longwave_flux)`: Tuple of shortwave and longwave fluxes [W/m²].
"""
function compute_earth_radiation_fluxes(
    body_albedo_model::UniformAlbedoModel{AT, ET}, 
    lat::Number, 
    lon::Number, 
    cos_solar_zenith::Number, 
    solar_flux_at_earth::Number,
    current_time::Number
) where {AT<:Number, ET<:Number}

    # For uniform model, use constant albedo and emissivity
    albedo = body_albedo_model.visible_albedo
    emissivity = body_albedo_model.infrared_emissivity
    
    # Shortwave flux (reflected solar radiation)
    # Only consider daytime (cos_solar_zenith > 0)
    if cos_solar_zenith > 0
        shortwave_flux = albedo * solar_flux_at_earth * cos_solar_zenith
    else
        shortwave_flux = 0.0
    end
    
    # Longwave flux (thermal emission) 
    # Simplified: constant emission proportional to average incoming solar flux
    # More sophisticated models would use temperature calculations
    average_solar_flux = solar_flux_at_earth * 0.25  # 1/4 factor for spherical averaging
    longwave_flux = emissivity * average_solar_flux
    
    return shortwave_flux + longwave_flux
end