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

export AlbedoAstroModel, albedo_accel
export AbstractAlbedoModel, UniformAlbedoModel

# Set up integration domain in radians
# Latitude: -π/2 to π/2, Longitude: -π to π
const ALBEDO_INTEGRATION_DOMAIN = (SVector{2}(-π/2, -π), SVector{2}(π/2, π))

"""
Abstract type for albedo radiation models used in albedo force calculations.
"""
abstract type AbstractAlbedoModel{AT<:Number,ET<:Number} end

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
Base.@kwdef struct UniformAlbedoModel{AT,ET} <: AbstractAlbedoModel{AT,ET}
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
- `eop_data::EopIau1980`: Earth orientation parameters.
- `solar_irradiance::Number`: Solar irradiance at 1 AU [W/m²].
- `speed_of_light::Number`: Speed of light [km/s].
- `lebedev_order::Int`: Order of Lebedev quadrature for spherical integration.
"""
mutable struct AlbedoAstroModel{ST,SDT,EAT,EoT,SFT,CT,AUT,PT,WT,TT,PhT} <:
               AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel,
    EoT<:EopIau1980,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    PT<:AbstractVector,
    WT<:AbstractVector,
    TT<:AbstractVector,
    PhT<:AbstractVector,
}
    satellite_shape_model::ST
    sun_data::SDT
    body_albedo_model::EAT
    eop_data::EoT

    solar_irradiance::SFT
    speed_of_light::CT
    AU::AUT

    # Lebedev quadrature points and weights for unit sphere
    lebedev_points::PT    # 3D points on unit sphere
    lebedev_weights::WT   # Integration weights

    # Pre-computed spherical coordinates for efficiency
    theta_coords::TT      # Colatitude coordinates [0, π]
    phi_coords::PhT       # Azimuth coordinates [0, 2π]  
    weights::WT           # Scaled integration weights
end

# Constructor
function AlbedoAstroModel(;
    satellite_shape_model,
    sun_data,
    body_albedo_model,
    eop_data,
    solar_irradiance=SOLAR_IRRADIANCE,
    speed_of_light=SPEED_OF_LIGHT,
    AU=ASTRONOMICAL_UNIT / 1E3,
    lebedev_order::Int=125,  # Order of Lebedev quadrature
)
    # Generate Lebedev quadrature nodes and weights for unit sphere
    x_coords, y_coords, z_coords, lebedev_weights = lebedev_by_order(lebedev_order)

    # Combine coordinates into points vector
    n_points = length(lebedev_weights)
    lebedev_points = [
        SVector{3,Float64}(x_coords[i], y_coords[i], z_coords[i]) for i in 1:n_points
    ]

    # Convert Lebedev (x,y,z) points to spherical coordinates (θ, φ)
    # θ ∈ [0, π] (colatitude), φ ∈ [0, 2π] (azimuth)
    theta_coords = Vector{Float64}(undef, n_points)
    phi_coords = Vector{Float64}(undef, n_points)

    @inbounds for i in 1:n_points
        x, y, z = lebedev_points[i]
        theta_coords[i] = acos(clamp(z, -1.0, 1.0))  # θ = acos(z)
        phi_coords[i] = atan(y, x)  # φ = atan2(y, x)
        # Ensure φ ∈ [0, 2π]
        phi_coords[i] = rem2pi(phi_coords[i], RoundDown)
    end

    # Lebedev weights already include 4π normalization for unit sphere
    # Scale by 4π to get proper integration weights
    scaled_weights = 4π .* lebedev_weights

    return AlbedoAstroModel{
        typeof(satellite_shape_model),
        typeof(sun_data),
        typeof(body_albedo_model),
        typeof(eop_data),
        typeof(solar_irradiance),
        typeof(speed_of_light),
        typeof(AU),
        typeof(lebedev_points),
        typeof(lebedev_weights),
        typeof(theta_coords),
        typeof(phi_coords),
    }(
        satellite_shape_model,
        sun_data,
        body_albedo_model,
        eop_data,
        solar_irradiance,
        speed_of_light,
        AU,
        lebedev_points,
        lebedev_weights,
        theta_coords,
        phi_coords,
        scaled_weights,
    )
end

function albedo_integrand(
    x::AbstractVector{XT},
    sat_pos::AbstractVector{ST},
    R_ECEF2ECI::DCM{RT},
    sun_pos::AbstractVector{SUT},
    RC::RCT,
    body_albedo_model::AM,
    current_time::TT,
    solar_irradiance::SFT,
    speed_of_light::CT,
    AU::AUT,
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
    AM<:AbstractAlbedoModel{AT,ET},
}
    RET = promote_type(XT, ST, RT, SUT, RCT, TT, SFT, CT, AUT, AT, ET)

    theta_rad, phi_rad = x  # θ ∈ [0, π] (colatitude), φ ∈ [0, 2π] (azimuth)

    # Convert spherical coordinates to geodetic coordinates
    # θ = 0 is North pole, θ = π is South pole
    lat_rad = π/2 - theta_rad  # Convert colatitude to latitude: lat ∈ [-π/2, π/2]
    lon_rad = phi_rad > π ? phi_rad - 2π : phi_rad  # Convert to longitude: lon ∈ [-π, π]

    # Surface element position in ECEF coordinates using SatelliteToolboxTransformations
    surface_pos_ecef = geodetic_to_ecef(lat_rad, lon_rad, 0.0) ./ 1E3 # Convert meters to km

    # Transform surface position to ECI
    surface_pos = R_ECEF2ECI * surface_pos_ecef

    # Vector from surface element to satellite
    surface_to_sat = sat_pos - surface_pos
    distance = norm(surface_to_sat)

    # Angle between surface normal and satellite direction
    angle_from_nadir = angle_between_vectors(surface_pos, surface_to_sat)
    cos_alpha = cos(angle_from_nadir)

    # Skip if surface element is not visible (behind horizon)
    if cos_alpha <= 0 || angle_from_nadir > π / 2
        return SVector{3,RET}(0.0, 0.0, 0.0)
    end

    # Surface element area: R² (sin(θ) factor included in Lebedev weights)
    R_earth = norm(surface_pos)
    jacobian = R_earth^2  # km²

    # Vector from surface element to Sun
    surface_to_sun = sun_pos - surface_pos
    surface_to_sun_norm = norm(surface_to_sun)

    # Solar irradiance at Earth's distance [W/m²]
    irradiance_at_earth = solar_irradiance * (AU / surface_to_sun_norm)^2

    # Solar zenith angle at surface element
    solar_zenith_angle = angle_between_vectors(surface_pos, surface_to_sun)
    cos_solar_zenith = cos(solar_zenith_angle)

    # Compute shortwave and longwave fluxes [W/m²]
    total_radiance = compute_earth_radiation_fluxes(
        body_albedo_model,
        lat_rad,
        lon_rad,
        cos_solar_zenith,
        irradiance_at_earth,
        current_time,
    )

    # Convert irradiance [W/m²] to radiation pressure [N/m²] by dividing by c [m/s]
    # speed_of_light is in [km/s], multiply by 1E3 to get [m/s]
    pressure = total_radiance / (speed_of_light * 1E3)  # [N/m²]

    # Geometric view factor: solid angle of surface element as seen from satellite
    geometric_factor = jacobian * cos_alpha / (π * distance^2)  # km² / km² = dimensionless

    # Acceleration magnitude [km/s²]
    # RC [m²/kg] * pressure [N/m²] * geometric_factor [-] = [m/s²], / 1E3 → [km/s²]
    F_albedo = RC * pressure * geometric_factor / 1E3

    # Direction: surface element → satellite
    surface_to_sat_unit = surface_to_sat / distance
    return SVector{3,RET}(
        F_albedo * surface_to_sat_unit[1],
        F_albedo * surface_to_sat_unit[2],
        F_albedo * surface_to_sat_unit[3],
    )
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
    u::AbstractVector, p::ComponentVector, t::Number, albedo_model::AlbedoAstroModel
)
    jd = current_jd(p, t)

    # Compute the Sun's position (convert from m to km)
    sun_pos = albedo_model.sun_data(jd, Position()) ./ 1E3

    # Compute the reflectivity ballistic coefficient
    RC = reflectivity_ballistic_coefficient(u, p, t, albedo_model.satellite_shape_model)

    # Return the 3-dimensional Albedo Force
    return albedo_accel(
        u,
        sun_pos,
        RC,
        jd,
        albedo_model.body_albedo_model,
        albedo_model.eop_data,
        albedo_model.theta_coords,
        albedo_model.phi_coords,
        albedo_model.weights;
        solar_irradiance=albedo_model.solar_irradiance,
        AU=albedo_model.AU,
        speed_of_light=albedo_model.speed_of_light,
    )
end

"""
    albedo_accel(
        u::AbstractVector, 
        sun_pos::AbstractVector, 
        RC::Number, 
        current_time::Number,
        body_albedo_model::AbstractAlbedoModel,
        eop_data::EopIau1980,
        theta_coords::AbstractVector,
        phi_coords::AbstractVector,
        weights::AbstractVector;
        solar_irradiance::Number=SOLAR_IRRADIANCE,
        speed_of_light::Number=SPEED_OF_LIGHT,
        AU::Number=ASTRONOMICAL_UNIT / 1E3
    )

Compute the acceleration from Earth albedo radiation pressure using Lebedev quadrature.

Earth albedo radiation pressure arises from two sources:
1. Solar radiation reflected by Earth's surface (shortwave, albedo component)
2. Thermal radiation emitted by Earth (longwave, infrared component)

The total acceleration is computed by integrating over the Earth's surface visible 
to the satellite (field of view), considering both reflected and emitted radiation.

Mathematical formulation based on Knocke et al. (1988):
    a_albedo = RC * ∫∫ [F_SW + F_LW] * cos(α) * dΩ / (π * c * r²) / 1E3

Where:
- RC: Reflectivity ballistic coefficient [m²/kg]
- F_SW: Shortwave flux (reflected solar radiation) [W/m²]
- F_LW: Longwave flux (thermal emission) [W/m²]
- α: Angle between Earth surface element normal and satellite direction
- dΩ: Surface element area [km²]
- c: Speed of light [m/s]
- r: Distance from surface element to satellite [km]

# Arguments
- `u::AbstractVector`: The current state of the spacecraft in the central body inertial frame [km, km/s].
- `sun_pos::AbstractVector`: The current position of the Sun [km].
- `RC::Number`: The reflectivity ballistic coefficient of the satellite [m²/kg].
- `current_time::Number`: The current Julian date.
- `body_albedo_model::AbstractAlbedoModel`: Earth albedo radiation model.
- `eop_data::EopIau1980`: Earth orientation parameters.
- `theta_coords::AbstractVector`: Pre-computed colatitude coordinates for Lebedev quadrature.
- `phi_coords::AbstractVector`: Pre-computed azimuth coordinates for Lebedev quadrature.
- `weights::AbstractVector`: Pre-computed scaled integration weights.

# Keyword Arguments
- `solar_irradiance::Number`: Solar irradiance at 1 AU [W/m²]. Default: `SOLAR_IRRADIANCE`.
- `speed_of_light::Number`: Speed of light [km/s]. Default: `SPEED_OF_LIGHT`.
- `AU::Number`: Astronomical Unit [km]. Default: `ASTRONOMICAL_UNIT / 1E3`.

# Returns
- `SVector{3}`: Inertial acceleration from albedo radiation pressure [km/s²].
"""
function albedo_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector{ST},
    RC::RCT,
    current_time::TT,
    body_albedo_model::BAM,
    eop_data::EoT,
    theta_coords::TT_C,
    phi_coords::PT_C,
    weights::WT;
    solar_irradiance::SFT=SOLAR_IRRADIANCE,
    AU::AUT=ASTRONOMICAL_UNIT / 1E3,
    speed_of_light::CT=SPEED_OF_LIGHT,
) where {
    UT<:Number,
    ST<:Number,
    RCT<:Number,
    TT<:Number,
    BAM<:AbstractAlbedoModel,
    EoT<:EopIau1980,
    TT_C<:AbstractVector,
    PT_C<:AbstractVector,
    WT<:AbstractVector,
    SFT<:Number,
    AUT<:Number,
    CT<:Number,
}
    RT = promote_type(UT, ST, RCT, TT, SFT, AUT, CT)

    sat_pos = SVector{3,UT}(u[1], u[2], u[3])

    R_ECEF2ECI = r_ecef_to_eci(ITRF(), J2000(), current_time, eop_data)

    n_points = length(weights)

    result_x = zero(RT)
    result_y = zero(RT)
    result_z = zero(RT)

    @inbounds for i in 1:n_points
        theta = theta_coords[i]
        phi = phi_coords[i]
        weight = weights[i]

        x = SVector{2,Float64}(theta, phi)
        accel = albedo_integrand(
            x,
            sat_pos,
            R_ECEF2ECI,
            sun_pos,
            RC,
            body_albedo_model,
            current_time,
            solar_irradiance,
            speed_of_light,
            AU,
        )

        result_x += weight * accel[1]
        result_y += weight * accel[2]
        result_z += weight * accel[3]
    end

    return SVector{3,RT}(result_x, result_y, result_z)
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
    body_albedo_model::UniformAlbedoModel{AT,ET},
    lat::Number,
    lon::Number,
    cos_solar_zenith::Number,
    solar_flux_at_earth::Number,
    current_time::Number,
) where {AT<:Number,ET<:Number}

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
