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

Surface positions are pre-computed at construction time on a spherical Earth using Lebedev
quadrature points, eliminating expensive geodetic conversions from the integration loop.

# Fields
- `satellite_shape_model::AbstractSatelliteSRPModel`: The satellite shape model for computing the ballistic coefficient.
- `sun_data::ThirdBodyModel`: The data to compute the Sun's position.
- `body_albedo_model::AbstractAlbedoModel{<:Number, <:Number}`: The Earth albedo radiation model.
- `eop_data::EopIau1980`: Earth orientation parameters.
- `solar_irradiance::Number`: Solar irradiance at 1 AU [W/m²].
- `speed_of_light::Number`: Speed of light [km/s].
- `surface_positions_ecef::Vector{SVector{3,Float64}}`: Pre-computed surface element positions in ECEF [km].
- `weights::Vector{Float64}`: Scaled Lebedev quadrature weights.
"""
struct AlbedoAstroModel{ST,SDT,EAT,EoT,SFT,CT,AUT,PT,WT} <:
       AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel,
    EoT<:EopIau1980,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    PT<:AbstractVector{<:AbstractVector{<:Number}},
    WT<:AbstractVector{<:Number},
}
    satellite_shape_model::ST
    sun_data::SDT
    body_albedo_model::EAT
    eop_data::EoT

    solar_irradiance::SFT
    speed_of_light::CT
    AU::AUT

    # Pre-computed at construction time
    surface_positions_ecef::PT  # Surface element positions in ECEF [km]
    weights::WT                 # Scaled integration weights (4π × Lebedev weights)
end

# Constructor
function AlbedoAstroModel(
    satellite_shape_model::ST,
    sun_data::SDT,
    body_albedo_model::EAT,
    eop_data::EoT;
    solar_irradiance::SFT=SOLAR_IRRADIANCE,
    speed_of_light::CT=SPEED_OF_LIGHT,
    AU::AUT=ASTRONOMICAL_UNIT / 1E3,
    lebedev_order::Int=125,
    radius::T=R_EARTH,
) where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel,
    EoT<:EopIau1980,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    T<:Number,
}
    x_coords, y_coords, z_coords, lebedev_weights = lebedev_by_order(lebedev_order)
    n_points = length(lebedev_weights)

    # Pre-compute surface positions in ECEF using spherical Earth approximation.
    # Lebedev (x,y,z) are unit vectors on the sphere, so R_EARTH * point gives
    # the surface position directly -- no geodetic conversion needed at runtime.
    surface_positions_ecef = [
        SVector{3,T}(radius * x_coords[i], radius * y_coords[i], radius * z_coords[i]) for
        i in 1:n_points
    ]

    scaled_weights = 4π .* lebedev_weights

    return AlbedoAstroModel{
        ST,SDT,EAT,EoT,SFT,CT,AUT,typeof(surface_positions_ecef),typeof(scaled_weights)
    }(
        satellite_shape_model,
        sun_data,
        body_albedo_model,
        eop_data,
        solar_irradiance,
        speed_of_light,
        AU,
        surface_positions_ecef,
        scaled_weights,
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
        albedo_model.surface_positions_ecef,
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
        surface_positions_ecef::AbstractVector,
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
Surface element positions are pre-computed in ECEF at construction time and rotated
to ECI once per evaluation step, using dot products for all angle computations.

Mathematical formulation based on Knocke et al. (1988):
    a_albedo = RC * ∫∫ [F_SW + F_LW] * cos(α) * dΩ / (π * c * r²) / 1E3

# Arguments
- `u::AbstractVector`: The current state of the spacecraft in the central body inertial frame [km, km/s].
- `sun_pos::AbstractVector`: The current position of the Sun [km].
- `RC::Number`: The reflectivity ballistic coefficient of the satellite [m²/kg].
- `current_time::Number`: The current Julian date.
- `body_albedo_model::AbstractAlbedoModel`: Earth albedo radiation model.
- `eop_data::EopIau1980`: Earth orientation parameters.
- `surface_positions_ecef::AbstractVector`: Pre-computed ECEF surface element positions [km].
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
    surface_positions_ecef::SPT,
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
    SPT<:AbstractVector{<:AbstractVector{<:Number}},
    WT<:AbstractVector{<:Number},
    SFT<:Number,
    AUT<:Number,
    CT<:Number,
}
    RT = promote_type(UT, ST, RCT, TT, SFT, AUT, CT)

    sat_pos = SVector{3,RT}(u[1], u[2], u[3])

    R_ECEF2ECI = r_ecef_to_eci(ITRF(), J2000(), current_time, eop_data)

    inv_c_mps = 1 / (speed_of_light * 1E3)
    inv_pi = 1 / π
    inv_1E3 = 1 / 1E3

    n_points = length(weights)
    z = zero(RT)
    result = SVector{3,RT}(z, z, z)

    @inbounds for i in 1:n_points
        sp = R_ECEF2ECI * surface_positions_ecef[i]
        result += _albedo_surface_element(
            sat_pos,
            sp,
            sun_pos,
            RC,
            weights[i],
            body_albedo_model,
            solar_irradiance,
            AU,
            inv_c_mps,
            inv_pi,
            inv_1E3,
        )
    end

    return result
end

@inline function _albedo_surface_element(
    sat_pos::AbstractVector{<:Number},
    sp::AbstractVector{<:Number},
    sun_pos::AbstractVector{<:Number},
    RC::Number,
    weight::Number,
    body_albedo_model::AbstractAlbedoModel,
    solar_irradiance::Number,
    AU::Number,
    inv_c_mps::Number,
    inv_pi::Number,
    inv_1E3::Number,
)
    R_e = norm(sp)
    inv_R_e = 1 / R_e

    d = sat_pos - sp
    dist_sq = d[1] * d[1] + d[2] * d[2] + d[3] * d[3]
    dist = √(dist_sq)

    cos_alpha = (sp[1] * d[1] + sp[2] * d[2] + sp[3] * d[3]) * inv_R_e / dist
    if cos_alpha <= 0
        z = zero(dist)
        return SVector{3}(z, z, z)
    end

    s = sun_pos - sp
    dist_sun = √(s[1] * s[1] + s[2] * s[2] + s[3] * s[3])
    cos_zenith = (sp[1] * s[1] + sp[2] * s[2] + sp[3] * s[3]) * inv_R_e / dist_sun

    irradiance = solar_irradiance * (AU / dist_sun)^2

    total_flux = compute_earth_radiation_fluxes(body_albedo_model, cos_zenith, irradiance)

    coeff =
        weight * RC * total_flux * inv_c_mps * R_e^2 * cos_alpha * inv_pi / dist_sq *
        inv_1E3

    return (coeff / dist) * d
end

"""
    compute_earth_radiation_fluxes(
        body_albedo_model::UniformAlbedoModel, 
        cos_solar_zenith::Number, 
        irradiance_at_earth::Number
    )

Compute the total radiation flux (shortwave + longwave) from a surface element.

# Arguments
- `body_albedo_model::UniformAlbedoModel`: Earth albedo model with uniform coefficients.
- `cos_solar_zenith::Number`: Cosine of solar zenith angle at surface element.
- `irradiance_at_earth::Number`: Solar irradiance at Earth's distance [W/m²].

# Returns
- `Number`: Total outgoing flux [W/m²].
"""
@inline function compute_earth_radiation_fluxes(
    body_albedo_model::UniformAlbedoModel{AT,ET},
    cos_solar_zenith::Number,
    irradiance_at_earth::Number,
) where {AT<:Number,ET<:Number}
    albedo = body_albedo_model.visible_albedo
    emissivity = body_albedo_model.infrared_emissivity

    # Shortwave (reflected): only on sunlit side
    shortwave = if cos_solar_zenith > 0
        albedo * irradiance_at_earth * cos_solar_zenith
    else
        zero(irradiance_at_earth)
    end

    # Longwave (thermal): isotropic, scaled by spherical averaging factor 1/4
    longwave = emissivity * irradiance_at_earth * 0.25

    return shortwave + longwave
end
