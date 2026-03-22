# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Earth Albedo Radiation Pressure
#
#   NOTE: This model is inherently Earth-specific. The albedo coefficients, surface
#   positions, and quadrature are designed for Earth. Use only for Earth-orbiting spacecraft.
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
    UniformAlbedoModel

Simple uniform albedo model with constant reflection and emission coefficients.

# Fields
- `visible_albedo::Number`: Visible light albedo coefficient (0-1). Default: 0.3.
- `infrared_emissivity::Number`: Infrared emissivity coefficient (0-1). Default: 0.7.
"""
Base.@kwdef struct UniformAlbedoModel{AT,ET} <: AbstractAlbedoModel{AT,ET}
    visible_albedo::AT = 0.3
    infrared_emissivity::ET = 0.7
end

"""
    AlbedoAstroModel

Earth albedo radiation pressure force model.

Earth-only: requires ITRF frame in FrameSystem for ECEF↔ECI rotation.
Surface positions are pre-computed at construction time using Lebedev quadrature.

# Fields
- `satellite_shape_model`: Satellite SRP model for reflectivity coefficient.
- `sun_data::ThirdBodyModel`: Sun position model.
- `body_albedo_model`: Earth albedo radiation model.
- `body_fixed_frame::Symbol`: Body-fixed frame for rotating surface positions (e.g., `:ITRF`).
- `propagation_frame::Symbol`: Propagation frame (e.g., `:ICRF`).
- `solar_irradiance::Number`: Solar irradiance at 1 AU [W/m²].
- `speed_of_light::Number`: Speed of light [km/s].
- `AU::Number`: Astronomical Unit [km].
- `surface_positions_ecef`: Pre-computed ECEF surface positions [km].
- `weights`: Scaled Lebedev quadrature weights.
"""
struct AlbedoAstroModel{ST,SDT,EAT,SFT,CT,AUT,PT,WT,CR3} <:
       AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    PT<:AbstractVector{<:AbstractVector{<:Number}},
    WT<:AbstractVector{<:Number},
}
    satellite_shape_model::ST
    sun_data::SDT
    body_albedo_model::EAT
    body_fixed_frame::Symbol
    propagation_frame::Symbol

    solar_irradiance::SFT
    speed_of_light::CT
    AU::AUT

    # Pre-computed at construction time
    surface_positions_ecef::PT
    weights::WT
    compiled_rotation3::CR3
end

# Constructor
function AlbedoAstroModel(;
    satellite_shape_model::ST,
    sun_data::SDT,
    body_albedo_model::EAT,
    body_fixed_frame::Symbol=:ITRF,
    propagation_frame::Symbol=:ICRF,
    solar_irradiance::SFT=SOLAR_IRRADIANCE,
    speed_of_light::CT=SPEED_OF_LIGHT,
    AU::AUT=ASTRONOMICAL_UNIT / 1E3,
    lebedev_order::Int=125,
    radius::T=R_EARTH,
    frames=nothing,
) where {
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    EAT<:AbstractAlbedoModel,
    SFT<:Number,
    CT<:Number,
    AUT<:Number,
    T<:Number,
}
    x_coords, y_coords, z_coords, lebedev_weights = lebedev_by_order(lebedev_order)
    n_points = length(lebedev_weights)

    surface_positions_ecef = [
        SVector{3,T}(radius * x_coords[i], radius * y_coords[i], radius * z_coords[i]) for
        i in 1:n_points
    ]

    scaled_weights = 4π .* lebedev_weights

    cr3 = nothing
    if !isnothing(frames) && body_fixed_frame != propagation_frame
        cr3 = compile_rotation3(frames, body_fixed_frame, propagation_frame)
    end

    return AlbedoAstroModel{
        ST,
        SDT,
        EAT,
        SFT,
        CT,
        AUT,
        typeof(surface_positions_ecef),
        typeof(scaled_weights),
        typeof(cr3),
    }(
        satellite_shape_model,
        sun_data,
        body_albedo_model,
        body_fixed_frame,
        propagation_frame,
        solar_irradiance,
        speed_of_light,
        AU,
        surface_positions_ecef,
        scaled_weights,
        cr3,
    )
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, albedo_model::AlbedoAstroModel)

Computes the albedo acceleration acting on a spacecraft.

Earth-only. Requires ITRF frame in the FrameSystem.

# Returns
- `SVector{3}`: Albedo acceleration [km/s²].
"""
function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, albedo_model::AlbedoAstroModel
)
    t_ft = ft_time(p, t)

    # Get Sun position from FrameSystem
    sun_pos = get_position(
        albedo_model.sun_data.ephem_type,
        albedo_model.sun_data.body,
        p.frames,
        t_ft,
        albedo_model.sun_data.compiled_vector3,
    )

    # Compute the reflectivity ballistic coefficient
    RC = reflectivity_ballistic_coefficient(u, p, t, albedo_model.satellite_shape_model)

    # Return the 3-dimensional Albedo Force
    return albedo_accel(
        u,
        sun_pos,
        RC,
        p.frames,
        t_ft,
        albedo_model.body_albedo_model,
        albedo_model.body_fixed_frame,
        albedo_model.propagation_frame,
        albedo_model.surface_positions_ecef,
        albedo_model.weights;
        solar_irradiance=albedo_model.solar_irradiance,
        AU=albedo_model.AU,
        speed_of_light=albedo_model.speed_of_light,
        compiled_rotation3=albedo_model.compiled_rotation3,
    )
end

"""
    albedo_accel(u, sun_pos, RC, frames, t_ft, body_albedo_model, body_fixed_frame, propagation_frame, surface_positions_ecef, weights; kwargs...)

Compute the acceleration from Earth albedo radiation pressure using Lebedev quadrature.

Earth-only. Uses FrameSystem rotation instead of EOP data.

# Returns
- `SVector{3}`: Inertial acceleration from albedo radiation pressure [km/s²].
"""
function albedo_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector{ST},
    RC::RCT,
    frames,
    t_ft::TT,
    body_albedo_model::BAM,
    body_fixed_frame::Symbol,
    propagation_frame::Symbol,
    surface_positions_ecef::SPT,
    weights::WT;
    solar_irradiance::SFT=SOLAR_IRRADIANCE,
    AU::AUT=ASTRONOMICAL_UNIT / 1E3,
    speed_of_light::CT=SPEED_OF_LIGHT,
    compiled_rotation3=nothing,
) where {
    UT<:Number,
    ST<:Number,
    RCT<:Number,
    TT<:Number,
    BAM<:AbstractAlbedoModel,
    SPT<:AbstractVector{<:AbstractVector{<:Number}},
    WT<:AbstractVector{<:Number},
    SFT<:Number,
    AUT<:Number,
    CT<:Number,
}
    RT = promote_type(UT, ST, RCT, TT, SFT, AUT, CT)

    sat_pos = SVector{3,RT}(u[1], u[2], u[3])

    # Get rotation from body-fixed to propagation frame
    R_rot = if isnothing(compiled_rotation3)
        rotation3(frames, body_fixed_frame, propagation_frame, t_ft)
    else
        compiled_rotation3(t_ft)
    end
    R_ECEF2ECI = R_rot.m[1]

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
    compute_earth_radiation_fluxes(body_albedo_model::UniformAlbedoModel, cos_solar_zenith, irradiance_at_earth)

Compute the total radiation flux (shortwave + longwave) from a surface element.

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
