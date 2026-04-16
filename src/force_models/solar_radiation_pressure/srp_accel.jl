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
    SRPAstroModel

Solar Radiation Pressure force model.

# Fields
- `satellite_srp_model`: Satellite SRP model for reflectivity coefficient.
- `sun_data::ThirdBodyModel`: Sun position model (uses FrameEphemeris).
- `shadow_model::ShadowModelType`: Shadow model type — defaults to `Conical()`.
- `R_Sun::Number`: Radius of the Sun [km].
- `R_Occulting::Number`: Radius of the central (propagation-body) occulter [km].
  No default — must be specified.
- `additional_occulters::Tuple`: Optional tuple of [`OccultingBody`](@ref) providing
  shadowing from bodies other than the central body (e.g. the Moon when propagating
  around Earth, or Jupiter when propagating around a Jovian moon). Default: `()`
  (no additional occulters → zero overhead). Each occulter's shadow factor is
  multiplied into the primary body's, so any body that eclipses the Sun contributes.
- `Ψ::Number`: Solar flux constant at 1 AU [N/m²].
- `AU::Number`: Astronomical Unit [km].
"""
Base.@kwdef struct SRPAstroModel{
    ST<:AbstractSatelliteSRPModel,
    SDT<:ThirdBodyModel,
    SMT<:ShadowModelType,
    RST<:Number,
    ROT<:Number,
    AOT<:Tuple,
    PT<:Number,
    AUT<:Number,
} <: AbstractNonPotentialBasedForce
    satellite_srp_model::ST
    sun_data::SDT
    shadow_model::SMT = Conical()

    R_Sun::RST = R_SUN
    R_Occulting::ROT
    additional_occulters::AOT = ()
    Ψ::PT = SOLAR_FLUX
    AU::AUT = ASTRONOMICAL_UNIT / 1E3
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, srp_model::SRPAstroModel)

Computes the SRP acceleration acting on a spacecraft.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `srp_model::SRPAstroModel`: SRP model.

# Returns
- `SVector{3}`: SRP acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, srp_model::SRPAstroModel
)
    t_ft = ft_time(p, t)

    sun_pos = get_position(
        srp_model.sun_data.ephem_type,
        srp_model.sun_data.body,
        p.frames,
        t_ft,
        srp_model.sun_data.compiled_vector3,
    )

    # Resolve additional-occulter positions through the frame system once.
    # For the (common) empty-tuple case this compiles out to a no-op.
    extra_occulters = _resolve_occulters(srp_model.additional_occulters, p.frames, t_ft)

    RC = reflectivity_ballistic_coefficient(u, p, t, srp_model.satellite_srp_model)

    return srp_accel(
        u,
        sun_pos,
        RC;
        ShadowModel=srp_model.shadow_model,
        R_Sun=srp_model.R_Sun,
        R_Occulting=srp_model.R_Occulting,
        additional_occulters=extra_occulters,
        Ψ=srp_model.Ψ,
        AU=srp_model.AU,
    )
end

"""
    srp_accel(u::AbstractVector, sun_pos::AbstractVector, RC::Number;
              ShadowModel, R_Sun, R_Occulting, additional_occulters, Ψ, AU)

Compute the acceleration from Solar Radiation Pressure.

Radiation from the Sun reflects off the satellite's surface and transfers momentum perturbing the satellite's trajectory. This
force can be computed using a cannonball model with the following equation:

                𝐚 = F * RC * Ψ * (AU/(R_sc_Sun))^2 * R̂_sc_Sun

# Arguments

- `u::AbstractVector`: The current state of the spacecraft in the central body's inertial frame [km, km/s].
- `sun_pos::AbstractVector`: The current position of the Sun [km].
- `RC::Number`: The reflectivity ballistic coefficient of the satellite -- (Area/mass) * Reflectivity Coefficient [m^2/kg].

# Keyword Arguments

- `ShadowModel::ShadowModelType`: Shadow model to use. Default: `Conical()`.
- `R_Sun::Number`: The radius of the Sun [km]. Default: `R_SUN`.
- `R_Occulting::Number`: The radius of the central occulting body [km]. **Required — no default**
  (prevents silent misuse for non-Earth missions).
- `additional_occulters::Tuple`: Pre-resolved tuple of `(body_pos, radius)` pairs for
  additional shadowing bodies. Each body's shadow factor is multiplied into the
  central body's. Default: `()`.
- `Ψ::Number`: Solar radiation pressure at 1 AU [N/m^2]. Default: `SOLAR_FLUX`.
- `AU::Number`: Astronomical Unit [km]. Default: `ASTRONOMICAL_UNIT / 1E3`.

# Returns

- `SVector{3}`: Inertial acceleration from SRP [km/s^2].
"""
@inline function srp_accel(
    u::AbstractVector{UT},
    sun_pos::AbstractVector,
    RC::Number;
    ShadowModel::ShadowModelType=Conical(),
    R_Sun::Number=R_SUN,
    R_Occulting::Number,
    additional_occulters::Tuple=(),
    Ψ::Number=SOLAR_FLUX,
    AU::Number=ASTRONOMICAL_UNIT / 1E3,
) where {UT}
    sat_pos = SVector{3,UT}(u[1], u[2], u[3])

    # Shadow factor from the central (propagation-body) occulter ...
    F_primary = shadow_model(
        sat_pos, sun_pos, ShadowModel; R_Sun=R_Sun, R_Occulting=R_Occulting
    )
    # ... multiplied by the shadow factors from every additional occulter.
    # Empty-tuple path compiles to `one(eltype(u))`, no overhead.
    F = F_primary * _shadow_prod(sat_pos, sun_pos, ShadowModel, R_Sun, additional_occulters)

    # Compute the Vector Between the Satellite and Sun
    R_spacecraft_Sun = sat_pos - sun_pos
    R_sc_sun = norm(R_spacecraft_Sun)

    F_srp = F * RC * Ψ * (AU / R_sc_sun)^2 / 1E3

    #Compute the SRP Force
    return SVector{3}(
        F_srp * R_spacecraft_Sun[1] / R_sc_sun,
        F_srp * R_spacecraft_Sun[2] / R_sc_sun,
        F_srp * R_spacecraft_Sun[3] / R_sc_sun,
    )
end
