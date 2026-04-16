# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Different Shadow Models used Mainly in SRP Calculation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Montenbruck, O., & Gill, E. (2000). Satellite Orbits: Models, Methods, and Applications. Springer.
#   [2] Aziz, J., et al. (2019). "A Smoothed Eclipse Model for Solar Electric Propulsion Trajectory
#       Optimization." Trans. JSASS Aerospace Tech. Japan, 17(2), 181-188. doi:10.2322/tastj.17.181
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export ShadowModelType, Conical, Cylindrical, NoShadow, SmoothedConical
abstract type ShadowModelType end
struct Conical <: ShadowModelType end
struct Cylindrical <: ShadowModelType end
struct NoShadow <: ShadowModelType end

"""
    SmoothedConical{CST,CTT} <: ShadowModelType

A smoothed (differentiable) shadow model based on Aziz et al. [2] that uses a logistic
sigmoid function to provide a continuous transition between sunlight and shadow.

This model is fully compatible with automatic differentiation since it avoids
branching (if/else) on shadow state, making it ideal for gradient-based optimization.

# Fields
- `cs::CST`: Sharpness coefficient controlling the slope of the transition (default: 289.78)
- `ct::CTT`: Transition coefficient controlling where the curve is centered (default: 1.0)

The default values of `cs = 289.78` and `ct = 1.0` were found by Aziz et al. to minimize
error relative to the discontinuous eclipse model for geocentric orbits.
"""
struct SmoothedConical{CST<:Number,CTT<:Number} <: ShadowModelType
    cs::CST
    ct::CTT
end
SmoothedConical() = SmoothedConical(289.78, 1.0)

# Zero-position default for the central-body (body-at-origin) shadow path.
# Exposed as an internal constant so callers can pass custom body positions
# for additional occulting bodies (see `OccultingBody`).
const _ORIGIN_SV3 = SVector{3,Float64}(0.0, 0.0, 0.0)

export shadow_model
"""
    shadow_model(sat_pos, sun_pos, ShadowModel; R_Sun, R_Occulting, body_pos)

Compute the lighting factor of the Sun occluded by a body.

By default the occulting body is assumed to sit at the origin of the propagation frame
(i.e. the central body). Passing `body_pos` places the occulting body at an arbitrary
inertial position, which is used by the multi-occulter path in `SRPAstroModel` and
`ThermalEmissionAstroModel`.

# Arguments
- `sat_pos::AbstractVector`: Spacecraft position in the propagation frame [km].
- `sun_pos::AbstractVector`: Sun position in the propagation frame [km].
- `ShadowModel::ShadowModelType`: Shadow model variant. One of `Cylindrical`,
  `Conical`, `SmoothedConical`, `NoShadow`.

# Keyword Arguments
- `R_Sun::Number`: Radius of the Sun [km]. Default: `R_SUN`.
- `R_Occulting::Number`: Radius of the occulting body [km]. Default: `R_EARTH`.
- `body_pos::AbstractVector`: Inertial position of the occulting body in the propagation
  frame [km]. Default: zero (body at origin = central body).

# Returns
- `Number`: Shadow factor between 0.0 (full shadow) and 1.0 (full sunlight).
"""
@inline function shadow_model(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    ShadowModel::Cylindrical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    body_pos::AbstractVector=_ORIGIN_SV3,
)
    # Spacecraft position relative to the occulting body
    r_rel = SVector{3}(
        sat_pos[1] - body_pos[1], sat_pos[2] - body_pos[2], sat_pos[3] - body_pos[3]
    )
    # Sun direction as seen from the occulting body (defines the shadow axis)
    sun_from_body = SVector{3}(
        sun_pos[1] - body_pos[1], sun_pos[2] - body_pos[2], sun_pos[3] - body_pos[3]
    )
    sun_direction = sun_from_body / norm(sun_from_body)

    dp_sun_sat = dot(sun_direction, r_rel)

    if dp_sun_sat >= 0.0 || norm(r_rel - dp_sun_sat * sun_direction) > R_Occulting
        shadow_factor = 1.0
    else
        shadow_factor = 0.0
    end

    return shadow_factor
end

@inline function shadow_model(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    ShadowModel::Conical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    body_pos::AbstractVector=_ORIGIN_SV3,
)

    # Montenbruck, Oliver, Eberhard Gill, and F. H. Lutze. "Satellite orbits: models, methods, and applications." Appl. Mech. Rev. 55.2 (2002): B27-B28.
    # https://link.springer.com/book/10.1007/978-3-642-58351-3
    # Section 3.4.2
    #
    # Generalised for a body at `body_pos`: replace `sat_pos` with `sat_pos - body_pos`
    # in the occluding-body geometry.

    R_spacecraft_Sun = SVector{3}(
        sat_pos[1] - sun_pos[1], sat_pos[2] - sun_pos[2], sat_pos[3] - sun_pos[3]
    )
    r_rel = SVector{3}(
        sat_pos[1] - body_pos[1], sat_pos[2] - body_pos[2], sat_pos[3] - body_pos[3]
    )

    a = asin(R_Sun / norm(R_spacecraft_Sun))
    b = asin(R_Occulting / norm(r_rel))

    c = angle_between_vectors(R_spacecraft_Sun, r_rel)

    if c ≥ (b + a)
        shadow_factor = 1.0
    elseif c < (b - a)
        shadow_factor = 0.0
    elseif c < (a - b)
        shadow_factor = 1.0 - (b^2.0) / (a^2.0)
    else
        x = (c^2.0 + a^2.0 - b^2.0) / (2.0 * c)
        y = √(a^2.0 - x^2.0)
        area = a^2.0 * acos(x / a) + b^2.0 * acos((c - x) / b) - c * y
        shadow_factor = 1.0 - area / (π * a^2.0)
    end

    return shadow_factor
end

"""
    shadow_model(sat_pos, sun_pos, ::SmoothedConical; R_Sun, R_Occulting, body_pos)

Compute the sunlight fraction using the smoothed eclipse model from Aziz et al. [2].

The model computes apparent angular semi-diameters of the Sun and occulting body as seen
from the spacecraft, then applies a logistic sigmoid to produce a smooth, differentiable
transition between full sunlight (γ ≈ 1) and full shadow (γ ≈ 0).

# Equations (Aziz et al. Eq. 16-19)
- `a_SR = asin(R_Sun / ||r_sun/sc||)` — apparent angular radius of the Sun
- `a_BR = asin(R_B / ||r_B/sc||)` — apparent angular radius of the occulting body
- `a_D = acos(r̂_B/sc · r̂_sun/sc)` — angular separation
- `γ = 1 / (1 + exp(-cs * [a_D - ct * (a_SR + a_BR)]))` — sunlight fraction
"""
@inline function shadow_model(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    model::SmoothedConical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    body_pos::AbstractVector=_ORIGIN_SV3,
)
    # Position of the Sun relative to the spacecraft
    r_sun_sc = SVector{3}(
        sun_pos[1] - sat_pos[1], sun_pos[2] - sat_pos[2], sun_pos[3] - sat_pos[3]
    )
    # Position of the occulting body relative to the spacecraft
    r_body_sc = SVector{3}(
        body_pos[1] - sat_pos[1], body_pos[2] - sat_pos[2], body_pos[3] - sat_pos[3]
    )

    # Apparent angular semi-diameters as seen from the spacecraft
    a_SR = asin(R_Sun / norm(r_sun_sc))
    a_BR = asin(R_Occulting / norm(r_body_sc))

    # Angular separation between Sun and occulting body
    a_D = angle_between_vectors(r_body_sc, r_sun_sc)

    # Smoothed sunlight fraction via logistic sigmoid
    γ = 1 / (1 + exp(-model.cs * (a_D - model.ct * (a_SR + a_BR))))

    return γ
end

@inline function shadow_model(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    ShadowModel::NoShadow;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
    body_pos::AbstractVector=_ORIGIN_SV3,
)
    return one(eltype(sat_pos))
end

# ==========================================================================================
# Multi-body shadow support
# ==========================================================================================

export OccultingBody

"""
    OccultingBody{TBT<:ThirdBodyModel, RT<:Number}

A lightweight wrapper pairing a `ThirdBodyModel` (which supplies the body's position via
its `FrameEphemeris` + compiled translation closure) with its physical radius.

Used as the element type of the `additional_occulters` tuple on
[`SRPAstroModel`](@ref) and [`ThermalEmissionAstroModel`](@ref) to add shadowing from
bodies other than the central (propagation) body — e.g. the Moon when propagating in
Earth-centered frames, or Jupiter when propagating relative to a Jovian moon.

# Fields
- `body::ThirdBodyModel`: Provides the occulter's position at evaluation time.
- `radius::Number`: Physical radius of the occulter [km].

# Example
```julia
moon_occulter = OccultingBody(
    ThirdBodyModel(;
        body=MoonBody(),
        ephem_type=FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF),
        frames=my_frames,
    ),
    AstroForceModels.R_MOON,
)
```
"""
struct OccultingBody{TBT<:ThirdBodyModel,RT<:Number}
    body::TBT
    radius::RT
end

# -- Tuple-recursive resolver: converts a tuple of OccultingBody into a tuple of
#    pre-resolved (body_pos, radius) pairs, to decouple ephemeris lookup from the
#    purely-numeric shadow product helper. --

@inline function _resolve_occulters(occs::Tuple, frames, t_ft)
    first_occ = first(occs)
    body_pos = get_position(
        first_occ.body.ephem_type,
        first_occ.body.body,
        frames,
        t_ft,
        first_occ.body.compiled_vector3,
    )
    return (
        (body_pos, first_occ.radius), _resolve_occulters(Base.tail(occs), frames, t_ft)...
    )
end

@inline _resolve_occulters(::Tuple{}, frames, t_ft) = ()

# -- Tuple-recursive shadow-factor product over an already-resolved occulter tuple. --

@inline function _shadow_prod(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    shadow_m::ShadowModelType,
    R_Sun::Number,
    occs::Tuple,
)
    body_pos, R_occ = first(occs)
    F = shadow_model(
        sat_pos, sun_pos, shadow_m; R_Sun=R_Sun, R_Occulting=R_occ, body_pos=body_pos
    )
    return F * _shadow_prod(sat_pos, sun_pos, shadow_m, R_Sun, Base.tail(occs))
end

@inline function _shadow_prod(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    shadow_m::ShadowModelType,
    R_Sun::Number,
    ::Tuple{},
)
    return one(eltype(sat_pos))
end
