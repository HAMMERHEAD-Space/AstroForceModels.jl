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

export shadow_model
"""
    shadow_model(sat_pos::AbstractVector, sun_pos::AbstractVector, R_Sun::Number, R_Occulting::Number, t::Number, ShadowModel::ShadowModelType)
Computes the Lighting Factor of the Sun occur from the Umbra and Prenumbra of Earth's Shadow

# Arguments

- `sat_pos::AbstractVector`: The current satellite position.
- `sun_pos::AbstractVector`: The current Sun position.
- `R_Sun::Number`: The radius of the Sun.
- `R_Occulting::Number`: The radius of the Occulting Body.
- `ShadowModel::ShadowModelType`: The Earth shadow model to use. Current Options -- Cylindrical, Conical, SmoothedConical, NoShadow

# Returns

- `Number`: Shadow factor between 0.0 (full shadow) and 1.0 (full sunlight)
"""
@inline function shadow_model(
    sat_pos::AbstractVector,
    sun_pos::AbstractVector,
    ShadowModel::Cylindrical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
)
    _sat_pos = SVector{3}(sat_pos[1], sat_pos[2], sat_pos[3])
    sun_direction = SVector{3}(normalize(sun_pos))

    # Compute dot product between sun and satellite positions
    dp_sun_sat = dot(sun_direction, _sat_pos)

    if dp_sun_sat >= 0.0 || norm(_sat_pos - dp_sun_sat * sun_direction) > R_Occulting
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
)

    # Montenbruck, Oliver, Eberhard Gill, and F. H. Lutze. "Satellite orbits: models, methods, and applications." Appl. Mech. Rev. 55.2 (2002): B27-B28.
    # https://link.springer.com/book/10.1007/978-3-642-58351-3
    # Section 3.4.2

    R_spacecraft_Sun = SVector{3}(
        sat_pos[1] - sun_pos[1], sat_pos[2] - sun_pos[2], sat_pos[3] - sun_pos[3]
    )

    a = asin(R_Sun / norm(R_spacecraft_Sun))
    b = asin(R_Occulting / norm(sat_pos))

    c = angle_between_vectors(R_spacecraft_Sun, sat_pos)

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
    shadow_model(sat_pos, sun_pos, ::SmoothedConical; R_Sun, R_Occulting)

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
)
    # Position of the Sun relative to the spacecraft
    r_sun_sc = SVector{3}(
        sun_pos[1] - sat_pos[1], sun_pos[2] - sat_pos[2], sun_pos[3] - sat_pos[3]
    )

    # Position of the occulting body (at origin) relative to the spacecraft
    r_body_sc = SVector{3}(-sat_pos[1], -sat_pos[2], -sat_pos[3])

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
)
    return one(eltype(sat_pos))
end
