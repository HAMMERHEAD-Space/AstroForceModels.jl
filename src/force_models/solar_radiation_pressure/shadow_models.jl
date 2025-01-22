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
#TODO: REFERENCE
#   [1]
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export Conical, Cylindrical, No_Shadow
abstract type ShadowModelType end
struct Conical <: ShadowModelType end
struct Cylindrical <: ShadowModelType end
struct No_Shadow <: ShadowModelType end

export shadow_model
"""
    shadow_model(sat_pos::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Occulting::Number, t::Number, ShadowModel::ShadowModelType)
Computes the Lighting Factor of the Sun occur from the Umbra and Prenumbra of Earth's Shadow

# Arguments

- `sat_pos::AbstractArray`: The current satellite position.
- `sun_pos::AbstractArray`: The current Sun position.
- `R_Sun::Number`: The radius of the Sun.
- `R_Occulting::Number`: The radius of the Occulting Body.
- `ShadowModel::ShadowModelType`: The Earth shadow model to use. Current Options -- Cylindrical, Conical, Conical_Simplified, None

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Cylindrical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
)
    sat_pos = SVector{3}(sat_pos[1:3])
    sun_direction = SVector{3}(normalize(sun_pos))

    # Compute dot product between sun and satellite positions
    dp_sun_sat = dot(sun_direction, sat_pos)

    if dp_sun_sat >= 0.0 || norm(sat_pos - dp_sun_sat * sun_direction) > R_Occulting
        shadow_factor = 1.0
    else
        shadow_factor = 0.0
    end

    return shadow_factor
end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Conical;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
)

    # Montenbruck, Oliver, Eberhard Gill, and F. H. Lutze. "Satellite orbits: models, methods, and applications." Appl. Mech. Rev. 55.2 (2002): B27-B28.
    # https://link.springer.com/book/10.1007/978-3-642-58351-3
    # Section 3.4.2

    R_spacecraft_Sun = SVector{3}(
        sat_pos[1] - sun_pos[1],
        sat_pos[2] - sun_pos[2],
        sat_pos[3] - sun_pos[3],
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

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::No_Shadow;
    R_Sun::Number=R_SUN,
    R_Occulting::Number=R_EARTH,
)
    return 1.0
end
