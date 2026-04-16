# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Third Body Model and Ephemeris Functions
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export FrameEphemeris

"""
    FrameEphemeris <: AbstractEphemerisType

Ephemeris type that uses FrameTransformations.jl to compute body positions and velocities.
Uses NAIF IDs and the FrameSystem to get body state vectors directly from SPK kernels.

# Fields
- `center_point::Int`: NAIF ID of the center body (e.g., 399 for Earth, 2000433 for Eros)
- `target_point::Int`: NAIF ID of the target body (e.g., 10 for Sun, 301 for Moon)
- `axes::Symbol`: Frame for output (e.g., `:ICRF`)

# Example
```julia
# Sun position relative to Earth in ICRF
sun_ephem = FrameEphemeris(center_point=399, target_point=10, axes=:ICRF)

# Moon position relative to Earth in ICRF
moon_ephem = FrameEphemeris(center_point=399, target_point=301, axes=:ICRF)

# Sun position relative to Eros in ICRF
sun_from_eros = FrameEphemeris(center_point=2000433, target_point=10, axes=:ICRF)
```
"""
Base.@kwdef struct FrameEphemeris <: AbstractEphemerisType
    center_point::Int
    target_point::Int
    axes::Symbol = :ICRF
end

export ThirdBodyModel

"""
    ThirdBodyModel{BT,EpT,CT3,CT6} <: AbstractNonPotentialBasedForce

Third body gravitational perturbation model.

# Constructor
    ThirdBodyModel(; body, ephem_type, frames=nothing)

# Arguments
- `body::CelestialBody`: Celestial body acting on the spacecraft (provides μ and radius).
- `ephem_type::FrameEphemeris`: Ephemeris type for computing body position via FrameTransformations.
  The `axes` field determines the output frame for the body's state vector.
- `frames`: Optional `FrameSystem`. When provided, pre-compiles the translation between
  `center_point` and `target_point` for allocation-free evaluation. Falls back to
  runtime `vector3`/`vector6` when `nothing` or when the point pair is not directly connected.

# Example
```julia
# Sun third-body perturbation for Earth-orbiting spacecraft
sun_model = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
    frames = my_frame_system,   # optional: enables allocation-free lookup
)

# Moon third-body perturbation
moon_model = ThirdBodyModel(
    body = MoonBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
)
```
"""
struct ThirdBodyModel{BT<:AbstractCelestialBody,EpT<:AbstractEphemerisType,CT3,CT6} <:
       AbstractNonPotentialBasedForce
    body::BT
    ephem_type::EpT
    compiled_vector3::CT3
    compiled_vector6::CT6
end

function ThirdBodyModel(;
    body::BT, ephem_type::EpT, frames=nothing
) where {BT<:AbstractCelestialBody,EpT<:AbstractEphemerisType}
    ct3 = nothing
    ct6 = nothing
    if !isnothing(frames) &&
        isa(ephem_type, FrameEphemeris) &&
        ephem_type.center_point != ephem_type.target_point
        # Pre-compile zero-overhead translations for allocation-free evaluation.
        # The new FrameTransformations API handles multi-hop paths internally, but
        # we still guard with try/catch in case the path cannot be resolved.
        try
            ct3 = compile_translation(
                frames,
                ephem_type.center_point,
                ephem_type.target_point,
                ephem_type.axes,
                Val(1),
            )
            ct6 = compile_translation(
                frames,
                ephem_type.center_point,
                ephem_type.target_point,
                ephem_type.axes,
                Val(2),
            )
        catch
            ct3 = nothing
            ct6 = nothing
        end
    end
    return ThirdBodyModel{BT,EpT,typeof(ct3),typeof(ct6)}(body, ephem_type, ct3, ct6)
end

"""
    get_position(ephem::FrameEphemeris, body::CelestialBody, frames, t_j2000)

Compute the position of a celestial body using the FrameSystem.

# Arguments
- `ephem::FrameEphemeris`: Ephemeris configuration (center, target, axes).
- `body::CelestialBody`: Celestial body (unused directly, kept for dispatch).
- `frames`: FrameSystem from FrameTransformations.jl.
- `t_j2000::Number`: Time in seconds since J2000 TDB.

# Returns
- `SVector{3}`: Position vector [km] in the specified frame.
"""
@inline function get_position(
    ephem::FrameEphemeris, body::CelestialBody, frames, t_j2000, compiled_vector3=nothing
)
    if !isnothing(compiled_vector3)
        # CompiledTranslation{1} returns a flat SVector{3} directly
        v = compiled_vector3(t_j2000)
        return SVector{3}(v[1], v[2], v[3])
    end
    return vector3(frames, ephem.center_point, ephem.target_point, ephem.axes, t_j2000)
end

"""
    get_velocity(ephem::FrameEphemeris, body::CelestialBody, frames, t_j2000)

Compute the position and velocity of a celestial body using the FrameSystem.

# Arguments
- `ephem::FrameEphemeris`: Ephemeris configuration (center, target, axes).
- `body::CelestialBody`: Celestial body (unused directly, kept for dispatch).
- `frames`: FrameSystem from FrameTransformations.jl.
- `t_j2000::Number`: Time in seconds since J2000 TDB.

# Returns
- `Tuple{SVector{3}, SVector{3}}`: (position [km], velocity [km/s]) in the specified frame.
"""
@inline function get_velocity(
    ephem::FrameEphemeris, body::CelestialBody, frames, t_j2000, compiled_vector6=nothing
)
    if !isnothing(compiled_vector6)
        # CompiledTranslation{2} returns a flat SVector{6} = [pos; vel]
        sv = compiled_vector6(t_j2000)
        return SVector{3}(sv[1], sv[2], sv[3]), SVector{3}(sv[4], sv[5], sv[6])
    end
    sv = vector6(frames, ephem.center_point, ephem.target_point, ephem.axes, t_j2000)
    return SVector{3}(sv[1], sv[2], sv[3]), SVector{3}(sv[4], sv[5], sv[6])
end
