# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Reference Frame Definitions for Low-Thrust Acceleration
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications (4th ed.).
#       Microcosm Press, Chapter 3.
#   [2] Schaub, H. and Junkins, J. L. (2018). Analytical Mechanics of Space Systems
#       (4th ed.). AIAA Education Series, Chapter 1.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AbstractThrustFrame, InertialFrame, RTNFrame, VNBFrame
export transform_to_inertial

"""
    AbstractThrustFrame

Abstract base type for thrust reference frame specifications.

The frame determines how the 3-component acceleration vector produced by an
[`AbstractThrustModel`](@ref) is interpreted before being transformed into the
inertial (ECI) frame.

# Subtypes
- [`InertialFrame`](@ref): Components are inertial (ECI) — no rotation needed
- [`RTNFrame`](@ref): Radial / Transverse / Normal (orbit-fixed)
- [`VNBFrame`](@ref): Velocity / Normal / Binormal (velocity-fixed)
"""
abstract type AbstractThrustFrame end

"""
    InertialFrame <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the inertial (ECI/J2000)
frame. No rotation is applied.

This is the default frame for [`LowThrustAstroModel`](@ref).
"""
struct InertialFrame <: AbstractThrustFrame end

"""
    RTNFrame <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the RTN
(Radial–Transverse–Normal) frame, also known as RIC (Radial–In-track–Cross-track).

The RTN basis vectors are defined as:

- **R̂** (Radial): Along the position vector, away from the central body: `r / |r|`
- **N̂** (Normal): Along the orbital angular momentum: `(r × v) / |r × v|`
- **T̂** (Transverse): Completes the right-handed triad: `N̂ × R̂`

The acceleration components `[a_R, a_T, a_N]` are converted to inertial via:

    a_inertial = a_R R̂ + a_T T̂ + a_N N̂
"""
struct RTNFrame <: AbstractThrustFrame end

"""
    VNBFrame <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the VNB
(Velocity–Normal–Binormal) frame.

The VNB basis vectors are defined as:

- **V̂** (Velocity): Along the velocity vector: `v / |v|`
- **N̂** (Normal): In the orbital plane, perpendicular to V̂: `B̂ × V̂`
- **B̂** (Binormal): Along the orbital angular momentum: `(r × v) / |r × v|`

The acceleration components `[a_V, a_N, a_B]` are converted to inertial via:

    a_inertial = a_V V̂ + a_N N̂ + a_B B̂

!!! note
    For circular orbits the VNB frame coincides with the RTN frame (V̂ ≈ T̂).
    The distinction becomes important for eccentric orbits where the velocity
    vector deviates from the transverse direction.
"""
struct VNBFrame <: AbstractThrustFrame end

"""
    transform_to_inertial(a_local::SVector{3}, u::AbstractVector, frame::InertialFrame)

Identity transformation — the acceleration is already in the inertial frame.
"""
@inline function transform_to_inertial(
    a_local::SVector{3}, u::AbstractVector, ::InertialFrame
)
    return a_local
end

"""
    transform_to_inertial(a_local::SVector{3}, u::AbstractVector, frame::RTNFrame)

Transform a thrust acceleration from the RTN frame to the inertial frame.

The RTN basis is constructed from the spacecraft state `u = [r; v]` as:

    R̂ = r / |r|
    N̂ = (r × v) / |r × v|
    T̂ = N̂ × R̂

The inertial acceleration is: `a_R R̂ + a_T T̂ + a_N N̂`.
"""
@inline function transform_to_inertial(
    a_local::SVector{3,AT}, u::AbstractVector{UT}, ::RTNFrame
) where {AT,UT}
    RT = promote_type(AT, UT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    v = SVector{3,UT}(u[4], u[5], u[6])

    r_norm = norm(r)
    R̂ = SVector{3}(r[1] / r_norm, r[2] / r_norm, r[3] / r_norm)

    h = cross(r, v)
    h_norm = norm(h)
    N̂ = SVector{3}(h[1] / h_norm, h[2] / h_norm, h[3] / h_norm)

    T̂ = cross(N̂, R̂)

    return SVector{3,RT}(
        a_local[1] * R̂[1] + a_local[2] * T̂[1] + a_local[3] * N̂[1],
        a_local[1] * R̂[2] + a_local[2] * T̂[2] + a_local[3] * N̂[2],
        a_local[1] * R̂[3] + a_local[2] * T̂[3] + a_local[3] * N̂[3],
    )
end

"""
    transform_to_inertial(a_local::SVector{3}, u::AbstractVector, frame::VNBFrame)

Transform a thrust acceleration from the VNB frame to the inertial frame.

The VNB basis is constructed from the spacecraft state `u = [r; v]` as:

    V̂ = v / |v|
    B̂ = (r × v) / |r × v|
    N̂ = B̂ × V̂

The inertial acceleration is: `a_V V̂ + a_N N̂ + a_B B̂`.
"""
@inline function transform_to_inertial(
    a_local::SVector{3,AT}, u::AbstractVector{UT}, ::VNBFrame
) where {AT,UT}
    RT = promote_type(AT, UT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    v = SVector{3,UT}(u[4], u[5], u[6])

    v_norm = norm(v)
    V̂ = SVector{3}(v[1] / v_norm, v[2] / v_norm, v[3] / v_norm)

    h = cross(r, v)
    h_norm = norm(h)
    B̂ = SVector{3}(h[1] / h_norm, h[2] / h_norm, h[3] / h_norm)

    N̂ = cross(B̂, V̂)

    return SVector{3,RT}(
        a_local[1] * V̂[1] + a_local[2] * N̂[1] + a_local[3] * B̂[1],
        a_local[1] * V̂[2] + a_local[2] * N̂[2] + a_local[3] * B̂[2],
        a_local[1] * V̂[3] + a_local[2] * N̂[3] + a_local[3] * B̂[3],
    )
end
