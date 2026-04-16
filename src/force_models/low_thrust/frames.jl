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
export transform_thrust_to_state_frame

"""
    AbstractThrustFrame

Abstract base type for thrust reference frame specifications.

The frame determines how the 3-component acceleration vector produced by an
[`AbstractThrustModel`](@ref) is interpreted before being transformed into the
propagation (state) frame.

# Subtypes
- [`InertialFrame`](@ref): Components are in the inertial frame (e.g., ICRF)
- [`RTNFrame`](@ref): Radial / Transverse / Normal (orbit-fixed)
- [`VNBFrame`](@ref): Velocity / Normal / Binormal (velocity-fixed)
"""
abstract type AbstractThrustFrame end

"""
    InertialFrame{S} <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the inertial frame
named `S` (default `:ICRF`). When the propagation frame equals `S`, no rotation is
applied. When propagating in a different frame (e.g., a body-fixed rotating frame
or a different inertial frame), the `FrameSystem` on the `FrameAwareParams` is used
to rotate the thrust vector into the propagation frame at evaluation time.

The source frame `S` must be registered in `p.frames`.

# Constructors
```julia
InertialFrame()            # S = :ICRF (default)
InertialFrame(:J2000)      # or any other inertial frame registered in p.frames
```

This is the default frame for [`LowThrustAstroModel`](@ref).

!!! note
    The 3-argument `transform_thrust_to_state_frame` overload (used by impulsive
    maneuvers that do not carry a `FrameAwareParams`) is an identity. It is only
    correct when the propagation frame coincides with `S`. Use the 5-argument
    `FrameAwareParams` method when propagating in a different frame.
"""
struct InertialFrame{S} <: AbstractThrustFrame end

InertialFrame() = InertialFrame{:ICRF}()
InertialFrame(name::Symbol) = InertialFrame{name}()

@inline _inertial_frame_name(::InertialFrame{S}) where {S} = S

"""
    RTNFrame <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the RTN
(Radial–Transverse–Normal) frame, also known as RIC (Radial–In-track–Cross-track).

The RTN basis vectors are constructed from the spacecraft state in the propagation frame,
so this frame works correctly regardless of which reference frame is used for propagation.

- **R̂** (Radial): Along the position vector, away from the central body: `r / |r|`
- **N̂** (Normal): Along the orbital angular momentum: `(r × v) / |r × v|`
- **T̂** (Transverse): Completes the right-handed triad: `N̂ × R̂`

The acceleration components `[a_R, a_T, a_N]` are converted to the state frame via:

    a_state = a_R R̂ + a_T T̂ + a_N N̂
"""
struct RTNFrame <: AbstractThrustFrame end

"""
    VNBFrame <: AbstractThrustFrame

Indicates that the thrust acceleration vector is expressed in the VNB
(Velocity–Normal–Binormal) frame.

The VNB basis vectors are constructed from the spacecraft state in the propagation frame,
so this frame works correctly regardless of which reference frame is used for propagation.

- **V̂** (Velocity): Along the velocity vector: `v / |v|`
- **N̂** (Normal): In the orbital plane, perpendicular to V̂: `B̂ × V̂`
- **B̂** (Binormal): Along the orbital angular momentum: `(r × v) / |r × v|`

The acceleration components `[a_V, a_N, a_B]` are converted to the state frame via:

    a_state = a_V V̂ + a_N N̂ + a_B B̂

!!! note
    For circular orbits the VNB frame coincides with the RTN frame (V̂ ≈ T̂).
    The distinction becomes important for eccentric orbits where the velocity
    vector deviates from the transverse direction.
"""
struct VNBFrame <: AbstractThrustFrame end

# ── RTN / VNB: orbital frames derived from state ─────────────────────────────
# These are correct in any propagation frame because the basis vectors are
# constructed from the state vector, which is already in the propagation frame.

"""
    transform_thrust_to_state_frame(a_local::SVector{3}, u, p, t, ::RTNFrame)

Transform a thrust acceleration from the RTN frame to the propagation (state) frame.

The RTN basis is constructed from the spacecraft state `u = [r; v]`:

    R̂ = r / |r|,  N̂ = (r × v) / |r × v|,  T̂ = N̂ × R̂

Returns `a_R R̂ + a_T T̂ + a_N N̂` in the propagation frame.
"""
@inline function transform_thrust_to_state_frame(
    a_local::SVector{3,AT}, u::AbstractVector{UT}, p, t, ::RTNFrame
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
    transform_thrust_to_state_frame(a_local::SVector{3}, u, p, t, ::VNBFrame)

Transform a thrust acceleration from the VNB frame to the propagation (state) frame.

The VNB basis is constructed from the spacecraft state `u = [r; v]`:

    V̂ = v / |v|,  B̂ = (r × v) / |r × v|,  N̂ = B̂ × V̂

Returns `a_V V̂ + a_N N̂ + a_B B̂` in the propagation frame.
"""
@inline function transform_thrust_to_state_frame(
    a_local::SVector{3,AT}, u::AbstractVector{UT}, p, t, ::VNBFrame
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

# ── InertialFrame: rotate from source inertial frame to propagation frame ──

"""
    transform_thrust_to_state_frame(
        a_inertial::SVector{3}, u, p::FrameAwareParams, t, frame::InertialFrame{S}
    ) where {S}

Transform a thrust acceleration from the inertial frame `S` (default `:ICRF`) to
`p.propagation_frame`.

When `S === p.propagation_frame`, this short-circuits to the identity and adds no
overhead. Otherwise a single `rotation3` lookup against `p.frames` produces the
DCM `R_{S → p.propagation_frame}` and the result is `R * a_inertial`.
"""
@inline function transform_thrust_to_state_frame(
    a_inertial::SVector{3},
    u::AbstractVector,
    p::FrameAwareParams,
    t,
    frame::InertialFrame{S},
) where {S}
    if S === p.propagation_frame
        return a_inertial
    else
        t_ft = ft_time(p, t)
        R = rotation3(p.frames, S, p.propagation_frame, t_ft)
        return R.m[1] * a_inertial
    end
end

# ── 3-arg convenience for impulsive maneuvers (no p/t needed) ────────────────
# These are used by AstroPropagators impulsive_burn! where the ΔV is applied
# directly to the state and only the orbital-frame rotation matters.

"""
    transform_thrust_to_state_frame(
        a_local::SVector{3}, u::AbstractVector, frame::AbstractThrustFrame
    )

3-argument convenience method for transforming thrust/ΔV vectors when frame-aware
parameters are not available (e.g., impulsive maneuvers applied directly to a state
vector).

- `RTNFrame` / `VNBFrame`: always correct — the orthonormal basis is constructed
  from `u = [r; v]` in the propagation frame, so the result is automatically in
  the propagation frame regardless of which frame that is.
- `InertialFrame{S}`: identity operation. This is **only** correct when the
  propagation frame coincides with `S`. Use the 5-argument method with
  `FrameAwareParams` whenever the propagation frame may differ from `S`.
"""
@inline function transform_thrust_to_state_frame(
    a_local::SVector{3}, u::AbstractVector, frame::RTNFrame
)
    return transform_thrust_to_state_frame(a_local, u, nothing, nothing, frame)
end

@inline function transform_thrust_to_state_frame(
    a_local::SVector{3}, u::AbstractVector, frame::VNBFrame
)
    return transform_thrust_to_state_frame(a_local, u, nothing, nothing, frame)
end

@inline function transform_thrust_to_state_frame(
    a_local::SVector{3}, u::AbstractVector, ::InertialFrame
)
    return a_local
end
