# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Frame transformation utilities for AstroForceModels using FrameTransformations.jl
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export create_default_frames, create_frame_aware_params, add_small_body_rotating_frame!

"""
    create_default_frames(; order=2, numtype=Float64)

Create a default `FrameSystem` with standard celestial reference frames.

# Arguments
- `order::Int=2`: Order of frame transformations (1-4). Order 2 includes position and velocity derivatives.
- `numtype::Type=Float64`: Numeric type for computations.

# Returns
- `FrameSystem{order, numtype}`: Configured frame system.

# Frames Added
- `:ICRF` - International Celestial Reference Frame
- `:GCRF` - Geocentric Celestial Reference Frame  
- `:EME2000` - Earth Mean Equator 2000 / J2000

# Example
```julia
frames = create_default_frames()  # Basic inertial frames
frames_high_order = create_default_frames(order=4)  # Include up to jerk
```

# Notes
- For most astrodynamics applications, order=2 (position + velocity) is sufficient
- ICRF and GCRF are nearly identical for most practical purposes
- J2000/EME2000 is commonly used and approximately equal to ICRF
- For Earth-fixed frames (ITRF), use IERSConventions.jl add_axes_itrf! separately
"""
function create_default_frames(; order::Int=2, numtype::Type=Float64)
    frames = FrameSystem{order, numtype}()

    # Add standard inertial frames
    add_axes_icrf!(frames)
    add_axes_gcrf!(frames)
    add_axes_eme2000!(frames)

    return frames
end

"""
    create_frame_aware_params(base_params; frames=nothing, propagation_frame=:ICRF)

Create frame-aware parameters for force model propagation. The `Epoch` is automatically
derived from the `JD` field in `base_params`.

# Arguments
- `base_params::ComponentVector`: Base parameter vector (must include `JD` field).
- `frames::Union{FrameSystem,Nothing}=nothing`: Frame system (creates default if not provided).
- `propagation_frame::Symbol=:ICRF`: Frame in which state vectors are expressed.

# Returns
- `FrameAwareParams`: Wrapped parameter object with frame transformation capabilities.

# Example
```julia
using ComponentArrays

base = ComponentVector(JD = 2460000.5, μ = 398600.4415)
params = create_frame_aware_params(base)

# Access like normal ComponentVector
params.JD  # 2460000.5
params.μ   # 398600.4415

# Access frame info
params.frames
params.epoch
params.propagation_frame
```
"""
function create_frame_aware_params(
    base_params::ComponentVector;
    frames=nothing,
    propagation_frame::Symbol=:ICRF
)
    if frames === nothing
        frames = create_default_frames()
    end

    epoch = Epoch((base_params.JD - 2451545.0) * 86400.0, TDB)
    return FrameAwareParams(base_params, frames, epoch, propagation_frame)
end

"""
    add_small_body_rotating_frame!(
        frames::FrameSystem,
        name::Symbol,
        naif_id::Int,
        parent_frame_id::Int,
        rotation_axis::AbstractVector,
        rotation_period::Real;
        epoch_offset::Real=0.0
    )

Add a constantly rotating small body frame to the frame system.

# Arguments
- `frames::FrameSystem`: Frame system to add the frame to
- `name::Symbol`: Name for the new frame (e.g., `:ErosBodyFixed`)
- `naif_id::Int`: NAIF ID for the small body (e.g., 2000433 for Eros)
- `parent_frame_id::Int`: Parent frame ID (typically 1 for ICRF)
- `rotation_axis::AbstractVector`: Rotation axis unit vector in parent frame [x, y, z]
- `rotation_period::Real`: Rotation period in seconds
- `epoch_offset::Real=0.0`: Phase offset at J2000 epoch (radians)

# Returns
- Nothing (modifies `frames` in place)

# Example
```julia
frames = create_default_frames()

# Add Eros body-fixed frame (rotation period ≈ 5.27 hours)
add_small_body_rotating_frame!(
    frames,
    :ErosBodyFixed,
    2000433,  # NAIF ID for Eros
    1,        # Parent is ICRF
    [0.0, 0.0, 1.0],  # Rotation axis (example: aligned with Z)
    5.27 * 3600.0,    # Period in seconds
    epoch_offset = 0.0
)
```

# Notes
- This creates a simple constant rotation rate frame
- For more accurate small body orientations, use SPICE kernels via Ephemerides.jl
- The rotation axis should be normalized (unit vector)
- Positive rotation period indicates right-hand rotation about the axis
"""
function add_small_body_rotating_frame!(
    frames::FrameSystem{O,T},
    name::Symbol,
    naif_id::Int,
    parent_frame_id::Int,
    rotation_axis::AbstractVector,
    rotation_period::Real;
    epoch_offset::Real=0.0
) where {O,T}

    # Normalize rotation axis
    axis_norm = SVector{3,T}(rotation_axis ./ norm(rotation_axis))

    # Angular velocity (rad/s)
    ω = T(2π / rotation_period)
    ω_vec = ω .* axis_norm

    # Create rotation functions
    function rotation_function(ep)
        # ep is Tempo.Epoch, convert to seconds since J2000
        t = value(ep - Epoch(0.0, TDB))  # Seconds since J2000 TDB
        θ = ω * t + epoch_offset
        return angle_to_dcm(θ, axis_norm)
    end

    # Time derivative of DCM
    function drotation_function(ep)
        t = value(ep - Epoch(0.0, TDB))
        θ = ω * t + epoch_offset
        dcm = angle_to_dcm(θ, axis_norm)
        # ω̃ is the skew-symmetric matrix of angular velocity
        ω_skew = SMatrix{3,3,T}(
            0,        -ω_vec[3],  ω_vec[2],
            ω_vec[3],  0,        -ω_vec[1],
           -ω_vec[2],  ω_vec[1],  0
        )
        return ω_skew * dcm
    end

    # Second derivative (for order ≥ 3)
    function d2rotation_function(ep)
        if O < 3
            return zero(SMatrix{3,3,T})
        end
        t = value(ep - Epoch(0.0, TDB))
        θ = ω * t + epoch_offset
        dcm = angle_to_dcm(θ, axis_norm)
        ω_skew = SMatrix{3,3,T}(
            0,        -ω_vec[3],  ω_vec[2],
            ω_vec[3],  0,        -ω_vec[1],
           -ω_vec[2],  ω_vec[1],  0
        )
        return ω_skew * ω_skew * dcm
    end

    # Third derivative (for order = 4)
    function d3rotation_function(ep)
        if O < 4
            return zero(SMatrix{3,3,T})
        end
        t = value(ep - Epoch(0.0, TDB))
        θ = ω * t + epoch_offset
        dcm = angle_to_dcm(θ, axis_norm)
        ω_skew = SMatrix{3,3,T}(
            0,        -ω_vec[3],  ω_vec[2],
            ω_vec[3],  0,        -ω_vec[1],
           -ω_vec[2],  ω_vec[1],  0
        )
        return ω_skew * ω_skew * ω_skew * dcm
    end

    # Add the rotating frame
    add_axes_rotating!(
        frames,
        name,
        naif_id,
        parent_frame_id,
        rotation_function,
        drotation_function,
        d2rotation_function,
        d3rotation_function
    )

    return nothing
end

"""
    angle_to_dcm(angle::Real, axis::AbstractVector)

Convert rotation angle about an arbitrary axis to a Direction Cosine Matrix.

Uses Rodrigues' rotation formula:
    R = I + sin(θ)K + (1-cos(θ))K²
where K is the skew-symmetric matrix of the axis vector.
"""
function angle_to_dcm(angle::Real, axis::AbstractVector{T}) where T
    s = sin(angle)
    c = cos(angle)
    C = one(T) - c

    x, y, z = axis[1], axis[2], axis[3]

    return SMatrix{3,3,T}(
        c + x*x*C,     x*y*C - z*s,   x*z*C + y*s,
        y*x*C + z*s,   c + y*y*C,     y*z*C - x*s,
        z*x*C - y*s,   z*y*C + x*s,   c + z*z*C
    )
end
