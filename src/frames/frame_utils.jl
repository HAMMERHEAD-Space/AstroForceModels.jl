# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Frame transformation utilities for AstroForceModels using FrameTransformations.jl
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export create_default_frames, create_frame_aware_params, add_small_body_rotating_frame!
export setup_inertial_frames, setup_earth_propagation_frames
export setup_ephemeris_frames

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
    frames = FrameSystem{order,numtype}()

    # Add standard inertial frames
    add_axes_icrf!(frames)
    add_axes_gcrf!(frames)
    add_axes_eme2000!(frames)

    return frames
end

"""
    create_frame_aware_params(epoch, frames; propagation_frame=:ICRF, extra_params...)

Create frame-aware parameters for force model propagation.

# Arguments
- `epoch::Epoch`: Reference epoch. JD is derived from this automatically.
- `frames::FrameSystem`: Frame system for coordinate transformations.
- `propagation_frame::Symbol=:ICRF`: Frame in which state vectors are expressed.
- `extra_params...`: Any additional named parameters to include.

# Example
```julia
p = create_frame_aware_params(epoch, frames)
p.JD   # derived from epoch
```
"""
function create_frame_aware_params(
    epoch::Epoch, frames; propagation_frame::Symbol=:ICRF, extra_params...
)
    params = isempty(extra_params) ? ComponentVector() : ComponentVector(; extra_params...)
    return FrameAwareParams(params, frames, epoch, propagation_frame)
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
    epoch_offset::Real=0.0,
) where {O,T}

    # Normalize rotation axis
    axis_norm = SVector{3,T}(rotation_axis ./ norm(rotation_axis))

    # Angular velocity (rad/s)
    ω = T(2π / rotation_period)
    ω_vec = ω .* axis_norm

    # Create rotation functions
    # FrameTransformations passes t as Float64 seconds since J2000 TDB
    function rotation_function(t)
        θ = ω * t + epoch_offset
        return angleaxis_to_dcm(θ, axis_norm)
    end

    # Time derivative of DCM
    function drotation_function(t)
        θ = ω * t + epoch_offset
        dcm = angleaxis_to_dcm(θ, axis_norm)
        ω_skew = SMatrix{3,3,T}(
            0, -ω_vec[3], ω_vec[2], ω_vec[3], 0, -ω_vec[1], -ω_vec[2], ω_vec[1], 0
        )
        return ω_skew * dcm
    end

    # Second derivative (for order ≥ 3)
    function d2rotation_function(t)
        if O < 3
            return zero(SMatrix{3,3,T})
        end
        θ = ω * t + epoch_offset
        dcm = angleaxis_to_dcm(θ, axis_norm)
        ω_skew = SMatrix{3,3,T}(
            0, -ω_vec[3], ω_vec[2], ω_vec[3], 0, -ω_vec[1], -ω_vec[2], ω_vec[1], 0
        )
        return ω_skew * ω_skew * dcm
    end

    # Third derivative (for order = 4)
    function d3rotation_function(t)
        if O < 4
            return zero(SMatrix{3,3,T})
        end
        θ = ω * t + epoch_offset
        dcm = angleaxis_to_dcm(θ, axis_norm)
        ω_skew = SMatrix{3,3,T}(
            0, -ω_vec[3], ω_vec[2], ω_vec[3], 0, -ω_vec[1], -ω_vec[2], ω_vec[1], 0
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
        d3rotation_function,
    )

    return nothing
end

"""
    setup_inertial_frames(epoch; ephemeris=Vallado(), include_sun=true, include_moon=true, order=2, numtype=Float64)

One-call setup for Earth-centric inertial propagation. Returns a `FrameAwareParams`
ready to pass directly to `build_dynamics_model` or an ODE integrator.

Builds a `FrameSystem` with:
- `:ICRF` — root inertial axes
- `:Earth` (NAIF 399) — root point at the origin of `:ICRF`
- Body points registered via the selected ephemeris backend

# Ephemeris Backends
- `Vallado()` (default): Sun + Moon only, moderate accuracy. Requires `SatelliteToolboxCelestialBodies`.
- `Meeus()`: Sun + all 8 planets, ~1000 km accuracy.
- `Kepler()`: All bodies with J2000 Keplerian elements, low accuracy.

# Arguments
- `epoch::Epoch`: Propagation start epoch.
- `ephemeris::AbstractEphemerisType=Vallado()`: Ephemeris backend for body positions.
- `include_sun::Bool=true`: Add the Sun as a dynamical point.
- `include_moon::Bool=true`: Add the Moon as a dynamical point.
- `order::Int=2`: Derivative order of the frame system (2 = position + velocity).
- `numtype::Type=Float64`: Numeric type for the frame system.

# Returns
- `FrameAwareParams` with `propagation_frame = :ICRF`.

# Example
```julia
using AstroForceModels, Tempo

epoch = Epoch("2024-01-05T12:00:00 TDB")

# Default Vallado: Sun + Moon
p = setup_inertial_frames(epoch)

# Meeus: Sun + all planets
p = setup_inertial_frames(epoch; ephemeris=Meeus())

# Kepler: all bodies from J2000 elements
p = setup_inertial_frames(epoch; ephemeris=Kepler())

# Minimal: ICRF + Earth only
p = setup_inertial_frames(epoch; include_sun=false, include_moon=false)
```
"""
function setup_inertial_frames(
    epoch::Epoch;
    ephemeris::CelestialBodies.AbstractEphemerisType=Vallado(),
    include_sun::Bool=true,
    include_moon::Bool=true,
    order::Int=2,
    numtype::Type=Float64,
)
    frames = FrameSystem{order,numtype}()
    add_axes_icrf!(frames)

    if ephemeris isa Vallado
        # Earth-centric: Earth as root, Sun + Moon relative to Earth
        add_point!(frames, :Earth, 399, :ICRF)
        include_sun &&
            add_body_point!(frames, SunBody(), ephemeris; parent_id=399, axes=:ICRF)
        include_moon &&
            add_body_point!(frames, MoonBody(), ephemeris; parent_id=399, axes=:ICRF)
    elseif ephemeris isa Meeus
        # Heliocentric: Sun as root, planets relative to Sun, Moon relative to Earth
        add_point!(frames, :Sun, 10, :ICRF)
        add_body_point!(frames, EarthBody(), ephemeris; parent_id=10, axes=:ICRF)
        for body_fn in [
            MercuryBody,
            VenusBody,
            MarsBody,
            JupiterBody,
            SaturnBody,
            UranusBody,
            NeptuneBody,
        ]
            add_body_point!(frames, body_fn(), ephemeris; parent_id=10, axes=:ICRF)
        end
        include_moon &&
            add_body_point!(frames, MoonBody(), Vallado(); parent_id=399, axes=:ICRF)
    elseif ephemeris isa Kepler
        # Heliocentric: Sun as root, all bodies from J2000 Keplerian elements
        add_point!(frames, :Sun, 10, :ICRF)
        for body_fn in [
            MercuryKeplerianBody,
            VenusKeplerianBody,
            EarthKeplerianBody,
            MarsKeplerianBody,
            JupiterKeplerianBody,
            SaturnKeplerianBody,
            UranusKeplerianBody,
            NeptuneKeplerianBody,
            PlutoKeplerianBody,
        ]
            add_body_point!(frames, body_fn(), ephemeris; parent_id=10, axes=:ICRF)
        end
        include_moon && add_body_point!(
            frames, MoonKeplerianBody(), ephemeris; parent_id=399, axes=:ICRF
        )
    end

    return FrameAwareParams(frames, epoch, :ICRF)
end

"""
    setup_earth_propagation_frames(epoch, eop_data; include_sun=true, include_moon=true, order=2, numtype=Float64)

One-call setup for Earth-centric propagation with body-fixed (ITRF) frame. Returns a
`FrameAwareParams` ready to pass directly to `build_dynamics_model` or an ODE integrator.

Builds a `FrameSystem` with:
- `:ICRF` — root inertial axes
- `:ITRF` — Earth body-fixed rotating axes (via SatelliteToolboxTransformations EOP)
- `:Earth` (NAIF 399) — root point at the origin of `:ICRF`
- `:Sun` (NAIF 10) — Vallado analytical ephemeris, if `include_sun=true`
- `:Moon` (NAIF 301) — Vallado analytical ephemeris, if `include_moon=true`

This is the recommended entry point for LEO/MEO/GEO propagation that uses gravity
harmonics, atmospheric drag, albedo, or any other force model requiring the ITRF frame.

# Arguments
- `epoch::Epoch`: Propagation start epoch.
- `eop_data`: Earth Orientation Parameters from `SatelliteToolboxTransformations.fetch_iers_eop()`.
- `include_sun::Bool=true`: Add the Sun as a dynamical point.
- `include_moon::Bool=true`: Add the Moon as a dynamical point.
- `order::Int=2`: Derivative order of the frame system (2 = position + velocity).
- `numtype::Type=Float64`: Numeric type for the frame system.

# Returns
- `FrameAwareParams` with `propagation_frame = :ICRF`.

# Example
```julia
using AstroForceModels, SatelliteToolboxTransformations, SatelliteToolboxGravityModels, Tempo

epoch = Epoch("2024-01-05T12:00:00 TDB")
eop_data = fetch_iers_eop()
p = setup_earth_propagation_frames(epoch, eop_data)

# Everything needed for high-fidelity Earth propagation is ready
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
gravity = GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, body_fixed_frame=:ITRF, propagation_frame=:ICRF,
    order=36, degree=36, frames=p.frames,
)
sun  = ThirdBodyModel(; body=SunBody(),  ephem_type=FrameEphemeris(; center_point=399, target_point=10,  axes=:ICRF), frames=p.frames)
moon = ThirdBodyModel(; body=MoonBody(), ephem_type=FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF), frames=p.frames)
drag = DragAstroModel(; satellite_drag_model=CannonballFixedDrag(0.2), atmosphere_model=JB2008(), frames=p.frames)

dynamics = CentralBodyDynamicsModel(gravity, (sun, moon, drag))
```

# See Also
- [`setup_inertial_frames`](@ref): Inertial-only setup (no ITRF).
"""
function setup_earth_propagation_frames(
    epoch::Epoch,
    eop_data;
    include_sun::Bool=true,
    include_moon::Bool=true,
    order::Int=2,
    numtype::Type=Float64,
)
    frames = FrameSystem{order,numtype}()
    add_axes_icrf!(frames)

    # Add ITRF as rotating frame using SatelliteToolboxTransformations EOP data
    add_axes_rotating!(
        frames,
        :ITRF,
        2,
        1,
        t -> SatelliteToolboxTransformations.r_eci_to_ecef(
            SatelliteToolboxTransformations.DCM,
            SatelliteToolboxTransformations.J2000(),
            SatelliteToolboxTransformations.ITRF(),
            2451545.0 + t / 86400.0,
            eop_data,
        ),
        t -> begin
            jd = 2451545.0 + t / 86400.0
            dt = 0.01
            R1 = SatelliteToolboxTransformations.r_eci_to_ecef(
                SatelliteToolboxTransformations.DCM,
                SatelliteToolboxTransformations.J2000(),
                SatelliteToolboxTransformations.ITRF(),
                jd,
                eop_data,
            )
            R2 = SatelliteToolboxTransformations.r_eci_to_ecef(
                SatelliteToolboxTransformations.DCM,
                SatelliteToolboxTransformations.J2000(),
                SatelliteToolboxTransformations.ITRF(),
                jd + dt / 86400.0,
                eop_data,
            )
            return (R2 - R1) / dt
        end,
    )

    add_point!(frames, :Earth, 399, :ICRF)

    include_sun && add_body_point!(frames, SunBody(), Vallado(); parent_id=399, axes=:ICRF)
    include_moon &&
        add_body_point!(frames, MoonBody(), Vallado(); parent_id=399, axes=:ICRF)

    return FrameAwareParams(frames, epoch, :ICRF)
end

"""
    setup_ephemeris_frames(epoch, eph; kwargs...)

Build a `FrameSystem` and register all points and orientation axes from SPK/PCK kernels
loaded in an `EphemerisProvider`. Returns a `FrameAwareParams` ready for propagation.

This is the highest-fidelity setup: body positions come from JPL binary SPK kernels
(e.g., DE440) and body orientations from binary PCK kernels (e.g., `earth_latest_high_prec.bpc`,
`moon_pa_de440.bpc`).

!!! note
    Requires `Ephemerides.jl` to be loaded (`using Ephemerides`) since FrameTransformations
    uses it via a package extension.

# Arguments
- `epoch::Epoch`: Propagation start epoch.
- `eph`: Loaded SPK/PCK kernel data (from `Ephemerides.EphemerisProvider`).
  Can include both `.bsp` (position) and `.bpc` (orientation) files:
  `EphemerisProvider(["de440.bsp", "earth_latest_high_prec.bpc"])`.
- `eop_data=nothing`: If provided, adds ITRF via SatelliteToolboxTransformations EOP
  (same as `setup_earth_propagation_frames`). Takes priority over PCK-based Earth
  orientation and `include_itrf`.
- `book::Union{Dict{Int,Symbol},Nothing}=nothing`: Mapping from NAIF IDs to point names.
  If `nothing`, a default book for common solar system bodies is used.
- `axes_book::Union{Dict{Int,Tuple{Symbol,Symbol}},Nothing}=nothing`: Mapping from NAIF
  axes IDs to `(name, rotation_sequence)` pairs for PCK orientation data. If `nothing`,
  a default book for common body-fixed frames is used. Set to an empty `Dict` to skip
  PCK axes registration entirely.
- `include_itrf::Bool=false`: Add ITRF using IERSConventions (via
  `FrameTransformations.add_axes_itrf!`). Ignored if `eop_data` is provided.
- `propagation_frame::Symbol=:ICRF`: Frame in which state vectors are expressed.
- `order::Int=2`: Derivative order of the frame system.
- `numtype::Type=Float64`: Numeric type for the frame system.

# Returns
- `FrameAwareParams` with the specified `propagation_frame`.

# Default Point Book (NAIF IDs → Names)
When `book=nothing`:
- `0 => :SSB`, `1–9 => :*Barycenter`, `10 => :Sun`
- `199 => :Mercury`, `299 => :Venus`, `301 => :Moon`, `399 => :Earth`
- `499 => :Mars`, `599 => :Jupiter`, `699 => :Saturn`
- `799 => :Uranus`, `899 => :Neptune`, `999 => :Pluto`

# Default Axes Book (PCK Axes IDs → Names + Rotation Sequences)
When `axes_book=nothing`:
- `3000 => (:EarthPCK, :ZXZ)` — Earth body-fixed from high-precision PCK
- `31006 => (:MoonME, :ZXZ)` — Moon Mean Earth frame
- `31007 => (:MoonPA421, :ZXZ)` — Moon Principal Axes DE421
- `31008 => (:MoonPA440, :ZXZ)` — Moon Principal Axes DE440

Only axes IDs actually present in the loaded kernels are registered.

# Example
```julia
using AstroForceModels, Ephemerides, Tempo

# SPK only — planetary positions from DE440
eph = EphemerisProvider("de440s.bsp")
epoch = Epoch("2024-01-05T12:00:00 TDB")
p = setup_ephemeris_frames(epoch, eph)

# SPK + PCK — positions and Earth orientation from kernels
eph = EphemerisProvider(["de440s.bsp", "earth_latest_high_prec.bpc"])
p = setup_ephemeris_frames(epoch, eph)
# p.frames now has :EarthPCK axes from the PCK kernel

# SPK + PCK + ITRF from SatelliteToolbox EOP (overrides PCK Earth rotation)
using SatelliteToolboxTransformations
eop_data = fetch_iers_eop()
p = setup_ephemeris_frames(epoch, eph; eop_data=eop_data)

# SPK + Moon PA kernel for lunar missions
eph = EphemerisProvider(["de440s.bsp", "moon_pa_de440.bpc"])
p = setup_ephemeris_frames(epoch, eph)
# p.frames now has :MoonPA440 axes
```

# See Also
- [`setup_inertial_frames`](@ref): Minimal inertial setup with Vallado analytical ephemeris.
- [`setup_earth_propagation_frames`](@ref): Earth-centric with ITRF + Vallado ephemeris.
"""
function setup_ephemeris_frames(
    epoch::Epoch,
    eph;
    eop_data=nothing,
    book::Union{Dict{Int,Symbol},Nothing}=nothing,
    axes_book::Union{Dict{Int,Tuple{Symbol,Symbol}},Nothing}=nothing,
    include_itrf::Bool=false,
    propagation_frame::Symbol=:ICRF,
    order::Int=2,
    numtype::Type=Float64,
)
    frames = FrameSystem{order,numtype}()

    # Standard inertial axes — always needed as the root of the graph
    add_axes_icrf!(frames)
    add_axes_gcrf!(frames)
    add_axes_eme2000!(frames)
    add_axes_ecl2000!(frames)   # Ecliptic J2000 (NAIF 17) — required parent for many PCK files

    # Add ITRF if requested via EOP data
    if !isnothing(eop_data)
        add_axes_rotating!(
            frames,
            :ITRF,
            2,
            1,
            t -> SatelliteToolboxTransformations.r_eci_to_ecef(
                SatelliteToolboxTransformations.DCM,
                SatelliteToolboxTransformations.J2000(),
                SatelliteToolboxTransformations.ITRF(),
                2451545.0 + t / 86400.0,
                eop_data,
            ),
            t -> begin
                jd = 2451545.0 + t / 86400.0
                dt = 0.01
                R1 = SatelliteToolboxTransformations.r_eci_to_ecef(
                    SatelliteToolboxTransformations.DCM,
                    SatelliteToolboxTransformations.J2000(),
                    SatelliteToolboxTransformations.ITRF(),
                    jd,
                    eop_data,
                )
                R2 = SatelliteToolboxTransformations.r_eci_to_ecef(
                    SatelliteToolboxTransformations.DCM,
                    SatelliteToolboxTransformations.J2000(),
                    SatelliteToolboxTransformations.ITRF(),
                    jd + dt / 86400.0,
                    eop_data,
                )
                return (R2 - R1) / dt
            end,
        )
    elseif include_itrf
        add_axes_itrf!(frames, :ITRF, 1, 2)
    end

    # Register all PCK orientation axes found in the kernel
    available_axes = try
        ephem_get_axes(eph)
    catch
        Int[]
    end

    if !isempty(available_axes)
        # Build the axes name book — use provided or generate from NAIF IDs
        if isnothing(axes_book)
            axes_book = _default_axes_book(available_axes)
        end

        for axes_id in available_axes
            haskey(axes_book, axes_id) || continue
            name, rot_seq = axes_book[axes_id]
            has_axes(frames, axes_id) && continue  # Skip if already registered
            try
                add_axes_ephemeris!(frames, eph, name, axes_id, rot_seq)
            catch e
                @warn "Could not register PCK axes $name (ID=$axes_id)" exception = e
            end
        end
    end

    # Register all SPK points found in the kernel
    available_points = try
        ephem_get_points(eph)
    catch
        Int[]
    end

    if !isempty(available_points)
        if isnothing(book)
            book = _default_point_book(available_points)
        end
        add_point_ephemeris!(frames, eph, book)
    end

    return FrameAwareParams(frames, epoch, propagation_frame)
end

# ── Default NAIF ID name generators ─────────────────────────────────────────

const _NAIF_POINT_NAMES = Dict{Int,Symbol}(
    0 => :SSB,
    1 => :MercuryBarycenter,
    2 => :VenusBarycenter,
    3 => :EarthBarycenter,
    4 => :MarsBarycenter,
    5 => :JupiterBarycenter,
    6 => :SaturnBarycenter,
    7 => :UranusBarycenter,
    8 => :NeptuneBarycenter,
    9 => :PlutoBarycenter,
    10 => :Sun,
    199 => :Mercury,
    299 => :Venus,
    301 => :Moon,
    399 => :Earth,
    499 => :Mars,
    599 => :Jupiter,
    699 => :Saturn,
    799 => :Uranus,
    899 => :Neptune,
    999 => :Pluto,
)

const _NAIF_AXES_NAMES = Dict{Int,Tuple{Symbol,Symbol}}(
    3000 => (:EarthPCK, :ZXZ),
    31006 => (:MoonME, :ZXZ),
    31007 => (:MoonPA421, :ZXZ),
    31008 => (:MoonPA440, :ZXZ),
)

"""
    _default_point_book(available_ids) -> Dict{Int, Symbol}

Generate a point name book for the given NAIF IDs. Uses known names for standard
solar system bodies and generates `:Body_NNNNN` names for unknown IDs.
"""
function _default_point_book(available_ids)
    book = Dict{Int,Symbol}()
    for id in available_ids
        if haskey(_NAIF_POINT_NAMES, id)
            book[id] = _NAIF_POINT_NAMES[id]
        else
            book[id] = Symbol("Body_", id)
        end
    end
    return book
end

"""
    _default_axes_book(available_ids) -> Dict{Int, Tuple{Symbol, Symbol}}

Generate an axes name book for the given NAIF axes IDs. Uses known names for standard
body-fixed frames and generates `:Axes_NNNNN` names with `:ZXZ` rotation for unknown IDs.
"""
function _default_axes_book(available_ids)
    book = Dict{Int,Tuple{Symbol,Symbol}}()
    for id in available_ids
        if haskey(_NAIF_AXES_NAMES, id)
            book[id] = _NAIF_AXES_NAMES[id]
        else
            book[id] = (Symbol("Axes_", id), :ZXZ)
        end
    end
    return book
end

