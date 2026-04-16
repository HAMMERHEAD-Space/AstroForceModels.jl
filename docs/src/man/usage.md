# Usage

AstroForceModels.jl provides a comprehensive framework for modeling astrodynamics forces affecting orbital motion around any central body. The package uses FrameTransformations.jl for all coordinate transformations and ephemeris queries, making it fully frame-independent.

## Quick Start

`setup_inertial_frames` sets up a complete Earth-centric inertial frame system in one call and returns a `FrameAwareParams` ready to pass straight into `build_dynamics_model`:

```julia
using AstroForceModels
using Tempo

epoch = Epoch("2024-01-05T12:00:00 TDB")

# ICRF + Earth + Sun (Vallado) + Moon (Vallado) — enough for most LEO propagations
p = setup_inertial_frames(epoch)

# Build force models using the frame system from p
grav = KeplerianGravityAstroModel(; μ = 3.986004415e5)
sun  = ThirdBodyModel(; body=SunBody(),  ephem_type=FrameEphemeris(; center_point=399, target_point=10,  axes=:ICRF), frames=p.frames)
moon = ThirdBodyModel(; body=MoonBody(), ephem_type=FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF), frames=p.frames)

dynamics = CentralBodyDynamicsModel(grav, (sun, moon))

function satellite_ode!(du, u, p, t)
    du[1:3] = u[4:6]
    du[4:6] = build_dynamics_model(u, p, t, dynamics)
end
```

If you do not need third-body perturbations or SRP, omit the Sun and Moon to keep the frame system minimal:

```julia
p = setup_inertial_frames(epoch; include_sun=false, include_moon=false)
```

Keyword options:

| Keyword | Default | Description |
|---------|---------|-------------|
| `include_sun` | `true` | Add Sun (NAIF 10) via Vallado ephemeris |
| `include_moon` | `true` | Add Moon (NAIF 301) via Vallado ephemeris |
| `order` | `2` | Derivative order of the frame system |
| `numtype` | `Float64` | Numeric type for the frame system |

For Earth propagations that need ITRF (gravity harmonics, drag), use `setup_earth_propagation_frames` instead. For JPL SPK kernel-based ephemeris, use `setup_ephemeris_frames`. See [Ephemeris Sources](@ref) for a detailed comparison of all three tiers.

## Frame-Aware Parameters

All force models accept `FrameAwareParams` as the parameter type. This stores a `FrameSystem`, epoch, propagation frame, and an optional `ComponentVector` of user-defined parameters. JD is derived automatically from the epoch; μ is read from the dynamics model — neither belongs in `FrameAwareParams`:

```julia
using AstroForceModels
using FrameTransformations
using Tempo

# Create a FrameSystem with standard frames
frames = create_default_frames()

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)

# Minimal construction (no user params needed)
p = FrameAwareParams(frames, epoch, :ICRF)

# Access
p.JD          # Julian date (derived from epoch)
p.frames      # FrameSystem
p.epoch       # Tempo.Epoch
```

## Force Model Interface

Each force model implements the common `acceleration` interface:

```julia
acceleration(state, parameters, time, force_model)
```

Where:
- `state`: Current spacecraft state vector (position and velocity) [km, km/s]
- `parameters`: `FrameAwareParams` containing frame system, epoch, and numeric params
- `time`: Current simulation time (seconds since epoch)
- `force_model`: Specific force model instance

## Third Body Models (FrameEphemeris)

Third-body positions are obtained from the `FrameSystem` via `FrameEphemeris`, which uses NAIF IDs:

```julia
# Sun position relative to Earth in ICRF
sun_model = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
)

# Moon position relative to Earth
moon_model = ThirdBodyModel(
    body = MoonBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
)
```

## Compiled Transforms (Allocation-Free Hot Path)

When you pass `frames` to a model constructor, the rotation or translation closures are
extracted from the `FrameSystem` at construction time via `compile_rotation` /
`compile_translation` and stored in type-parameterized structs. This bypasses
`FunctionWrapper` dispatch and graph lookups on every call, giving allocation-free
evaluation that passes `@check_allocs`:

```julia
# Without compiled transforms (allocates on every call via FunctionWrapper)
grav = GravityHarmonicsAstroModel(;
    gravity_model=coeffs, body_fixed_frame=:ITRF, propagation_frame=:ICRF,
)

# With compiled transforms (allocation-free)
grav = GravityHarmonicsAstroModel(;
    gravity_model=coeffs, body_fixed_frame=:ITRF, propagation_frame=:ICRF,
    frames=my_frames,   # pre-compiles ICRF↔ITRF rotation
)
```

The compiled callables handle both direct parent-child pairs and multi-hop paths
through the frame graph, so any connected pair of points or axes works (e.g.
`ICRF→ITRF`, `ICRF→IAU_MARS`, `ICRF→ErosBodyFixed`, `SSB→Earth` via `EMB`, etc.).
Construction is guarded with a `try`/`catch`; if compilation fails for any reason
the model transparently falls back to runtime `rotation3`/`vector3`/`vector6` lookups.

## Combining Force Models

Multiple force models can be combined using `CentralBodyDynamicsModel`:

```julia
dynamics_model = CentralBodyDynamicsModel(
    gravity_model,
    (drag_model, srp_model, sun_model, moon_model)
)

# System dynamics function for ODE solvers
function orbital_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics_model)
end
```

## Low-Thrust Propulsion

Low-thrust models are frame-independent (RTN/VNB are orbital frames):

```julia
using AstroForceModels
using StaticArrays

# Constant tangential thrust in VNB frame (orbit raising)
lt_vnb = LowThrustAstroModel(;
    thrust_model = ConstantCartesianThrust(1e-7, 0.0, 0.0),  # [V, N, B]
    frame = VNBFrame(),
)

# Piecewise-constant thrust schedule in RTN
lt_rtn = LowThrustAstroModel(;
    thrust_model = PiecewiseConstantThrust(
        [0.0, 3600.0, 7200.0],
        [SVector{3}(0.0, 1e-7, 0.0),
         SVector{3}(0.0, 0.0, 0.0),
         SVector{3}(0.0, -1e-7, 0.0)],
    ),
    frame = RTNFrame(),
)

# 1 mN thrust on a 500 kg spacecraft
lt_simple = LowThrustAstroModel(;
    thrust_model = ConstantTangentialThrust(1e-3, 500.0),
)
```
