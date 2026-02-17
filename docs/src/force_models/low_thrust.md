# Low-Thrust Propulsion

Low-thrust propulsion models the continuous acceleration produced by electric propulsion systems (ion engines, Hall thrusters, etc.), solar sails, or any other source of small, sustained thrust. Unlike impulsive maneuvers that instantaneously change spacecraft velocity, low-thrust systems apply small forces over extended periods to achieve large cumulative velocity changes with high propellant efficiency.

## Physical Description

Low-thrust propulsion is characterized by:

- **Small acceleration magnitudes**: Typically 10⁻⁷ to 10⁻⁴ km/s² (milli-Newton to Newton level thrust)
- **Continuous operation**: Thrust applied over hours, days, or months
- **High specific impulse**: Electric propulsion systems achieve 1,000–10,000 s Isp vs. 300–450 s for chemical
- **Variable thrust direction**: Can be steered to optimize trajectory

The thrust acceleration acting on a spacecraft is given by:

```
a_thrust = (T / m) × d̂
```

Where:
- `T` is the thrust magnitude [N]
- `m` is the spacecraft mass [kg]
- `d̂` is the thrust direction unit vector

Note that the acceleration is divided by 1000 to convert from m/s² to km/s² for consistency with the library's unit convention.

## Reference Frames

Thrust can be specified in three reference frames via the `frame` field on `LowThrustAstroModel`:

### InertialFrame (default)

Components are expressed in the ECI/J2000 inertial frame. No rotation is applied.

### RTNFrame (Radial–Transverse–Normal)

Also known as RIC (Radial–In-track–Cross-track). Components `[a_R, a_T, a_N]` are defined by:

- **R̂** (Radial): `r / |r|` — away from the central body
- **T̂** (Transverse): `N̂ × R̂` — roughly along-track for near-circular orbits
- **N̂** (Normal): `(r × v) / |r × v|` — along the orbital angular momentum

This frame is intuitive for orbit maneuvers: transverse thrust raises/lowers the orbit, radial thrust changes eccentricity, and normal thrust changes inclination.

### VNBFrame (Velocity–Normal–Binormal)

Components `[a_V, a_N, a_B]` are defined by:

- **V̂** (Velocity): `v / |v|` — along the velocity vector
- **N̂** (Normal): `B̂ × V̂` — in the orbital plane, perpendicular to velocity
- **B̂** (Binormal): `(r × v) / |r × v|` — along the orbital angular momentum

This frame is natural for velocity-aligned thrust (orbit raising/lowering) and out-of-plane maneuvers.

!!! note
    For circular orbits, VNB and RTN coincide (V̂ ≈ T̂). The distinction matters for eccentric orbits where the velocity direction deviates from the transverse direction.

## Components

### LowThrustAstroModel

The main struct for low-thrust computations, wrapping an [`AbstractThrustModel`] and an [`AbstractThrustFrame`]:

- **thrust_model**: An instance of `AbstractThrustModel` specifying how the thrust acceleration is computed
- **frame**: An instance of `AbstractThrustFrame` specifying the reference frame (default: `InertialFrame()`)

### Thrust Models

Four concrete thrust model types are provided:

#### ConstantCartesianThrust

A fixed 3-component acceleration vector, interpreted in the frame specified by `LowThrustAstroModel`. Suitable for:
- Constant station-keeping thrust
- Fixed pointing in any frame (inertial, RTN, or VNB)
- Simple analysis and testing

**Constructors:**
- `ConstantCartesianThrust(ax, ay, az)`: Direct acceleration components [km/s²]
- `ConstantCartesianThrust(direction, thrust, mass)`: From direction vector, thrust [N], and mass [kg]

#### ConstantTangentialThrust

Constant-magnitude thrust directed along the velocity vector. Always returns the acceleration in the inertial frame (should be used with `InertialFrame()`). For velocity-aligned thrust in other frames, use `ConstantCartesianThrust(magnitude, 0, 0)` with `VNBFrame()`.

**Constructors:**
- `ConstantTangentialThrust(magnitude)`: Direct acceleration magnitude [km/s²]
- `ConstantTangentialThrust(thrust, mass)`: From thrust [N] and mass [kg]

#### StateThrustModel

A fully user-defined thrust model using a callable `f(u, p, t) → SVector{3}`. The output is interpreted in the model's frame. Supports:
- Feedback control laws (e.g., Q-Law, Lyapunov controllers)
- Scheduled thrust arcs (on/off profiles)
- Optimization-based thrust histories

#### PiecewiseConstantThrust

A piecewise-constant thrust schedule defined by a sequence of time-tagged arcs. Each arc
specifies a constant 3-component acceleration vector that is active from its start time
until the next arc begins.

- `times`: Strictly increasing vector of arc start times [s]
- `accelerations`: Corresponding acceleration vectors [km/s²] for each arc

This is the natural representation for Sims-Flanagan trajectory segments, finite-burn
maneuver schedules, and any thrust profile that changes discretely between arcs.

## Usage Examples

### Tangential Thrust (Orbit Raising)

```julia
using AstroForceModels
using ComponentArrays

thrust = ConstantTangentialThrust(1e-3, 500.0)  # 1 mN, 500 kg
lt_model = LowThrustAstroModel(; thrust_model=thrust)
```

### Transverse Thrust in RTN Frame

```julia
# Constant transverse thrust for orbit raising
lt_rtn = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0),  # [R, T, N]
    frame=RTNFrame(),
)
```

### Velocity-Aligned Thrust in VNB Frame

```julia
# Equivalent to ConstantTangentialThrust but explicitly in VNB
lt_vnb = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0),  # [V, N, B]
    frame=VNBFrame(),
)
```

### Piecewise-Constant Thrust Schedule (Sims-Flanagan Style)

```julia
using StaticArraysCore

# Three-arc Sims-Flanagan segment in RTN frame
lt_piecewise = LowThrustAstroModel(;
    thrust_model=PiecewiseConstantThrust(
        [0.0, 3600.0, 7200.0],                       # arc start times [s]
        [
            SVector{3}(0.0, 1e-7, 0.0),              # arc 1: transverse thrust
            SVector{3}(0.0, 0.0, 0.0),               # arc 2: coast
            SVector{3}(0.0, -1e-7, 0.0),             # arc 3: reverse transverse
        ],
    ),
    frame=RTNFrame(),
)
```

### Thrust-Coast-Thrust Maneuver

```julia
# Thrust for 1 orbit, coast, then thrust again
lt_tct = LowThrustAstroModel(;
    thrust_model=PiecewiseConstantThrust(
        [0.0, 5400.0, 10800.0],
        [
            SVector{3}(0.0, 1e-7, 0.0),              # burn 1
            SVector{3}(0.0, 0.0, 0.0),               # coast
            SVector{3}(0.0, 1e-7, 0.0),              # burn 2
        ],
    ),
    frame=RTNFrame(),
)
```

### Combining with Other Force Models

```julia
using AstroForceModels
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations

eop_data = fetch_iers_eop()
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
gravity = GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, eop_data=eop_data, order=20, degree=20
)

lt_model = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0),
    frame=RTNFrame(),
)

dynamics = CentralBodyDynamicsModel(gravity, (lt_model,))

function low_thrust_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]
    du[4:6] = build_dynamics_model(u, p, t, dynamics)
end
```

## Typical Thrust Magnitudes

| Propulsion System | Thrust Range | Specific Impulse | Acceleration (500 kg s/c) |
|---|---|---|---|
| Hall thruster | 40–600 mN | 1,000–3,000 s | 8×10⁻⁸ to 1.2×10⁻⁶ km/s² |
| Ion engine (NSTAR) | 20–90 mN | 1,000–3,100 s | 4×10⁻⁸ to 1.8×10⁻⁷ km/s² |
| Electrospray | 0.01–1 mN | 500–5,000 s | 2×10⁻¹¹ to 2×10⁻⁹ km/s² |
| Solar sail (1 AU) | ~9 μN/m² | ∞ | Varies with area |
| Pulsed plasma | 0.1–10 mN | 500–2,000 s | 2×10⁻¹⁰ to 2×10⁻⁸ km/s² |

## Implementation Details

The low-thrust acceleration computation follows two steps:

1. **Thrust model evaluation**: The `AbstractThrustModel` produces a 3-component acceleration vector in the model frame
2. **Frame transformation**: `transform_to_inertial` rotates the vector from the specified frame (Inertial, RTN, or VNB) to the ECI frame using basis vectors constructed from the spacecraft state

The frame transformation constructs orthonormal basis vectors directly from `r` and `v`, avoiding intermediate matrix construction for efficiency. All operations use `SVector{3}` for type stability and are marked `@inline` for performance. The `InertialFrame` transformation is a compile-time no-op.

## References

[1] Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications* (4th ed.). Microcosm Press.
[2] Conway, B. A. (Ed.) (2010). *Spacecraft Trajectory Optimization*. Cambridge University Press.
[3] Betts, J. T. (2010). *Practical Methods for Optimal Control and Estimation Using Nonlinear Programming* (2nd ed.). SIAM.
[4] Goebel, D. M. and Katz, I. (2008). *Fundamentals of Electric Propulsion: Ion and Hall Thrusters*. Wiley.
[5] Schaub, H. and Junkins, J. L. (2018). *Analytical Mechanics of Space Systems* (4th ed.). AIAA Education Series.
[6] Sims, J. A. and Flanagan, S. N. (1999). "Preliminary Design of Low-Thrust Interplanetary Missions." AAS/AIAA Astrodynamics Specialist Conference, AAS 99-338.
