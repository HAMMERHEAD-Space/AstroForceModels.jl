# API Reference

This page provides a comprehensive reference for the AstroForceModels.jl API. All force models implement the common `acceleration` interface and can be combined to create comprehensive orbital dynamics simulations.

## Core Interface

### acceleration Function

The primary interface for all force models:

```@docs
acceleration
```

## Abstract Types

AstroForceModels defines a hierarchy of abstract types to organize different force models:

- `AbstractAstroForceModel`: Base type for all force models
- `AbstractNonPotentialBasedForce`: Non-conservative forces (drag, SRP, etc.)  
- `AbstractPotentialBasedForce`: Conservative forces derivable from potential (gravity)
- `AbstractDynamicsModel`: Base type for dynamics model combinations

## Dynamics Models

The dynamics model system provides efficient ways to combine multiple force models:

```@docs
build_dynamics_model
```

The `CentralBodyDynamicsModel` type is documented in the [Library](@ref) section.

### Usage Example

```julia
# Create individual force models
gravity_model = GravityHarmonicsAstroModel(...)
drag_model = DragAstroModel(...)
srp_model = SRPAstroModel(...)

# Combine into dynamics model
dynamics_model = CentralBodyDynamicsModel(
    gravity_model,
    (drag_model, srp_model)
)

# Use in ODE system
function dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics_model)
end
```

## Force Models

AstroForceModels provides several categories of force models:

### Gravity Models

- `GravityHarmonicsAstroModel`: Spherical harmonics gravity model with EGM96/EGM2008 coefficients
- `KeplerianGravityAstroModel`: Simple point-mass gravity model

### Atmospheric Drag Models

- `DragAstroModel`: Atmospheric drag force model with various atmospheric density models
- `CannonballFixedDrag`: Fixed ballistic coefficient satellite shape model

### Solar Radiation Pressure Models

- `SRPAstroModel`: Solar radiation pressure force model with shadow modeling
- `CannonballFixedSRP`: Fixed reflectivity coefficient satellite shape model

### Albedo Radiation Pressure Models

- `AlbedoAstroModel`: Earth albedo radiation pressure force model with Lebedev quadrature integration
- `UniformAlbedoModel`: Uniform albedo model with constant visible albedo and infrared emissivity

### Third Body Models

- `ThirdBodyModel`: Third body ephemeris provider (Sun, Moon, planets)

### Relativistic Effects

- `RelativityModel`: General relativistic effects (Schwarzschild, Lense-Thirring, de Sitter)

### Solid Body Tides

- `SolidBodyTidesModel`: Tidal deformation perturbation from tide-raising bodies (IERS 2010, Step 1)

### Low-Thrust Propulsion

- `LowThrustAstroModel`: Low-thrust propulsion force model with frame support
- `ConstantCartesianThrust`: Fixed thrust acceleration vector (frame-dependent)
- `ConstantTangentialThrust`: Constant-magnitude velocity-aligned thrust
- `StateThrustModel`: User-defined state-dependent thrust profile
- `PiecewiseConstantThrust`: Piecewise-constant arc-based thrust schedule

### Reference Frames

- `InertialFrame`: ECI/J2000 inertial frame (default)
- `RTNFrame`: Radial–Transverse–Normal orbit-fixed frame
- `VNBFrame`: Velocity–Normal–Binormal velocity-fixed frame

All force models implement the common `acceleration(state, params, time, model)` interface.

## Type Definitions

### Component Arrays

The package uses ComponentArrays.jl for structured parameter handling:

```julia
using ComponentArrays
using SatelliteToolboxBase

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(; JD=JD)
```

### State Vectors

State vectors are typically 6-element arrays representing position and velocity:

```julia
# State vector format: [rx, ry, rz, vx, vy, vz]
state = [
    6.378137e3 + 400,    # x position [km]
    0.0,                 # y position [km]  
    0.0,                 # z position [km]
    0.0,                 # x velocity [km/s]
    7.660,               # y velocity [km/s]
    0.0                  # z velocity [km/s]
]
```

## Integration with SatelliteToolbox.jl

AstroForceModels is designed to be built off packages in the SatelliteToolbox.jl ecosystem:

### Required Dependencies

```julia
using SatelliteToolboxBase           # Base types and constants
using SatelliteToolboxGravityModels  # Gravity field models
using SatelliteToolboxAtmosphericModels  # Atmospheric density models
using SatelliteToolboxCelestialBodies    # Celestial body ephemeris
using SatelliteToolboxTransformations    # Coordinate transformations
```

### Earth Orientation Parameters

```julia
# Fetch latest EOP data
eop_data = fetch_iers_eop()

# Use with force models
gravity_model = GravityHarmonicsAstroModel(
    gravity_model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96)),
    eop_data = eop_data,
    order = 20,
    degree = 20
)
```

## Examples and Tutorials

For practical examples and step-by-step tutorials, see:

- [Usage Guide](usage.md): Comprehensive usage examples
- Test suite: Working examples in the `/test` directory
- Benchmarks: Performance comparisons in the `/benchmark` directory