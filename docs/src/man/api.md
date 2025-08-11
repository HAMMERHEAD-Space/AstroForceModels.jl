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

See the main [AstroForceModels](@ref) module documentation for detailed descriptions of these abstract types.

## Gravity Models

### Types

```
GravityHarmonicsAstroModel
KeplerianGravityAstroModel
```

### Functions

```@docs  
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::GravityHarmonicsAstroModel)
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::KeplerianGravityAstroModel)
```

## Atmospheric Drag Models

### Types

```@docs
DragAstroModel
```

### Satellite Shape Models

```@docs
AbstractSatelliteDragModel
CannonballDragModel
```

### Functions

```@docs
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::DragAstroModel)
density_calculator
ballistic_coefficient
```

## Solar Radiation Pressure Models

### Types

```@docs
SRPAstroModel
```

### Satellite Shape Models

```@docs
AbstractSatelliteSRPModel
CannonballFixedSRP
StateSRPModel
```

### Shadow Models

```@docs
ShadowModelType
Conical
Cylindrical
```

### Functions

```@docs
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::SRPAstroModel)
shadow_function
radiation_pressure_coefficient
```

## Third Body Gravity Models

### Types

```@docs
ThirdBodyAstroModel
ThirdBodyModel
CelestialBody
```

### Functions

```@docs
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::ThirdBodyAstroModel)
third_body_position
```

## Relativistic Effects Models

### Types

```@docs
RelativisticAstroModel
```

### Functions

```@docs
acceleration(::AbstractArray, ::ComponentVector, ::Number, ::RelativisticAstroModel)
schwarzschild_acceleration  
lense_thirring_acceleration
de_sitter_acceleration
```

## Dynamics Models

The dynamics model system provides efficient ways to combine multiple force models:

### Types

```@docs
CentralBodyDynamicsModel
```

### Functions

```@docs
build_dynamics_model
```

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

## Type Definitions

### Component Arrays

The package uses ComponentArrays.jl for structured parameter handling:

```julia
using ComponentArrays

JD = 2.460310e6  # Julian Date for 2024-01-05 12:00:00
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
- [Force Model Documentation](../force_models/): Detailed descriptions of each force model
- Test suite: Working examples in the `/test` directory
- Benchmarks: Performance comparisons in the `/benchmark` directory

## Support and Development

### Contributing

Contributions are welcome! Please see the development guidelines:

1. **Fork the repository** and create a feature branch
2. **Write tests** for new functionality
3. **Follow coding style** conventions
4. **Document new features** with docstrings
5. **Submit a pull request** with a clear description

### Reporting Issues

Please report bugs and feature requests on the GitHub issue tracker:

- **Provide a minimal working example** that demonstrates the issue
- **Include version information** for Julia and all packages
- **Describe expected vs. actual behavior**
- **Include error messages** and stack traces if applicable

## Version History

See [CHANGELOG.md](../../../CHANGELOG.md) for detailed version history and breaking changes.