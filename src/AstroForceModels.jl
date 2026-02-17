"""
    AstroForceModels

A comprehensive Julia package for modeling astrodynamics forces affecting satellite orbital motion.

This package provides a unified framework for computing accelerations from various perturbative 
forces including atmospheric drag, solar radiation pressure, third-body gravity, relativistic 
effects, and gravitational harmonics. All force models implement a common `acceleration` interface
and can be efficiently combined using the `CentralBodyDynamicsModel` system.

# Key Features

- **Unified Interface**: All force models implement `acceleration(state, params, time, model)`
- **Efficient Combination**: Use `CentralBodyDynamicsModel` and `build_dynamics_model` for optimal performance
- **SatelliteToolbox Integration**: Built on the SatelliteToolbox.jl ecosystem
- **Automatic Differentiation**: All models support ForwardDiff.jl for gradient-based optimization
- **High Performance**: Type-stable implementations with compile-time optimizations

# Force Models Available

- **Gravity Models**: Keplerian and spherical harmonics (via SatelliteToolboxGravityModels)
- **Atmospheric Drag**: Multiple atmospheric models (JR1971, JB2008, NRLMSISE00, etc.)
- **Solar Radiation Pressure**: With shadow modeling (conical, cylindrical)
- **Third-Body Gravity**: Sun, Moon, and planetary perturbations
- **Relativistic Effects**: Schwarzschild, Lense-Thirring, and de Sitter effects
- **Low-Thrust Propulsion**: Constant, tangential, and user-defined thrust profiles

# Example Usage

```julia
using AstroForceModels

# Create individual force models
gravity = GravityHarmonicsAstroModel(...)
drag = DragAstroModel(...)
srp = SRPAstroModel(...)

# Combine efficiently
dynamics = CentralBodyDynamicsModel(gravity, (drag, srp))

# Use in ODE integration
function satellite_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics)
end
```

# See Also

- [Usage Guide](https://astroforcemodels.jl.org/dev/man/usage/): Comprehensive examples
- [API Reference](https://astroforcemodels.jl.org/dev/man/api/): Complete function documentation
- [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl): Ecosystem foundation
"""
module AstroForceModels

using ComponentArrays, StaticArraysCore
using LinearAlgebra
using SatelliteToolboxBase
using SatelliteToolboxCelestialBodies
using SatelliteToolboxGravityModels
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxTransformations
using SpaceIndices

"""
    AbstractAstroForceModel

Abstract base type for all astrodynamics force models.

All force models must implement the `acceleration(state, params, time, model)` interface
to compute the 3-dimensional acceleration vector acting on a spacecraft.

# Subtypes
- [`AbstractPotentialBasedForce`](@ref): Forces derivable from a potential (gravity)
- [`AbstractNonPotentialBasedForce`](@ref): Non-conservative forces (drag, SRP, etc.)
"""
abstract type AbstractAstroForceModel end

"""
    AbstractNonPotentialBasedForce <: AbstractAstroForceModel

Abstract type for non-conservative force models that cannot be derived from a potential.

Examples include atmospheric drag, solar radiation pressure, and relativistic effects.
These forces typically depend on velocity and other non-conservative factors.
"""
abstract type AbstractNonPotentialBasedForce <: AbstractAstroForceModel end

"""
    AbstractPotentialBasedForce <: AbstractAstroForceModel

Abstract type for conservative force models derivable from a gravitational potential.

Examples include point-mass gravity, gravitational harmonics, and third-body gravity.
These forces depend only on position and can be derived from a potential function.
"""
abstract type AbstractPotentialBasedForce <: AbstractAstroForceModel end

"""
    AbstractDynamicsModel

Abstract base type for dynamics models that combine multiple force models.

Dynamics models provide efficient ways to compute total acceleration from multiple
force sources. The primary implementation is [`CentralBodyDynamicsModel`](@ref).
"""
abstract type AbstractDynamicsModel end

include("constants.jl")
include("utils.jl")

include("force_models/drag/satellite_shape_model.jl")
include("force_models/drag/density_calculator.jl")
include("force_models/drag/drag_accel.jl")

include("force_models/third_body/celestial_body.jl")
include("force_models/third_body/third_body_model.jl")
include("force_models/third_body/third_body_accel.jl")

include("force_models/relativity/relativity_accel.jl")

include("force_models/solar_radiation_pressure/satellite_shape_model.jl")
include("force_models/solar_radiation_pressure/shadow_models.jl")
include("force_models/solar_radiation_pressure/srp_accel.jl")

include("force_models/gravity/utils.jl")
include("force_models/gravity/gravity_accel.jl")

include("force_models/low_thrust/frames.jl")
include("force_models/low_thrust/thrust_model.jl")
include("force_models/low_thrust/low_thrust_accel.jl")

include("dynamics_builder.jl")

export acceleration,
    potential,
    potential_time_derivative,
    AbstractAstroForceModel,
    AbstractDynamicsModel,
    AbstractNonPotentialBasedForce,
    AbstractPotentialBasedForce

end
