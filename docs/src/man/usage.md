# Usage

AstroForceModels.jl provides a comprehensive framework for modeling the dominant astrodynamics forces affecting satellite orbital motion. The package is designed to integrate seamlessly with the SatelliteToolbox.jl ecosystem and modern Julia differential equation solvers.

## Basic Concepts

### Force Models

Each force model in AstroForceModels implements the common `acceleration` interface:

```julia
acceleration(state, parameters, time, force_model)
```

Where:
- `state`: Current spacecraft state vector (position and velocity)
- `parameters`: Additional parameters (spacecraft properties, etc.)
- `time`: Current simulation time
- `force_model`: Specific force model instance

### Combining Force Models

Multiple force models can be combined using the `CentralBodyDynamicsModel` structure and `build_dynamics_model` function for efficient dynamics computation:

```julia
using AstroForceModels
using ComponentArrays
using DifferentialEquations

# Create individual force models
gravity_model = KeplerianGravityAstroModel()  # or GravityHarmonicsAstroModel
drag_model = DragAstroModel(...)
srp_model = SRPAstroModel(...)
third_body = ThirdBodyModel(; body=SunBody(), eop_data=fetch_iers_eop())

# Combine into a dynamics model
dynamics_model = CentralBodyDynamicsModel(
    gravity_model,
    (drag_model, srp_model, third_body)
)

# System dynamics function using the dynamics model
function orbital_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics_model)
end

# Alternative: Manual combination (less efficient)
function orbital_dynamics_manual!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = sum(acceleration(u, p, t, force) for force in [gravity_model, drag_model, srp_model, third_body_model])
end
```

## Complete Example: LEO Satellite

Here's a comprehensive example modeling a Low Earth Orbit satellite with multiple perturbations:

```julia
using AstroForceModels
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxGravityModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxBase
using ComponentArrays
using DifferentialEquations

# Initial conditions (ISS-like orbit)
r0 = [6378.137 + 408, 0.0, 0.0]  # Position [km]
v0 = [0.0, 7.660, 0.0]           # Velocity [km/s]
u0 = [r0; v0]

# Time span (24 hours)
tspan = (0.0, 86400.0)

# Spacecraft parameters
spacecraft_mass = 450000.0  # kg (ISS mass)
spacecraft_params = ComponentArray(
    mass = spacecraft_mass,
    Cd = 2.2,    # drag coefficient
    area = 1500.0, # cross-sectional area [m²]
    CR = 1.3     # radiation pressure coefficient
)

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(; JD=JD)

## 1. Gravity Model (with harmonics)
gravity_data = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
eop_data = fetch_iers_eop()

gravity_model = GravityHarmonicsAstroModel(
    gravity_model = gravity_data,
    eop_data = eop_data,
    order = 20,
    degree = 20
)

## 2. Atmospheric Drag Model
satellite_drag = CannonballFixedDrag(
    area = spacecraft_params.area,
    drag_coeff = spacecraft_params.Cd,
    mass = spacecraft_params.mass
)

drag_model = DragAstroModel(
    satellite_drag_model = satellite_drag,
    atmosphere_model = JR1971(),
    eop_data = eop_data
)

## 3. Solar Radiation Pressure Model
satellite_srp = CannonballFixedSRP(
    radius = sqrt(spacecraft_params.area / π),  # Convert area to radius
    mass = spacecraft_params.mass,
    reflectivity_coeff = spacecraft_params.CR
)

# Create Sun model
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

srp_model = SRPAstroModel(
    satellite_srp_model = satellite_srp,
    sun_data = sun_model,
    eop_data = eop_data,
    shadow_model = Conical()
)

## 4. Third Body Gravity (Sun and Moon)
# Sun model already created above for SRP
# Create Moon model  
moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

## 5. Relativistic Effects
relativity_model = RelativityModel(eop_data)

# Create a comprehensive dynamics model
dynamics_model = CentralBodyDynamicsModel(
    gravity_model,
    (drag_model, srp_model, sun_model, moon_model, relativity_model)
)

# System dynamics function for ODE solver
function satellite_dynamics!(du, u, p, t)
    du[1:3] = u[4:6]  # velocity
    du[4:6] = build_dynamics_model(u, p, t, dynamics_model)
end

# Solve the orbital dynamics
prob = ODEProblem(satellite_dynamics!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12)

# Alternative: Manual acceleration computation
total_accel = build_dynamics_model(u0, p, 0.0, dynamics_model)
```

## Troubleshooting

### Common Issues

1. **Coordinate System Consistency**: Ensure all models use the same reference frame
2. **Time System Handling**: Be consistent with time scales (UTC, TT, etc.)
3. **Unit Consistency**: Verify all quantities use consistent units (SI recommended)
4. **Numerical Stability**: Use appropriate integration tolerances
5. **Earth Orientation Parameters**: Keep EOP data current for high-precision work