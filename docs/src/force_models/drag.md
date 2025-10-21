# Drag

The drag force model accounts for the deceleration experienced by a spacecraft due to atmospheric resistance. As a spacecraft orbits Earth, it encounters atmospheric particles that exert a retarding force opposite to the spacecraft's velocity vector. This force is particularly significant for low Earth orbit (LEO) satellites, where atmospheric density is higher.

## Physical Description

Atmospheric drag arises from the collision of atmospheric molecules with the spacecraft's surface. The magnitude of the drag force depends on several factors:

- **Atmospheric density**: Varies with altitude, solar activity, and geomagnetic conditions
- **Spacecraft cross-sectional area**: The effective area presented to the velocity vector
- **Ballistic coefficient**: A measure combining spacecraft mass, drag coefficient, and cross-sectional area
- **Relative velocity**: The speed of the spacecraft relative to the atmosphere

The drag acceleration is given by:

```
a_drag = -0.5 * ρ * (v_rel · v_rel) * (A/m) * C_D * (v_rel / |v_rel|)
```

Where:
- `ρ` is atmospheric density
- `v_rel` is the spacecraft velocity relative to the atmosphere
- `A` is the cross-sectional area
- `m` is the spacecraft mass
- `C_D` is the drag coefficient

## Components

The drag force model in AstroForceModels consists of several key components:

### DragAstroModel

The main struct that encapsulates all drag-related parameters:

- **satellite_drag_model**: Defines the spacecraft's physical properties (area, drag coefficient, mass)
- **atmosphere_model**: Specifies which atmospheric model to use for density calculations
- **eop_data**: Earth Orientation Parameters for accurate atmospheric modeling

### Satellite Shape Models

Different representations of spacecraft geometry:

- **CannonballDragModel**: Simplified spherical model with constant drag coefficient

### Atmospheric Models

AstroForceModels supports five different atmospheric models for calculating atmospheric density:

#### JB2008 - Jacchia-Bowman 2008
- **Type**: Semi-empirical thermospheric model
- **Altitude Range**: Surface to 1,000 km
- **Features**: 
  - Includes solar activity effects (F10.7, Ap index)
  - Temperature and density variations
  - Based on satellite drag data
  - High accuracy for thermospheric modeling
- **Primary Use**: LEO satellites, precise orbit determination
- **Data Sources**: Satellite drag measurements, solar indices

#### JR1971 - Jacchia-Roberts 1971  
- **Type**: Static diffusion model with empirical temperature profiles
- **Altitude Range**: 90 to 2,500 km
- **Features**:
  - Includes latitudinal, seasonal, geomagnetic, and solar effects
  - Statistical accuracy of ~15%
  - Assumption of rigid body rotation with Earth
  - Temperature-dependent composition
- **Primary Use**: General orbital mechanics, extended altitude range
- **Data Sources**: Spacecraft drag data, atmospheric measurements

#### MSIS2000 - NRL MSIS 2000
- **Type**: Empirical global reference atmospheric model  
- **Altitude Range**: Surface to 1,000 km
- **Features**:
  - Mass spectrometer and incoherent scatter radar data
  - Comprehensive species composition (N₂, O₂, O, He, H, Ar, N)
  - Solar and geomagnetic activity dependence
  - International standard for space research
- **Primary Use**: High-precision applications, research
- **Data Sources**: Mass spectrometer, incoherent scatter radar, satellite data

#### ExpAtmo - Exponential Atmosphere
- **Type**: Simple exponential decay model
- **Altitude Range**: All altitudes (no upper limit)
- **Features**:
  - Assumes constant temperature with altitude
  - Exponential pressure/density decay: ρ = ρ₀ × e^(-h/H)
  - Scale height H ≈ 8.42 km
  - Very fast computation
- **Primary Use**: Quick analysis, initial estimates, educational purposes
- **Mathematical Form**: ρ(h) = ρ₀ × exp(-h/8420) where h is in meters

#### None - No Atmosphere
- **Type**: Vacuum model
- **Altitude Range**: All altitudes
- **Features**:
  - Returns zero atmospheric density
  - Eliminates drag effects entirely
  - Useful for comparative studies
- **Primary Use**: High-altitude missions, interplanetary trajectories, model comparison

## Model Selection Guidelines

### Choosing the Right Atmospheric Model

The choice of atmospheric model depends on your mission requirements:

**For High-Precision LEO Missions (< 600 km altitude):**
- **Recommended**: JB2008 or MSIS2000
- **Rationale**: Most accurate for thermospheric conditions
- **Applications**: ISS, precise orbit determination, collision avoidance

**For General LEO Applications (200-1000 km altitude):**
- **Recommended**: JB2008 or JR1971
- **Rationale**: Good balance of accuracy and computational efficiency
- **Applications**: CubeSats, standard orbital propagation

**For Extended Altitude Range (up to 2500 km):**
- **Recommended**: JR1971
- **Rationale**: Only model extending beyond 1000 km altitude
- **Applications**: HEO missions, Molniya orbits, disposal orbits

**For Rapid Analysis or Educational Use:**
- **Recommended**: ExpAtmo
- **Rationale**: Fastest computation, simple physics
- **Applications**: Mission planning, parametric studies, teaching

**For High-Altitude or Interplanetary Missions:**
- **Recommended**: None
- **Rationale**: Negligible atmospheric effects
- **Applications**: GEO satellites, lunar missions, deep space

### Computational Performance Comparison

| Model | Accuracy | Altitude Range | Memory Usage |
|-------|----------|----------------|--------------|
| JB2008 | Very High | 0-1000 km | Low |
| JR1971 | High | 90-2500 km | Low |
| MSIS2000 | Very High | 0-1000 km | Medium |
| ExpAtmo | Low | All altitudes | Minimal |
| None | Perfect* | All altitudes | Minimal |

*Perfect for vacuum conditions

## Usage Examples

### Basic Drag Model Setup

```julia
using AstroForceModels
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

# Define spacecraft properties
satellite_model = CannonballDragModel(
    area = 1.2,      # m²
    drag_coeff = 2.2, # dimensionless
    mass = 100.0     # kg
)

# Load Earth orientation parameters
eop_data = fetch_iers_eop()

# Create drag model with different atmospheric models
drag_model_jb2008 = DragAstroModel(
    satellite_drag_model = satellite_model,
    atmosphere_model = JB2008(),
    eop_data = eop_data
)

# Compute acceleration (typically called within integrator)
acceleration(state, parameters, time, drag_model_jb2008)
```

### Comparing Different Atmospheric Models

```julia
using AstroForceModels
using ComponentArrays

# Test state (ISS-like orbit)
state = [6378.137 + 408.0, 0.0, 0.0, 0.0, 7.6600, 0.0]  # [km, km/s]
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(; JD=JD)
time = 0.0

# Create models with different atmospheres
models = [
    DragAstroModel(satellite_model, JB2008(), eop_data),
    DragAstroModel(satellite_model, JR1971(), eop_data), 
    DragAstroModel(satellite_model, MSIS2000(), eop_data),
    DragAstroModel(satellite_model, ExpAtmo(), eop_data),
    DragAstroModel(satellite_model, None(), eop_data)
]

model_names = ["JB2008", "JR1971", "MSIS2000", "ExpAtmo", "None"]

# Compare accelerations
for (i, model) in enumerate(models)
    accel = acceleration(state, p, time, model)
    println("$(model_names[i]): $(norm(accel)) m/s²")
end
```

## Implementation Details

The drag acceleration computation involves:

1. **Density Calculation**: Using the specified atmospheric model and current spacecraft position
2. **Relative Velocity**: Computing spacecraft motion relative to the rotating atmosphere
3. **Cross-sectional Area**: Determining effective area based on velocity direction
4. **Force Calculation**: Computing the drag force and converting to acceleration

## References

[1] Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[2] Montenbruck, Oliver, and Eberhard Gill. "Satellite Orbits: Models, Methods and Applications." Springer, 2000.