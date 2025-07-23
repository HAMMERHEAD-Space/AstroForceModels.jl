# Solar Radiation Pressure

Solar Radiation Pressure (SRP) represents the momentum transfer from solar photons to a spacecraft's surface. When sunlight strikes a spacecraft, photons impart momentum that creates a small but persistent force. For spacecraft with large surface areas relative to their mass, or for missions requiring high precision, SRP can be a significant perturbing force.

## Physical Description

Solar radiation pressure arises when electromagnetic radiation from the Sun interacts with spacecraft surfaces. The interaction depends on the surface material properties:

- **Absorption**: Photons are absorbed, transferring all momentum to the surface
- **Specular reflection**: Photons reflect like a mirror, reversing momentum direction  
- **Diffuse reflection**: Photons scatter in all directions from the surface

The magnitude of SRP force depends on:

- **Solar flux**: Intensity of solar radiation at the spacecraft's location
- **Cross-sectional area**: Effective area exposed to solar radiation
- **Surface properties**: Reflectivity, absorptivity, and specularity coefficients
- **Sun-spacecraft distance**: Follows inverse square law with distance

The basic SRP acceleration is given by:

```
a_srp = -(Ψ/c) * (A/m) * C_R * (s/|s|)
```

Where:
- `Ψ` is solar flux at spacecraft distance
- `c` is speed of light  
- `A` is cross-sectional area perpendicular to Sun direction
- `m` is spacecraft mass
- `C_R` is radiation pressure coefficient (1 ≤ C_R ≤ 2)
- `s` is unit vector from spacecraft to Sun

## Shadow Effects

Spacecraft experience reduced or eliminated SRP when in shadow:

### Eclipse Types
- **Umbra**: Complete shadow (no direct sunlight)
- **Penumbra**: Partial shadow (partially blocked sunlight)
- **No shadow**: Full illumination

### Shadow Models

AstroForceModels implements three shadow models for eclipse calculations:

**Conical Shadow Model**: Sophisticated model accounting for Sun's finite size
- Models both umbra (complete shadow) and penumbra (partial shadow) effects
- Accounts for angular sizes of both Sun and Earth as seen from spacecraft
- Based on Montenbruck & Gill algorithm with geometric shadow calculations
- Most accurate shadow modeling available
- Computationally efficient despite sophistication

**Cylindrical Shadow Model**: Simplified model assuming parallel solar rays
- Sharp transition between sunlight (factor = 1.0) and complete shadow (factor = 0.0)
- Assumes Sun rays are parallel (infinite distance approximation)
- Faster computation than conical model
- Suitable for missions where high shadow accuracy isn't critical

**No_Shadow Model**: Disables all eclipse effects
- Always returns shadow factor = 1.0 (full sunlight)
- Eliminates computational overhead of shadow calculations
- Useful for high-altitude missions where eclipses are rare
- Essential for comparative studies and debugging

### Performance Comparison

| Shadow Model | Accuracy | Memory Usage | Best For |
|--------------|----------|--------------|----------|
| Conical | Highest | Low | LEO precision missions |
| Cylindrical | Medium | Minimal | Fast analysis, CubeSats |
| No_Shadow | N/A* | Minimal | High altitudes, debugging |

*Perfect accuracy for missions without eclipses

## Components

The SRP force model in AstroForceModels includes several components:

### SRPAstroModel

The main struct that encapsulates SRP parameters:

- **satellite_srp_model**: Spacecraft optical and geometric properties
- **sun_data**: Solar position calculation methods  
- **eop_data**: Earth orientation parameters
- **shadow_model**: Eclipse/shadow calculation method
- **R_Sun**: Solar radius (default: 695,700 km)
- **R_Occulting**: Occulting body radius (default: Earth radius)
- **Ψ**: Solar flux constant (default: 1361 W/m²)
- **AU**: Astronomical unit for distance scaling

### Satellite Shape Models

AstroForceModels provides two satellite shape models for SRP calculations:

**CannonballFixedSRP**: Spherical spacecraft model with fixed reflectivity coefficient
- Assumes spacecraft is a perfect sphere with radius `r`  
- Cross-sectional area: A = π × r²
- Fixed reflectivity coefficient throughout mission
- Reflectivity ballistic coefficient: RC = C_R × A / m
- Most common model for preliminary analysis and small satellites

### Surface Optical Properties

Key parameters defining surface interactions:

- **Absorptivity (α)**: Fraction of incident energy absorbed
- **Specularity (ρ_s)**: Fraction of reflected energy that reflects specularly  
- **Diffusivity (ρ_d)**: Fraction of reflected energy that reflects diffusely
- **Emissivity (ε)**: Ability to emit thermal radiation

Relationship: α + ρ_s + ρ_d = 1

## Usage Example

```julia
using AstroForceModels
using SatelliteToolboxCelestialBodies

# Define spacecraft optical properties using fixed reflectivity coefficient
satellite_model = CannonballFixedSRP(
    radius = 0.89,     # m (equivalent to 2.5 m² area)
    mass = 150.0,      # kg
    reflectivity_coeff = 1.3  # radiation pressure coefficient
)

# Load Earth orientation parameters  
eop_data = fetch_iers_eop()

# Create Sun position model
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

# Create SRP models with different shadow models

# High-precision model with penumbra effects
srp_model_conical = SRPAstroModel(
    satellite_srp_model = satellite_model,
    sun_data = sun_model,
    eop_data = eop_data,
    shadow_model = Conical(),    # Most accurate shadow model
)

# Compute acceleration (typically called within integrator)
acceleration(state, parameters, time, srp_model_conical)
```

## Implementation Details

The SRP acceleration computation involves:

1. **Solar Position**: Computing Sun's position relative to spacecraft
2. **Shadow Check**: Determining eclipse status using shadow model
3. **Solar Flux**: Calculating actual solar flux at spacecraft distance
4. **Surface Interactions**: Computing force based on optical properties
5. **Area Calculation**: Determining effective cross-sectional area
6. **Force Vector**: Computing acceleration vector in appropriate coordinate frame

## Radiation Pressure Coefficient

The radiation pressure coefficient C_R characterizes overall surface interaction:

- **C_R = 1**: Perfect absorption (black body)
- **1 < C_R < 2**: Mixed absorption and reflection (typical spacecraft)
- **C_R = 2**: Perfect specular reflection (ideal mirror)

Common values:
- Solar panels: C_R ≈ 1.3-1.5
- Spacecraft bodies: C_R ≈ 1.2-1.4  
- Thermal blankets: C_R ≈ 1.1-1.3

## Typical Magnitudes

SRP acceleration varies significantly with mission characteristics:

- **Large GEO satellites**: ~10⁻⁸ to 10⁻⁷ m/s²
- **Small satellites/CubeSats**: ~10⁻⁶ to 10⁻⁵ m/s²
- **Solar sails**: ~10⁻⁴ to 10⁻³ m/s²
- **GPS satellites**: ~10⁻⁸ m/s²

For comparison, other perturbations:
- Earth gravity (LEO): ~10⁻¹ m/s²
- Atmospheric drag (LEO): ~10⁻⁶ to 10⁻⁴ m/s²
- Third-body gravity: ~10⁻⁶ to 10⁻⁵ m/s²

## References

[1] Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[2] Montenbruck, Oliver, and Eberhard Gill. "Satellite Orbits: Models, Methods and Applications." Springer, 2000.
[3] AI Solutions. "Solar Radiation Pressure." FreeFlyer User Guide.
[4] Milani, A., et al. "Non-gravitational perturbations and satellite geodesy." 1987.
[5] Knocke, P. C., et al. "Earth radiation pressure effects on satellites." 1988. 