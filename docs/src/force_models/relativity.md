# Relativistic Effects

Relativistic effects account for corrections to Newtonian mechanics required by Einstein's theories of special and general relativity. While these effects are typically small for most spacecraft applications, they become increasingly important for high-precision missions, navigation satellites, and scientific applications requiring extreme accuracy.

## Physical Description

Relativistic effects arise from two main sources:

### Special Relativity
- **Time dilation**: Moving clocks run slower
- **Length contraction**: Objects contract along direction of motion  
- **Mass-energy equivalence**: E = mc²
- **Velocity dependence**: Effects scale with v²/c²

### General Relativity  
- **Gravitational time dilation**: Clocks run slower in stronger gravitational fields
- **Curved spacetime**: Gravity as geometry rather than force
- **Frame dragging**: Rotation of massive bodies affects nearby spacetime
- **Geodetic precession**: Precession of gyroscopes in curved spacetime

## Relativistic Acceleration Terms

The post-Newtonian relativistic acceleration includes several terms:

```
a_rel = a_schwarzschild + a_lense_thirring + a_de_sitter + a_geodetic
```

### Schwarzschild Term (Static Gravitational Field)
The primary relativistic correction for a static, spherically symmetric gravitational field:

```
a_schw = (μ/c²r³) * [4μ/r - v² + 4(r·v)²/r²] * r + 4μ(r·v)/(c²r³) * v
```

Where:
- `μ` is the gravitational parameter
- `c` is the speed of light
- `r` is position vector  
- `v` is velocity vector
- `·` denotes dot product

### Lense-Thirring Effect (Frame Dragging)
Caused by the rotation of the central body (Earth):

```
a_LT = (2μ/c²r³) * (3(J×r)/r² × v + (J×v)/r²)
```

Where:
- `J` is the angular momentum vector of the central body
- `×` denotes cross product

### De Sitter Precession
Due to the curvature of spacetime caused by external masses (primarily the Sun):

```
a_dS = (μ_sun/c²r_sun³) * [3(r_sun·r_sat)/r_sun² * r_sun - r_sat]
```

## Components

The relativistic force model in AstroForceModels provides:

### RelativityModel

The main struct for relativistic computations:

- **schwarzschild_effect**: Include primary relativistic correction (default: true)
- **lense_thirring_effect**: Include frame dragging effects (default: true)
- **de_Sitter_effect**: Include external body effects (default: true)
- **central_body**: ThirdBodyModel for the central body (default: Earth)
- **sun_body**: ThirdBodyModel for the Sun (for de Sitter effects)
- **c**: Speed of light [km/s]
- **γ**, **β**: PPN parameters (default: 1.0)

## Usage Example

```julia
using AstroForceModels
using SatelliteToolboxTransformations

# Load EOP data once and share across sub-models
eop_data = fetch_iers_eop()

# Create relativistic model with all effects using convenience constructor
rel_model = RelativityModel(eop_data)

# Or selectively enable effects
rel_schwarzschild_only = RelativityModel(eop_data;
    schwarzschild_effect = true,
    lense_thirring_effect = false,
    de_Sitter_effect = false
)

# Compute acceleration (typically called within integrator)
acceleration(state, parameters, time, rel_model)
```

## Magnitude of Effects

Relativistic corrections vary significantly with orbital characteristics:

### Low Earth Orbit (LEO, ~400 km)
- **Schwarzschild**: ~10⁻⁹ m/s² (small but measurable)
- **Lense-Thirring**: ~10⁻¹¹ m/s² (very small)
- **Applications**: Gravity Probe B, precise orbit determination

### Medium Earth Orbit (MEO, ~20,000 km) 
- **Schwarzschild**: ~10⁻¹⁰ m/s² (significant for precise timing)
- **Applications**: GPS satellites (μs timing accuracy required)

### Geostationary Orbit (GEO, ~35,786 km)
- **Schwarzschild**: ~10⁻¹⁰ m/s² (long-term orbital evolution)
- **Applications**: Communication satellites, long-term tracking

## Implementation Details

Relativistic acceleration computation involves:

1. **Reference Frame**: Ensure consistent coordinate system
2. **Time Systems**: Proper handling of different time scales
3. **Velocity Calculation**: Accurate spacecraft velocity determination
4. **Central Body Properties**: Mass, angular momentum, multipole moments
5. **External Body Positions**: For third-body relativistic effects
6. **Numerical Precision**: High-precision arithmetic may be required

## References

[1] Will, Clifford M. "Theory and Experiment in Gravitational Physics." Cambridge University Press, 2018.
[2] Soffel, Michael H., et al. "The IAU 2000 Resolutions for Astrometry, Celestial Mechanics, and Metrology in the Relativistic Framework." Astronomical Journal, 2003.
[3] Petit, Gérard, and Brian Luzum, eds. "IERS Conventions (2010)." IERS Technical Note 36, 2010.
[4] Ashby, Neil. "Relativity in the Global Positioning System." Living Reviews in Relativity, 2003.
[5] Kopeikin, Sergei, et al. "Post-Newtonian approximations for the propagation of light in the field of moving bodies." Physical Review D, 2002. 