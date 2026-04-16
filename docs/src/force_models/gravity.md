# Gravity / Spherical Harmonics

The gravity force model accounts for the gravitational acceleration experienced by a spacecraft due to the non-uniform mass distribution of the central body (typically Earth). While a point-mass gravity model assumes the central body is perfectly spherical and uniform, real celestial bodies have irregular mass distributions that create additional gravitational perturbations.

## Physical Description

The Earth's gravitational field deviates from that of a perfect sphere due to:

- **Oblateness**: The Earth is flattened at the poles due to rotation
- **Mass irregularities**: Variations in density and topography
- **Tidal deformation**: Deformation due to external gravitational forces

These irregularities are mathematically represented using spherical harmonics, which decompose the gravitational potential into a series expansion. The most significant terms are the zonal harmonics (particularly J₂), which account for the Earth's oblateness.

The gravitational acceleration including harmonics is:

```
a_grav = -∇U = -∇(μ/r + ΔU)
```

Where:
- `μ` is the gravitational parameter
- `r` is the distance from the center of mass
- `ΔU` represents the harmonic corrections to the potential

## Harmonic Coefficients

The spherical harmonic expansion uses coefficients Cₙₘ and Sₙₘ where:

- **n** is the degree (order of the harmonic)
- **m** is the order (number of nodal lines)
- **Zonal harmonics** (m=0): Symmetric about the rotation axis
- **Tesseral harmonics** (m≠0): Asymmetric terms

## Components

The gravity force model in AstroForceModels provides two main implementations:

### GravityHarmonicsAstroModel

A comprehensive model that includes spherical harmonic terms:

- **gravity_model**: Contains the harmonic coefficients and reference data
- **body_fixed_frame**: The body-fixed frame symbol (e.g., `:ITRF` for Earth, `:ErosBodyFixed` for a small body)
- **propagation_frame**: The propagation/inertial frame symbol (default: `:ICRF`)
- **order**: Maximum order of harmonics to compute (-1 for maximum available)
- **degree**: Maximum degree of harmonics to compute (-1 for maximum available)
- **P**, **dP**: Optional pre-allocated Legendre polynomial buffers

### KeplerianGravityAstroModel

A simplified point-mass gravity model for comparison or computational efficiency:

- **μ**: Gravitational parameter of the central body [km³/s²] (required, no default)
- Assumes perfectly spherical, uniform central body
- Only includes the μ/r² term
- Useful for initial orbit determination or when high precision isn't required

## Usage Example

```julia
using AstroForceModels
using SatelliteToolboxGravityModels

# Load a gravity model (e.g., EGM96)
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

# Create high-fidelity gravity model (Earth)
gravity_model = GravityHarmonicsAstroModel(
    gravity_model = grav_coeffs,
    body_fixed_frame = :ITRF,
    propagation_frame = :ICRF,
    order = 20,    # Use up to degree/order 20
    degree = 20
)

# Create Keplerian (point-mass) gravity model
kepler_gravity = KeplerianGravityAstroModel(μ = 398600.4415)

# Set up frame-aware parameters (see examples for full frame setup)
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = FrameAwareParams(frames, epoch, :ICRF)

# Compute acceleration (typically called within integrator)
acceleration(state, p, time, gravity_model)
```

## Gravity Models Available

Common gravity models supported:

- **EGM96**: Earth Gravitational Model 1996 (360×360)
- **EGM2008**: High-resolution model (2190×2190)
- **GGM03C**: GRACE-based model
- **JGM-3**: Joint Gravity Model 3
- **WGS84**: World Geodetic System 1984 reference

Any ICGEM format is accepted, see SatelliteToolboxGravityModels.jl for details.

## Implementation Details

The harmonic gravity computation involves:

1. **Coordinate Transformation**: Converting position to body-fixed coordinates
2. **Legendre Polynomials**: Computing associated Legendre functions
3. **Harmonic Summation**: Evaluating the spherical harmonic series
4. **Gradient Calculation**: Computing spatial derivatives for acceleration
5. **Coordinate Transformation**: Converting back to inertial frame

## References

[1] Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[2] Montenbruck, Oliver, and Eberhard Gill. "Satellite Orbits: Models, Methods and Applications." Springer, 2000.
[3] Hofmann-Wellenhof, B., and H. Moritz. "Physical Geodesy." 2nd ed., Springer, 2006.
[4] Lemoine, F. G., et al. "The Development of the Joint NASA GSFC and the National Imagery and Mapping Agency (NIMA) Geopotential Model EGM96." NASA Technical Publication, 1998. 