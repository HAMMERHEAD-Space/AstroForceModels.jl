# Earth Albedo Radiation Pressure

Earth Albedo Radiation Pressure (ERP) represents the momentum transfer from Earth-reflected and Earth-emitted radiation to a spacecraft. Unlike direct solar radiation pressure, albedo forces originate from the Earth's surface and are generally smaller in magnitude but can be significant for precision orbit determination, especially for low-Earth orbit (LEO) missions such as GRACE.

## Physical Description

Earth albedo radiation pressure arises from two distinct sources:

- **Shortwave (reflected) radiation**: Solar radiation reflected by the Earth's surface and atmosphere. The reflected fraction depends on the local albedo coefficient, which varies with surface type (ocean, desert, ice, clouds).
- **Longwave (thermal) radiation**: Infrared radiation emitted by the Earth due to its thermal equilibrium. This emission is approximately isotropic and present on both the sunlit and dark sides.

The total acceleration on the spacecraft is computed by integrating the radiation contributions over the visible portion of the Earth's surface (field of view), using Lebedev quadrature for efficient spherical integration.

The mathematical formulation follows Knocke et al. (1988):

```
a_albedo = RC * integral over visible surface of [F_total * cos(alpha) / (pi * c * r^2)] dOmega
```

Where:
- `RC` is the reflectivity ballistic coefficient (area-to-mass ratio times reflectivity)
- `F_total` is the sum of shortwave and longwave fluxes from each surface element
- `alpha` is the angle between the surface element normal and the satellite direction
- `c` is the speed of light
- `r` is the distance from the surface element to the satellite
- `dOmega` is the surface element area

## Albedo Models

### UniformAlbedoModel

The simplest model, assuming constant albedo and emissivity across the entire Earth's surface.

**Parameters:**
- `visible_albedo`: Fraction of solar radiation reflected (default: 0.3)
- `infrared_emissivity`: Fraction of thermal radiation emitted (default: 0.7)

The default values represent commonly accepted Earth-averaged values from the radiation budget literature.

## Components

### AlbedoAstroModel

The main struct that encapsulates albedo radiation pressure parameters:

- **satellite_shape_model**: Spacecraft optical and geometric properties (uses same `AbstractSatelliteSRPModel` as SRP)
- **sun_data**: Solar position calculation via `ThirdBodyModel`
- **body_albedo_model**: Earth albedo model (e.g., `UniformAlbedoModel`)
- **eop_data**: Earth orientation parameters for ECEF/ECI transformations
- **solar_flux**: Solar flux at 1 AU [W/m^2]
- **speed_of_light**: Speed of light [km/s]
- **lebedev_order**: Order of Lebedev quadrature (default: 125). Higher orders give more integration points and better accuracy at the cost of computation time.

### Lebedev Quadrature Orders

The integration accuracy and computational cost scale with the Lebedev quadrature order:

| Order | Points | Accuracy | Best For |
|-------|--------|----------|----------|
| 21    | 170    | Low      | Fast preliminary analysis |
| 59    | 590    | Medium   | General-purpose analysis |
| 125   | 1202   | High     | Precision orbit determination |

## Usage Example

```julia
using AstroForceModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxTransformations

# Load Earth orientation parameters
eop_data = fetch_iers_eop()

# Create Sun position model
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

# Define spacecraft shape model (same as SRP)
satellite_model = CannonballFixedSRP(0.030)  # m^2/kg

# Create uniform albedo model
uniform_albedo = UniformAlbedoModel(; visible_albedo=0.3, infrared_emissivity=0.7)

# Create albedo force model
albedo_model = AlbedoAstroModel(;
    satellite_shape_model=satellite_model,
    sun_data=sun_model,
    body_albedo_model=uniform_albedo,
    eop_data=eop_data,
    lebedev_order=125,  # High accuracy
)

# Compute acceleration (typically called within an integrator)
acceleration(state, parameters, time, albedo_model)
```

## Combining with Other Forces

The albedo model integrates naturally with the `CentralBodyDynamicsModel`:

```julia
gravity = GravityHarmonicsAstroModel(...)
srp = SRPAstroModel(...)
albedo = AlbedoAstroModel(...)

dynamics = CentralBodyDynamicsModel(gravity, (srp, albedo))
accel = build_dynamics_model(u, p, t, dynamics)
```

## Typical Magnitudes

Albedo acceleration is typically 1-2 orders of magnitude smaller than direct SRP:

- **LEO satellites (GRACE-like)**: ~10^-11 km/s^2 (~10^-8 m/s^2)
- **Higher altitudes**: Decreases approximately with inverse square of altitude

For comparison:
- Direct SRP: ~10^-9 to 10^-7 km/s^2
- Atmospheric drag (LEO): ~10^-9 to 10^-7 km/s^2
- Third-body gravity: ~10^-9 to 10^-8 km/s^2

## References

[1] Knocke, P. C., Ries, J. C., and Tapley, B. D. (1988). "Earth radiation pressure effects on satellites." Proceedings of AIAA/AAS Astrodynamics Conference, pp. 577-587.

[2] Borderies, N., & Longaretti, P. Y. (1990). "A new treatment of the albedo radiation pressure in the case of a uniform albedo and of a spherical satellite." Celestial Mechanics and Dynamical Astronomy, 49(1), 69-98.

[3] Rubincam, D. P., & Weiss, N. R. (1986). "Earth albedo and the orbit of Lageos." Celestial Mechanics, 38(3), 233-296.

[4] Vielberg, K., & Kusche, J. (2020). "Extended forward and inverse modeling of radiation pressure accelerations for LEO satellites." Journal of Geodesy, 94(4), 1-29.
