# Third Body Gravity

Third body gravity accounts for the gravitational perturbations on a spacecraft caused by celestial bodies other than the primary body (typically Earth). The most significant third bodies are the Sun and Moon, but other planets can also contribute to orbital perturbations, especially for high-altitude satellites or interplanetary missions.

## Physical Description

In the two-body problem, only the gravitational attraction between the spacecraft and the primary body (Earth) is considered. However, real spacecraft are also influenced by the gravitational fields of other massive bodies in the solar system. These third body perturbations become increasingly important as the spacecraft's distance from Earth increases.

The third body acceleration is given by:

```
a_third = μ_third * [(r_third - r_sat)/|r_third - r_sat|³ - r_third/|r_third|³]
```

Where:
- `μ_third` is the gravitational parameter of the third body
- `r_third` is the position vector of the third body relative to Earth's center
- `r_sat` is the position vector of the spacecraft relative to Earth's center
- The first term represents direct attraction to the third body
- The second term represents the differential effect (Earth's acceleration toward the third body)

## Key Third Bodies

### Sun
- **Gravitational parameter**: μ☉ = 1.327×10¹¹ km³/s²
- **Distance from Earth**: ~1 AU (varies ±3.4% due to eccentricity)
- **Primary effects**: Long-period orbital variations, secular drift
- **Most significant for**: High-altitude and interplanetary missions

### Moon  
- **Gravitational parameter**: μ🌙 = 4.903×10³ km³/s²
- **Distance from Earth**: ~384,400 km (varies ±11% due to eccentricity)
- **Primary effects**: Short-period oscillations, resonance phenomena
- **Most significant for**: High Earth orbits, lunar missions

### Other Planets
For specialized missions, other planets may be significant:
- **Jupiter**: μ♃ = 1.267×10⁸ km³/s² (interplanetary missions)
- **Venus**: μ♀ = 3.249×10⁵ km³/s² (inner solar system missions)
- **Mars**: μ♂ = 4.283×10⁴ km³/s² (Mars missions and transfers)

## Components

The third body force model in AstroForceModels includes:

### ThirdBodyModel

The main struct for third body computations and ephemeris provider:

- **body**: `CelestialBody` instance (e.g., `SunBody()`, `MoonBody()`)
- **ephem_type**: `FrameEphemeris` specifying how to compute body position via FrameTransformations.jl (center point, target point, and axes)

### CelestialBody

Defines properties of celestial bodies using a `Name` type parameter for compile-time dispatch:

- **Name** (type parameter): Body identifier symbol (e.g., `:Sun`, `:Moon`, `:Earth`)
- **central_body**: Central body symbol
- **jpl_code**: NAIF ID code
- **μ**: Gravitational parameter [km³/s²]
- **Req**: Equatorial radius [km]

Use `nameof(body)` to retrieve the body name. Pre-built constructors: `SunBody()`, `MoonBody()`, `EarthBody()`.

### FrameEphemeris

Ephemeris type that uses FrameTransformations.jl to compute body positions and velocities via SPK kernels:

- **center_point**: NAIF ID of the center body (e.g., 399 for Earth, 2000433 for Eros)
- **target_point**: NAIF ID of the target body (e.g., 10 for Sun, 301 for Moon)
- **axes**: Frame for output (default: `:ICRF`)

## Usage Example

```julia
using AstroForceModels
using Tempo, ComponentArrays

# Create third body models using FrameEphemeris
sun_model = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=10, axes=:ICRF)
)
moon_model = ThirdBodyModel(
    body = MoonBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=301, axes=:ICRF)
)

# Set up frame-aware parameters (see examples for full frame setup)
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = FrameAwareParams(frames, epoch, :ICRF)

# These can be combined with other force models
acceleration(state, p, time, sun_model)
acceleration(state, p, time, moon_model)
```

## Implementation Details

Third body acceleration computation involves:

1. **Ephemeris Calculation**: Computing third body position at current time
2. **Coordinate Transformation**: Converting to appropriate reference frame
3. **Distance Calculation**: Computing distances from spacecraft and Earth to third body
4. **Differential Acceleration**: Computing the difference in gravitational acceleration
5. **Vector Operations**: Determining acceleration direction and magnitude

## Magnitude of Effects

Third body perturbations vary significantly with orbital characteristics:

### Low Earth Orbit (LEO, ~400 km)
- **Sun**: ~10⁻⁹ m/s² (negligible)
- **Moon**: ~10⁻⁹ m/s² (negligible)
- **Dominates**: Atmospheric drag, Earth oblateness

### Medium Earth Orbit (MEO, ~20,000 km) 
- **Sun**: ~10⁻⁷ m/s² (significant)
- **Moon**: ~10⁻⁶ m/s² (very significant)
- **Applications**: GPS constellation, navigation satellites

### Geostationary Orbit (GEO, ~35,786 km)
- **Sun**: ~10⁻⁶ m/s² (very significant)  
- **Moon**: ~10⁻⁶ m/s² (very significant)
- **Effects**: Station-keeping requirements, orbital drift

### High Earth Orbit (HEO, >35,786 km)
- **Sun**: ~10⁻⁶ to 10⁻⁵ m/s² (dominant perturbation)
- **Moon**: ~10⁻⁶ to 10⁻⁵ m/s² (dominant perturbation)
- **Applications**: Scientific missions, space telescopes

## References

[1] Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[2] Montenbruck, Oliver, and Eberhard Gill. "Satellite Orbits: Models, Methods and Applications." Springer, 2000.
[3] Battin, Richard H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA, 1999.
[4] Kaplan, Marshall H. "Modern Spacecraft Dynamics and Control." Wiley, 1976.
[5] Cappellari, J. O., et al. "Mathematical Theory of the Goddard Trajectory Determination System." NASA, 1976. 