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

### ThirdBodyAstroModel

The main struct for third body computations:

- **third_body_model**: Contains celestial body properties and ephemeris data
- **eop_data**: Earth orientation parameters for coordinate transformations

### ThirdBodyModel

Contains information about specific celestial bodies:

- **body**: Celestial body identifier (SUN, MOON, etc.)
- **ephemeris_type**: Method for computing body position (:analytical, :de430, etc.)
- **gravitational_parameter**: μ value for the celestial body

### CelestialBody

Defines properties of celestial bodies:

- **name**: Body identifier  
- **μ**: Gravitational parameter
- **radius**: Physical radius
- **ephemeris_data**: Orbital elements or ephemeris coefficients

## Usage Example

```julia
using AstroForceModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxBase

# Load Earth orientation parameters
eop_data = fetch_iers_eop()

# Create Sun model
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
# Create Moon model  
moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

# Create third body models
sun_perturbation = ThirdBodyAstroModel(
    third_body_model = sun_model,
    eop_data = eop_data
)

moon_perturbation = ThirdBodyAstroModel(
    third_body_model = moon_model, 
    eop_data = eop_data
)

# These can be combined with other force models
acceleration(state, parameters, time, sun_perturbation)
acceleration(state, parameters, time, moon_perturbation)
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