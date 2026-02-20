# Solid Body Tides

Solid body tides model the perturbation to the geopotential caused by the tidal deformation of a central body due to the gravitational attraction of nearby massive bodies. For Earth-orbiting satellites, the dominant tide-raising bodies are the Moon and Sun.

## Physical Description

When the Moon and Sun exert differential gravitational forces across the volume of the Earth, they induce small elastic deformations (tidal bulges). These mass redistributions alter the gravitational potential experienced by an orbiting satellite. The effect, while small (on the order of 10⁻¹⁰ to 10⁻⁹ km/s² in LEO), is important for high-precision orbit determination and geodetic applications.

The formulation follows IERS Conventions (2010), Section 6.2, Step 1 — using frequency-independent Love numbers and the Legendre polynomial addition theorem for efficient computation directly in the inertial frame.

### Tidal Perturbation Potential

The change in geopotential at satellite position **r** from a tide-raising body *j* at degree *n* is:

```
ΔV_j^(n) = k_n * GM_j * R_e^(2n+1) / (r_j^(n+1) * r^(n+1)) * P_n(cos γ_j)
```

Where:
- `k_n` is the degree-*n* Love number characterizing the elastic response
- `GM_j` is the gravitational parameter of body *j*
- `R_e` is the central body's equatorial radius
- `r_j` is the distance to body *j*
- `γ_j` is the geocentric angle between the satellite and body *j*
- `P_n` is the Legendre polynomial of degree *n*

### Acceleration

Taking the gradient of the potential and applying the Legendre addition theorem yields closed-form acceleration expressions for each degree, summed over all tide-raising bodies.

**Degree 2** (dominant term):
```
a_j^(2) = 3 k₂ GM_j R_e⁵ / (r_j³ r⁴) * [(1 - 5ξ²)/2 r̂ + ξ r̂_j]
```

**Degree 3** (secondary correction):
```
a_j^(3) = k₃ GM_j R_e⁷ / (2 r_j⁴ r⁵) * [(15ξ - 35ξ³) r̂ + (15ξ² - 3) r̂_j]
```

Where `ξ = r̂ · r̂_j` is the cosine of the geocentric angle.

## Components

### SolidBodyTidesModel

The main struct for solid body tides computations:

- **tide_raising_bodies**: Tuple of `ThirdBodyModel`s for each tide-raising body. Gravitational parameters are obtained from each body's `CelestialBody` definition.
- **k2**: Degree-2 Love number (default: 0.30190, IERS 2010 anelastic)
- **k3**: Degree-3 Love number (default: 0.093, IERS 2010)
- **k2_plus**: Degree-2 to degree-4 coupling Love number (default: 0.0)
- **R_e**: Central body equatorial radius [km] (default: Earth)
- **include_degree_3**: Include degree-3 contribution (default: true)

The model accepts any number of tide-raising bodies as a tuple of `ThirdBodyModel`s, enabling use with any central body — not just Earth.

## Magnitude of Effects

For Earth-orbiting satellites:

| Orbit       | Altitude   | Tidal Acceleration  |
|-------------|------------|---------------------|
| LEO         | ~400 km    | ~10⁻⁹ km/s²        |
| MEO         | ~20,000 km | ~10⁻¹¹ km/s²       |
| GEO         | ~35,786 km | ~10⁻¹² km/s²       |

The Moon contributes roughly 2–3 times the tidal acceleration of the Sun at most geometries, though the exact ratio varies with the instantaneous positions.

## Future Plans

This implementation covers IERS 2010 Step 1 (frequency-independent Love numbers, degrees 2 and 3). The full IERS solid tides model also includes frequency-dependent corrections (Step 2), permanent tide handling, ocean loading tides, and pole tides. These additional components will likely be spun out into a dedicated tidal models package with a more complete implementation. The `SolidBodyTidesModel` interface here is designed to remain stable through that transition.

## References

[1] Petit, G. and Luzum, B. (eds.), "IERS Conventions (2010)," IERS Technical Note 36, Section 6.2.

[2] Montenbruck, O. and Gill, E., "Satellite Orbits: Models, Methods, and Applications," Springer, 2000, Section 3.2.5.
