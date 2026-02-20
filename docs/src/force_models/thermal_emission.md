# Thermal Emission

Spacecraft thermal emission (thermal re-radiation) models the perturbative force arising from anisotropic thermal photon emission from spacecraft surfaces. When a spacecraft absorbs solar radiation, the surface heats up and re-emits thermal infrared radiation. If the emission is not isotropic — for example, because solar panels have different thermal emissivities on front and back surfaces — the net photon recoil creates a small but persistent acceleration.

## Physical Description

The dominant source of imbalanced thermal radiation on typical spacecraft is the solar panels, due to their large exposed area and low heat capacity. The front (Sun-facing) side absorbs solar radiation and re-emits thermally, while the back side also emits thermally but at a potentially different rate governed by its emissivity. The recoil force from this asymmetric emission acts along the Sun-spacecraft line.

For GNSS satellites, thermal emission forces can perturb orbits by 1-3 cm, particularly during eclipse season transitions when spacecraft undergo periodic thermal cycling. The effect has been shown to be a limiting factor in precise orbit determination at the centimeter level.

### Equilibrium Temperature Model

In the steady-state (equilibrium) approximation, the surface temperature is determined by balancing absorbed solar radiation against thermal emission from both sides:

```
(1 - η) · α · S₀ = (ε_f + ξ · ε_b) · σ · T_eq⁴
```

where:
- `η` is the electric conversion efficiency of solar panels
- `α` is the surface absorptivity
- `S₀` is the solar irradiance at 1 AU
- `ε_f` and `ε_b` are the front and back emissivities
- `ξ` is the temperature ratio factor between front and back sides
- `σ` is the Stefan-Boltzmann constant

### Thermal Emission Acceleration

The net thermal acceleration from the emissivity imbalance is:

```
a_thm = ν · (2/3) · (A / (M · c)) · (ε_f - ξ · ε_b) · σ · T_eq⁴ · d̂_sun
```

Substituting the equilibrium temperature, this simplifies to:

```
a_thm = ν · C_thm · Ψ · (AU / d_sun)² · d̂_sun
```

where the **thermal emission coefficient** is:

```
C_thm = (2/3) · (A/M) · [(ε_f - ξ·ε_b) / (ε_f + ξ·ε_b)] · (1 - η) · α
```

and:
- `ν` is the shadow factor (0 in full shadow, 1 in full sunlight)
- `Ψ` is the solar radiation pressure at 1 AU [N/m²]
- `d_sun` is the Sun-spacecraft distance
- `d̂_sun` is the unit vector from Sun to spacecraft

Note that when `ε_f > ξ·ε_b`, the coefficient is positive and the force is directed away from the Sun; when `ε_f < ξ·ε_b`, the force is toward the Sun.

## Components

### Satellite Thermal Models

Two satellite thermal model types are provided:

- **`FixedThermalEmission(C_thm)`**: Directly specify the thermal emission coefficient [m²/kg].
- **`FlatPlateThermalModel`**: Compute the coefficient from physical surface properties (area, mass, emissivities, absorptivity, efficiency).

### ThermalEmissionAstroModel

The main force model struct:

- **`satellite_thermal_model`**: Satellite thermal model providing the emission coefficient
- **`sun_data`**: `ThirdBodyModel` for computing the Sun's position
- **`eop_data`**: Earth Orientation Parameters
- **`shadow_model`**: Shadow model type (default: `Conical()`)
- **`R_Sun`**, **`R_Occulting`**, **`Ψ`**, **`AU`**: Physical constants with sensible defaults

## Magnitude of Effects

| Satellite Type | Thermal Acceleration | Orbit Impact |
|----------------|---------------------|--------------|
| GPS Block II   | ~1.0 nm/s²          | >10 m / 7 days |
| Galileo IOV    | ~1-3 nm/s²          | ~2-3 cm (eclipse season) |
| Galileo FOC    | ~1-2 nm/s²          | ~1-2 cm (eclipse season) |

## Limitations

This implementation uses the **steady-state equilibrium** temperature model, which accurately captures the constant thermal acceleration in sunlight. The transient thermal response during eclipse transitions (cooling and heating phases), which requires solving the thermal ODE, is not currently modeled. For most applications outside of eclipse-season GNSS orbit determination at sub-centimeter precision, the equilibrium model is sufficient.

## References

[1] Duan, B. & Hugentobler, U. (2022). "Estimating surface optical properties and thermal thrust for Galileo satellite body and solar panels." GPS Solutions, 26, 135.

[2] Vigue, Y., Schutz, B. E. & Abusali, P. (1994). "Thermal force modeling for global positioning system using the finite element method." Journal of Spacecraft and Rockets, 31(5), 855-859.

[3] Milani, A., Nobili, A. M. & Farinella, P. (1987). "Non-gravitational Perturbations and Satellite Geodesy." Adam Hilger, Bristol.
