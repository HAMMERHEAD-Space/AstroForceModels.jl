# Geomagnetic Lorentz Force

The geomagnetic Lorentz force model computes the perturbative acceleration experienced by a charged spacecraft moving through Earth's magnetic field. When a spacecraft accumulates electrostatic charge (from the space plasma environment or intentional charging), it interacts with the geomagnetic field to produce a velocity-dependent force perpendicular to both the velocity and magnetic field vectors.

## Physical Description

A spacecraft in Earth orbit naturally acquires surface charge through interactions with the ambient plasma environment (photoelectric emission, electron/ion collection). The resulting charge creates a Lorentz force when the spacecraft moves relative to the co-rotating geomagnetic field:

```
𝐅 = q(𝐯_rel × 𝐁)
```

where:
- `q` is the spacecraft charge [C]
- `𝐯_rel = 𝐯 - 𝛚×𝐫` is the velocity relative to the co-rotating magnetic field [m/s]
- `𝐁` is the geomagnetic field vector [T]
- `𝛚` is Earth's angular velocity vector

The acceleration is:

```
𝐚 = (q/m)(𝐯_rel × 𝐁)
```

### Key Properties

- **Perpendicular force**: Always acts perpendicular to both the relative velocity and the magnetic field
- **Non-dissipative**: Does not add or remove energy in the rotating frame
- **Charge-dependent**: Linear in the charge-to-mass ratio, enabling propellantless control via active charging
- **Strongest at poles**: The Lorentz force effect is most pronounced in polar orbits where the magnetic field is strongest

## Geomagnetic Field Models

Two geomagnetic field models are available:

### IGRF (International Geomagnetic Reference Field)

The IGRF v14 model uses a spherical harmonic expansion up to degree 13 for high-fidelity geomagnetic field computation. Valid for dates between 1900 and 2035. This is the recommended model for quantitative analysis.

### Simplified Dipole

A simplified dipole model that assumes Earth's magnetic field is a perfect tilted dipole. Less accurate than IGRF but computationally cheaper and sufficient for preliminary analysis where uncertainties are high.

## Components

### Spacecraft Charge Models

Two charge-to-mass ratio models are provided:

- **`FixedChargeMassRatio(q_over_m)`**: Fixed charge-to-mass ratio [C/kg]. Use for constant-charge analysis.
- **`StateChargeModel(f)`**: State-dependent model where `f(u, p, t)` returns q/m [C/kg]. Use for time-varying or active charge control scenarios.

### MagneticFieldAstroModel

The main force model struct:

- **`spacecraft_charge_model`**: Model providing the charge-to-mass ratio q/m [C/kg]
- **`geomagnetic_field_model`**: `IGRFField()` (default: `DipoleMagneticField()`) for field computation
- **`eop_data`**: Earth Orientation Parameters for coordinate transformations
- **`max_degree`**: Maximum spherical harmonic degree for IGRF (default: 13, ignored for dipole)
- **`P`**, **`dP`**: Optional pre-allocated Legendre polynomial buffers for allocation-free IGRF evaluation

## Typical Charge-to-Mass Ratios

| Scenario | q/m [C/kg] | Notes |
|----------|-----------|-------|
| Natural LEO charging | 1e-6 to 1e-4 | Passive surface charging from plasma environment |
| GEO charging events | 1e-5 to 1e-3 | Geomagnetic storm surface charging |
| Active Lorentz propulsion | ~0.01 to 0.03 | Proposed intentional charging for propellantless maneuvers |

## Magnitude of Effects

For a spacecraft with q/m = 1e-3 C/kg in LEO (B ~ 30 μT, v_rel ~ 7.6 km/s):

```
|a| ~ (q/m) × v_rel × B ~ 1e-3 × 7600 × 3e-5 ≈ 2.3e-4 m/s² ≈ 2.3e-7 km/s²
```

This is a second-order perturbation for typical natural charging levels, but becomes significant for proposed active Lorentz-augmented orbit applications.

## References

[1] Peck, M. A. (2005). "Prospects and Challenges for Lorentz-Augmented Orbits." AIAA Guidance, Navigation, and Control Conference. AIAA 2005-5995.

[2] Streetman, B. & Peck, M. A. (2007). "New Synchronous Orbits Using the Geomagnetic Lorentz Force." Journal of Guidance, Control, and Dynamics, 30(6), 1677-1690.

[3] Khalil, K. I. & Abdel-Aziz, Y. A. (2014). "Electromagnetic effects on the orbital motion of a charged spacecraft." Research in Astronomy and Astrophysics, 14(5), 589.
