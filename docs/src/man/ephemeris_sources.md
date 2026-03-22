# Ephemeris Sources

AstroForceModels supports three tiers of ephemeris data for celestial body positions and
Earth orientation. Each tier trades off simplicity against fidelity. All three are
interchangeable — force models do not care which source populated the `FrameSystem`.

## Overview

| Tier | Setup Function | Ephemeris Source | ITRF Source | Typical Accuracy |
|------|---------------|-----------------|-------------|-----------------|
| **Analytical** | `setup_inertial_frames` | Vallado analytical (Sun/Moon) | — | ~100–500 km Sun, ~10–50 km Moon |
| **Earth HF** | `setup_earth_propagation_frames` | Vallado analytical (Sun/Moon) | SatelliteToolboxTransformations EOP | Same as above + accurate Earth rotation |
| **JPL Kernels** | `setup_ephemeris_frames` | SPK binary kernels (DE440, etc.) | SatelliteToolboxTransformations EOP or IERSConventions | Sub-km (limited by kernel precision) |

## Tier 1: Vallado Analytical (`setup_inertial_frames`)

Uses closed-form expressions from Vallado's *Fundamentals of Astrodynamics and
Applications* for the Sun and Moon positions. The Sun model uses the low-precision solar
ephemeris in the MOD frame, rotated to J2000. The Moon model uses Vallado's analytical
lunar ephemeris with a finite-difference velocity.

**Pros:**
- Zero external data files — works out of the box
- Fast to initialize (no file I/O)
- Sufficient for most LEO/MEO mission analysis

**Cons:**
- Sun position accuracy ~100–500 km (varies with epoch)
- Moon position accuracy ~10–50 km
- No planetary positions beyond Sun and Moon
- Moon velocity is finite-differenced (not analytical)

**When to use:**
- Quick analyses, trade studies, algorithm development
- LEO/MEO propagation where third-body perturbations are secondary
- Anywhere the ~0.001% position error in Sun/Moon doesn't matter

```julia
using AstroForceModels, Tempo

epoch = Epoch("2024-01-05T12:00:00 TDB")
p = setup_inertial_frames(epoch)
```

## Tier 2: Vallado + ITRF (`setup_earth_propagation_frames`)

Same Sun/Moon ephemeris as Tier 1, but adds the ITRF (International Terrestrial Reference
Frame) for Earth body-fixed computations. The ITRF rotation uses Earth Orientation
Parameters from SatelliteToolboxTransformations.

**Pros:**
- Accurate Earth rotation (EOP-based precession, nutation, polar motion)
- Enables gravity harmonics, atmospheric drag, albedo, magnetic field models
- One-line setup for the full Earth propagation environment

**Cons:**
- Requires `fetch_iers_eop()` call (downloads EOP data on first use)
- Same Sun/Moon accuracy limitations as Tier 1
- Earth-specific — not useful for Mars/asteroid missions

**When to use:**
- Any Earth-centric propagation requiring body-fixed frames
- High-fidelity LEO with gravity harmonics and drag
- The recommended default for most Earth satellite work

```julia
using AstroForceModels, SatelliteToolboxTransformations, Tempo

epoch = Epoch("2024-01-05T12:00:00 TDB")
eop_data = fetch_iers_eop()
p = setup_earth_propagation_frames(epoch, eop_data)
```

## Tier 3: JPL SPK Kernels (`setup_ephemeris_frames`)

Uses JPL binary SPK (Spacecraft and Planet Kernel) files for body positions. These are
the same ephemeris files used by NASA/JPL mission design tools (SPICE). Body positions are
computed via Chebyshev polynomial interpolation of numerically integrated trajectories.

Requires [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) to
be loaded (`using Ephemerides`).

**Pros:**
- Sub-km accuracy for all solar system bodies
- Consistent velocities (no finite differencing)
- Full solar system coverage (planets, barycenters, Moon, Pluto)
- Can load mission-specific kernels (asteroid, comet, spacecraft)
- Used by NASA/ESA for mission design — traceable to the same data

**Cons:**
- Requires SPK kernel files (e.g., `de440.bsp`, ~120 MB for DE440)
- Slightly slower initialization (file parsing)
- Kernels have finite time spans (DE440 covers 1550–2650)

**When to use:**
- High-fidelity propagation matching flight dynamics standards
- Interplanetary missions (Mars, asteroids, etc.)
- Validation against JPL/GSFC tools
- Any analysis where sub-km ephemeris accuracy matters

```julia
using AstroForceModels, Ephemerides, Tempo

# SPK only — planetary positions
eph = EphemerisProvider("de440.bsp")
epoch = Epoch("2024-01-05T12:00:00 TDB")
p = setup_ephemeris_frames(epoch, eph)
```

With ITRF via EOP:

```julia
using SatelliteToolboxTransformations

eop_data = fetch_iers_eop()
p = setup_ephemeris_frames(epoch, eph; eop_data=eop_data)
```

With PCK body orientation kernels:

```julia
# SPK + Earth high-precision PCK — both positions and Earth rotation from kernels
eph = EphemerisProvider(["de440.bsp", "earth_latest_high_prec.bpc"])
p = setup_ephemeris_frames(epoch, eph)
# p.frames now has :EarthPCK axes from the PCK kernel

# SPK + Moon PA kernel — for lunar missions
eph = EphemerisProvider(["de440.bsp", "moon_pa_de440.bpc"])
p = setup_ephemeris_frames(epoch, eph)
# p.frames now has :MoonPA440 axes
```

### PCK (Planetary Constants Kernel) Support

Binary PCK files (`.bpc`) provide high-precision body orientation data — the rotation
from an inertial frame to a body-fixed frame as a function of time. When loaded alongside
SPK files, `setup_ephemeris_frames` automatically registers the orientation axes.

Common PCK files:

| Kernel | Axes ID | Axes Name | Description |
|--------|---------|-----------|-------------|
| `earth_latest_high_prec.bpc` | 3000 | `:EarthPCK` | Earth body-fixed (high-precision, relative to ECL2000) |
| `moon_pa_de440.bpc` | 31008 | `:MoonPA440` | Moon Principal Axes DE440 |
| `moon_pa_de421.bpc` | 31007 | `:MoonPA421` | Moon Principal Axes DE421 |

The PCK axes are registered under their NAIF axes IDs. If the kernel contains axes IDs
not in the default book, they are automatically named `:Axes_NNNNN`.

!!! note
    Earth PCK files express rotation relative to the Ecliptic J2000 frame (NAIF ID 17).
    `setup_ephemeris_frames` automatically registers ICRF, GCRF, EME2000, and ECL2000 as
    standard inertial axes so that PCK parent chains resolve correctly.

## Obtaining SPK/PCK Kernels

JPL distributes ephemeris kernels at:

- **Planetary SPK (DE440):** [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/)
- **Satellite SPK:** [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/)
- **Binary PCK (body orientation):** [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/)
- **Small body kernels:** Generated via [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/)

Common kernels:

| Kernel | Size | Coverage | Contents |
|--------|------|----------|----------|
| `de440s.bsp` | ~3 MB | 1849–2150 | Short-span planetary ephemeris |
| `de440.bsp` | ~120 MB | 1550–2650 | Full planetary ephemeris |
| `de441.bsp` | ~3.1 GB | −13200–17191 | Extended planetary ephemeris |

## NAIF ID Reference

The default body name book used by `setup_ephemeris_frames`:

| NAIF ID | Name | Description |
|---------|------|-------------|
| 0 | `:SSB` | Solar System Barycenter |
| 1–9 | `:MercuryBarycenter` ... `:PlutoBarycenter` | Planet-system barycenters |
| 10 | `:Sun` | Sun |
| 199 | `:Mercury` | Mercury |
| 299 | `:Venus` | Venus |
| 301 | `:Moon` | Moon |
| 399 | `:Earth` | Earth |
| 499 | `:Mars` | Mars |

You can provide a custom `book` dictionary for non-standard bodies:

```julia
custom_book = Dict{Int,Symbol}(
    0 => :SSB, 10 => :Sun, 399 => :Earth,
    2000433 => :Eros,  # Custom asteroid
)
p = setup_ephemeris_frames(epoch, eph; book=custom_book)
```

## Mixing Sources

The three tiers are not mutually exclusive. You can start with Tier 1 for prototyping
and switch to Tier 3 for validation without changing any force model code:

```julia
# Force models are identical regardless of ephemeris source
sun = ThirdBodyModel(;
    body=SunBody(),
    ephem_type=FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
    frames=p.frames,  # ← only this changes between tiers
)
```

The `FrameEphemeris` type queries positions from the `FrameSystem` graph using NAIF IDs.
It does not care whether those positions came from Vallado formulas or JPL kernels — the
force model sees the same interface.

## Earth Orientation: SatelliteToolboxTransformations vs IERSConventions

Two options exist for ITRF (Earth rotation):

| Feature | SatelliteToolboxTransformations | IERSConventions (via FrameTransformations) |
|---------|-------------------------------|------------------------------------------|
| Setup | `fetch_iers_eop()` → `eop_data` | `eop_load_data!(iers2010b, file)` |
| Model | FK5/IAU-76/80 | IERS 2010 (CIO-based) |
| Derivative | Finite-difference (0.01 s step) | Analytical |
| Convenience | `setup_earth_propagation_frames(epoch, eop_data)` | `setup_ephemeris_frames(epoch, eph; include_itrf=true)` |
| External files | Auto-downloaded on first use | Requires manual EOP data setup |

For most applications, SatelliteToolboxTransformations is simpler. IERSConventions is
more rigorous (IERS 2010 conventions with analytical derivatives) but requires more setup.
