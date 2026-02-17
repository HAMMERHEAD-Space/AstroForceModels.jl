AstroForceModels.jl
================================

This package contains the dominant astrodynamics forces affecting the orbital trajectory of a satellite for the **SatelliteToolbox.jl** ecosystem. Currently this package implements:
- [x] Spherical Harmonics
- [x] Solar Radiation Pressure
- [x] Drag
- [x] Third Body Gravity
- [x] Relativistic
- [x] Low Thrust
- [ ] Albedo
- [ ] Solid Tides
- [ ] Spacecraft thermal emission
- [ ] Magnetic Field Effects

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroForceModels")
```

## Citing

If you use `AstroForceModels.jl` in your work, please consider citing it.

```bibtex
@software{jordan_murphy_2025_16954386,
  author       = {Jordan Murphy},
  title        = {HAMMERHEAD-Space/AstroForceModels.jl},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.16954386},
  url          = {https://doi.org/10.5281/zenodo.16954386},
}
```
