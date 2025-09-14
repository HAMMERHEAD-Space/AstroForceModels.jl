# AstroForceModels.jl

[![CI](https://github.com/HAMMERHEAD-Space/AstroForceModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/AstroForceModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/HAMMERHEAD-Space/AstroForceModels.jl/graph/badge.svg?token=OYYBK1VZ6C)](https://codecov.io/gh/HAMMERHEAD-Space/AstroForceModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/762543633.svg)](https://doi.org/10.5281/zenodo.16954385)

This package contains the dominant astrodynamics forces affecting the orbital trajectory of a satellite. Currently this package implements:
- [x] Zonal Harmonics
- [x] Solar Radiation Pressure
- [x] Drag
- [x] Third Body Gravity
- [x] Relativistic
- [ ] Albedo
- [ ] Solid Tides
- [ ] Spacecraft thermal emission
- [ ] Magnetic Field Effects

The plan is likely to merge this effort into Julia Space Mission Design's AstroModels -- https://github.com/JuliaSpaceMissionDesign/AstroModels.jl, but this is still under discussion. In the meantime interfaces to the models there will be added.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroForceModels")
```

## Documentation

For more information, see the [documentation][docs-dev-url].

## Citing

If you use `AstroForceModels.jl` in your work, please consider citing it.

```bibtex
@software{jordan_murphy_2025_16954386,
  author       = {Jordan Murphy},
  title        = {HAMMERHEAD-Space/AstroForceModels.jl: v0.3.13},
  month        = aug,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v0.3.13},
  doi          = {10.5281/zenodo.16954386},
  url          = {https://doi.org/10.5281/zenodo.16954386},
}
```

[docs-dev-url]: https://HAMMERHEAD-Space.github.io/AstroForceModels.jl/stable/
[docs-stable-url]: https://HAMMERHEAD-Space.github.io/AstroForceModels.jl/stable/
