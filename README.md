# AstroForceModels.jl

[![CI](https://github.com/jmurphy6895/AstroForceModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroForceModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/AstroForceModels.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/AstroForceModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package contains the dominant astrodynamics forces affecting the orbital trajectory of a satellite. Currently this package implements:
- [x] Zonal Harmonics
- [x] Solar Radiation Pressure
- [x] Drag
- [x] Third Body Gravity
- [x] Relativistic
- [ ] Albedo
- [ ] Solid Tides

This first implementation is just to get a set of working models, the next steps are to get everything allocation-free and working with automatic differentiation. Furthermore, the plan is likely to merge this effort into Julia Space Mission Design's AstroModels -- https://github.com/JuliaSpaceMissionDesign/AstroModels.jl. But plans are still under disccusion there.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroForceModels")
```

## Documentation

For more information, see the [documentation][docs-dev-url].

[docs-dev-url]: https://jmurphy6895.github.io/AstroForceModels.jl/stable/
[docs-stable-url]: https://jmurphy6895.github.io/AstroForceModels.jl/stable/
