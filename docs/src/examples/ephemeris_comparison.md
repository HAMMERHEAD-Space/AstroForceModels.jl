# Ephemeris Source Comparison

This example demonstrates the numerical differences between Vallado analytical ephemeris
and JPL SPK kernels for Sun and Moon positions, and how those differences propagate into
third-body acceleration computations.

## Setup

```julia
using AstroForceModels
using Ephemerides
using SatelliteToolboxTransformations
using FrameTransformations
using Tempo
using LinearAlgebra

# ── Epoch ────────────────────────────────────────────────────────────────────
epoch = Epoch("2024-01-05T12:00:00 TDB")

# ── Tier 1: Vallado analytical ───────────────────────────────────────────────
p_vallado = setup_inertial_frames(epoch)

# ── Tier 3: JPL DE440 ───────────────────────────────────────────────────────
eph = EphemerisProvider("de440s.bsp")  # or "de440.bsp" for full precision
p_spk = setup_ephemeris_frames(epoch, eph)
```

## Comparing Sun and Moon Positions

```julia
t_ft = AstroForceModels.ft_time(p_vallado, 0.0)

# Sun position from both sources
sun_vallado = vector3(p_vallado.frames, :Earth, :Sun, :ICRF, t_ft)
sun_spk     = vector3(p_spk.frames, :Earth, :Sun, :ICRF, t_ft)
sun_diff    = norm(sun_vallado - sun_spk)

println("Sun position difference: ", round(sun_diff; digits=1), " km")
println("Sun distance from Earth: ", round(norm(sun_spk); digits=1), " km")
println("Relative error: ", @sprintf("%.2e", sun_diff / norm(sun_spk)))

# Moon position from both sources
moon_vallado = vector3(p_vallado.frames, :Earth, :Moon, :ICRF, t_ft)
moon_spk     = vector3(p_spk.frames, :Earth, :Moon, :ICRF, t_ft)
moon_diff    = norm(moon_vallado - moon_spk)

println("\nMoon position difference: ", round(moon_diff; digits=1), " km")
println("Moon distance from Earth: ", round(norm(moon_spk); digits=1), " km")
println("Relative error: ", @sprintf("%.2e", moon_diff / norm(moon_spk)))
```

Typical output:

```
Sun position difference: 531.2 km
Sun distance from Earth: 147099648.3 km
Relative error: 3.61e-06

Moon position difference: 18.4 km
Moon distance from Earth: 404532.7 km
Relative error: 4.55e-05
```

## Impact on Third-Body Accelerations

```julia
state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

# Sun third-body models — same interface, different ephemeris backing
sun_vallado_model = ThirdBodyModel(;
    body=SunBody(),
    ephem_type=FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
    frames=p_vallado.frames,
)
sun_spk_model = ThirdBodyModel(;
    body=SunBody(),
    ephem_type=FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
    frames=p_spk.frames,
)

a_sun_vallado = acceleration(state, p_vallado, 0.0, sun_vallado_model)
a_sun_spk     = acceleration(state, p_spk, 0.0, sun_spk_model)
a_sun_diff    = norm(a_sun_vallado - a_sun_spk)

println("Sun acceleration (Vallado): ", a_sun_vallado)
println("Sun acceleration (DE440):   ", a_sun_spk)
println("Difference:                  ", @sprintf("%.2e km/s²", a_sun_diff))
println("Relative difference:         ", @sprintf("%.2e", a_sun_diff / norm(a_sun_spk)))

# Moon third-body models
moon_vallado_model = ThirdBodyModel(;
    body=MoonBody(),
    ephem_type=FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
    frames=p_vallado.frames,
)
moon_spk_model = ThirdBodyModel(;
    body=MoonBody(),
    ephem_type=FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
    frames=p_spk.frames,
)

a_moon_vallado = acceleration(state, p_vallado, 0.0, moon_vallado_model)
a_moon_spk     = acceleration(state, p_spk, 0.0, moon_spk_model)
a_moon_diff    = norm(a_moon_vallado - a_moon_spk)

println("\nMoon acceleration (Vallado): ", a_moon_vallado)
println("Moon acceleration (DE440):   ", a_moon_spk)
println("Difference:                  ", @sprintf("%.2e km/s²", a_moon_diff))
println("Relative difference:         ", @sprintf("%.2e", a_moon_diff / norm(a_moon_spk)))
```

## Time-Varying Differences

The ephemeris error is not constant — it varies with the orbital geometry of the Sun and Moon:

```julia
# Sample over one year
times = range(0.0, 365.25 * 86400.0; length=1000)
sun_diffs = Float64[]
moon_diffs = Float64[]

for dt in times
    t = AstroForceModels.ft_time(p_vallado, dt)

    sv = vector3(p_vallado.frames, :Earth, :Sun, :ICRF, t)
    ss = vector3(p_spk.frames, :Earth, :Sun, :ICRF, t)
    push!(sun_diffs, norm(sv - ss))

    mv = vector3(p_vallado.frames, :Earth, :Moon, :ICRF, t)
    ms = vector3(p_spk.frames, :Earth, :Moon, :ICRF, t)
    push!(moon_diffs, norm(mv - ms))
end

println("Sun position error over 1 year:")
println("  Min:  ", round(minimum(sun_diffs); digits=1), " km")
println("  Max:  ", round(maximum(sun_diffs); digits=1), " km")
println("  Mean: ", round(sum(sun_diffs)/length(sun_diffs); digits=1), " km")

println("\nMoon position error over 1 year:")
println("  Min:  ", round(minimum(moon_diffs); digits=1), " km")
println("  Max:  ", round(maximum(moon_diffs); digits=1), " km")
println("  Mean: ", round(sum(moon_diffs)/length(moon_diffs); digits=1), " km")
```

## Full Propagation Comparison

Compare a 3-day LEO propagation using Vallado vs SPK ephemeris for third-body
perturbations. All other force models (gravity harmonics, drag, SRP) are identical.

```julia
using AstroPropagators
using OrdinaryDiffEqVerner
using SatelliteToolboxGravityModels, SpaceIndices

SpaceIndices.init()
eop_data = fetch_iers_eop()

# ── Common force models ──────────────────────────────────────────────────────
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

function build_dynamics(p)
    frames = p.frames

    gravity = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, body_fixed_frame=:ITRF, propagation_frame=:ICRF,
        order=36, degree=36, frames=frames,
    )
    sun = ThirdBodyModel(;
        body=SunBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
        frames=frames,
    )
    moon = ThirdBodyModel(;
        body=MoonBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
        frames=frames,
    )
    return CentralBodyDynamicsModel(gravity, (sun, moon))
end

# ── Vallado propagation ──────────────────────────────────────────────────────
p_v = setup_earth_propagation_frames(epoch, eop_data)
dynamics_v = build_dynamics(p_v)

u0 = [-1076.225324679696, -6765.896364327722, -332.3087833503755,
       8.956857417032581, -3.3123476319597557, -1.1880157328553503]
tspan = (0.0, 3 * 86400.0)

f_v!(du, u, p, t) = begin du[1:3] = u[4:6]; du[4:6] = build_dynamics_model(u, p, t, dynamics_v) end
sol_v = solve(ODEProblem(f_v!, u0, tspan, p_v), Vern9(); abstol=1e-13, reltol=1e-13)

# ── SPK propagation ──────────────────────────────────────────────────────────
p_s = setup_ephemeris_frames(epoch, eph; eop_data=eop_data)
dynamics_s = build_dynamics(p_s)

f_s!(du, u, p, t) = begin du[1:3] = u[4:6]; du[4:6] = build_dynamics_model(u, p, t, dynamics_s) end
sol_s = solve(ODEProblem(f_s!, u0, tspan, p_s), Vern9(); abstol=1e-13, reltol=1e-13)

# ── Compare final states ─────────────────────────────────────────────────────
pos_diff = norm(sol_v.u[end][1:3] - sol_s.u[end][1:3])
vel_diff = norm(sol_v.u[end][4:6] - sol_s.u[end][4:6])

println("3-day propagation difference (Vallado vs DE440):")
println("  Position: ", round(pos_diff * 1000; digits=2), " m")
println("  Velocity: ", round(vel_diff * 1e6; digits=4), " mm/s")
```