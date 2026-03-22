# Small Body (Eros) Orbiter

Spacecraft orbiting asteroid 433 Eros with Keplerian gravity, a constant-rate
body-fixed rotating frame, Sun third-body perturbation, and SRP with Eros as the
occulting body.

This example is fully self-contained and can be run as-is. It uses Keplerian
(point-mass) gravity. For higher fidelity, load a shape-based gravity model
(e.g., from NEAR Shoemaker data) via the SmallBodyGravity package and replace
the `KeplerianGravityAstroModel` with a `GravityHarmonicsAstroModel`.

```julia
using AstroForceModels
using ComponentArrays
using FrameTransformations
using Tempo
using StaticArrays

# ── Epoch ────────────────────────────────────────────────────────────────────
JD = date_to_jd(2035, 3, 1, 0, 0, 0.0)

# ── Eros physical constants ─────────────────────────────────────────────────
μ_eros      = 4.463e-4          # km³/s² (GM of Eros)
R_eros      = 8.42              # km (mean radius)
eros_period = 5.27 * 3600.0     # s (rotation period ≈ 5.27 hours)

# Pole orientation in ICRF (RA=11.35°, Dec=17.22° from IAU 2000 report)
eros_pole_ra  = deg2rad(11.35)
eros_pole_dec = deg2rad(17.22)
eros_axis = SVector{3}(
    cos(eros_pole_dec) * cos(eros_pole_ra),
    cos(eros_pole_dec) * sin(eros_pole_ra),
    sin(eros_pole_dec),
)

# ── Frame System ─────────────────────────────────────────────────────────────
frames = FrameSystem{2, Float64}()
add_axes_icrf!(frames)

# Eros body-fixed rotating frame (direct child of ICRF)
add_small_body_rotating_frame!(
    frames, :ErosFixed, 2000433, 1,
    eros_axis,
    eros_period,
)

# Eros as root point, Sun as dynamical child
add_point!(frames, :Eros, 2000433, :ICRF)

# Approximate Sun position relative to Eros via Vallado Earth-Sun vector.
# In production, use SPK kernels: add_point_ephemeris!(frames, eph, :Sun, 10, 2000433)
function sun_from_eros_state(t)
    # Scale Vallado Earth-Sun to approximate Eros-Sun distance (~1.46 AU)
    earth_sun = vallado_sun_state(t)
    return vcat(earth_sun[1:3] .* 1.46, earth_sun[4:6] .* 1.46)
end
add_point_dynamical!(frames, :Sun, 10, 2000433, :ICRF, sun_from_eros_state)

# ── Parameters ───────────────────────────────────────────────────────────────
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = FrameAwareParams(frames, epoch, :ICRF)

# ── Force Models ─────────────────────────────────────────────────────────────

# 1. Eros gravity (Keplerian point-mass)
eros_gravity = KeplerianGravityAstroModel(μ=μ_eros)

# 2. Sun third-body perturbation relative to Eros
sun_tb = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=2000433, target_point=10, axes=:ICRF),
    frames = frames,
)

# 3. SRP (Eros as occulting body)
srp = SRPAstroModel(;
    satellite_srp_model = CannonballFixedSRP(0.05),
    sun_data     = sun_tb,
    shadow_model = Conical(),
    R_Occulting  = R_eros,
)

# ── Combined Dynamics ────────────────────────────────────────────────────────
dynamics = CentralBodyDynamicsModel(
    eros_gravity,
    (sun_tb, srp),
)

# ── Compute acceleration at close retrograde orbit (35 km altitude) ──────────
v_circ = sqrt(μ_eros / (R_eros + 35.0))
u0 = [R_eros + 35.0, 0.0, 0.0, 0.0, 0.0, -v_circ]  # km, km/s (retrograde)
total_accel = build_dynamics_model(u0, p, 0.0, dynamics)
println("Total acceleration: ", total_accel, " km/s²")
```
