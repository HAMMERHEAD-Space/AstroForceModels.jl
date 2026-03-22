# Mars Orbiter

Mars orbiter with Keplerian gravity, Sun third-body perturbation, and SRP.
Uses a constant-rate body-fixed rotating frame based on Mars IAU pole and spin rate.

This example is fully self-contained and can be run as-is. It uses Keplerian
(point-mass) gravity. For higher fidelity, load a Mars gravity model (e.g., GMM-3,
MRO120D) in ICGEM format via `GravityModels.load(IcgemFile, "path/to/file.gfc")` and
replace the `KeplerianGravityAstroModel` with a `GravityHarmonicsAstroModel`.

```julia
using AstroForceModels
using ComponentArrays
using FrameTransformations
using Tempo
using StaticArrays

# ── Epoch ────────────────────────────────────────────────────────────────────
JD = date_to_jd(2028, 6, 15, 0, 0, 0.0)

# ── Mars physical constants ──────────────────────────────────────────────────
μ_mars     = 42828.375816   # km³/s²
R_mars     = 3396.19        # km (equatorial radius)
mars_period = 24.6229 * 3600.0  # s (sidereal rotation period)

# Mars pole orientation in ICRF (IAU 2015 report)
mars_pole_ra  = deg2rad(317.68143)
mars_pole_dec = deg2rad(52.88650)
mars_axis = SVector{3}(
    cos(mars_pole_dec) * cos(mars_pole_ra),
    cos(mars_pole_dec) * sin(mars_pole_ra),
    sin(mars_pole_dec),
)

# ── Frame System ─────────────────────────────────────────────────────────────
frames = FrameSystem{2, Float64}()
add_axes_icrf!(frames)

# Mars as root point (NAIF ID 499)
add_point!(frames, :Mars, 499, :ICRF)

# Mars body-fixed rotating frame (direct child of ICRF)
add_small_body_rotating_frame!(
    frames, :IAU_MARS, 10499, 1,
    mars_axis,
    mars_period,
)

# Approximate Sun position relative to Mars via Vallado Earth-Sun vector.
# In production, use SPK kernels: add_point_ephemeris!(frames, eph, :Sun, 10, 499)
function sun_from_mars_state(t)
    # Negate and scale Vallado Earth-Sun to approximate Mars-Sun
    earth_sun = vallado_sun_state(t)
    return vcat(-earth_sun[1:3] .* 1.524, -earth_sun[4:6] .* 1.524)
end
add_point_dynamical!(frames, :Sun, 10, 499, :ICRF, sun_from_mars_state)

# ── Parameters ───────────────────────────────────────────────────────────────
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = FrameAwareParams(frames, epoch, :ICRF)

# ── Force Models ─────────────────────────────────────────────────────────────

# 1. Mars gravity (Keplerian point-mass)
mars_gravity = KeplerianGravityAstroModel(μ=μ_mars)

# 2. Sun third-body perturbation relative to Mars
sun_tb = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=499, target_point=10, axes=:ICRF),
    frames = frames,
)

# 3. SRP (Mars as occulting body)
srp = SRPAstroModel(;
    satellite_srp_model = CannonballFixedSRP(0.05),
    sun_data     = sun_tb,
    shadow_model = Conical(),
    R_Occulting  = R_mars,
)

# ── Combined Dynamics ────────────────────────────────────────────────────────
dynamics = CentralBodyDynamicsModel(
    mars_gravity,
    (sun_tb, srp),
)

# ── Compute acceleration at circular Mars orbit (300 km altitude) ────────────
v_circ = sqrt(μ_mars / (R_mars + 300.0))
u0 = [R_mars + 300.0, 0.0, 0.0, 0.0, v_circ, 0.0]  # km, km/s
total_accel = build_dynamics_model(u0, p, 0.0, dynamics)
println("Total acceleration: ", total_accel, " km/s²")
```
