# Heliocentric Asteroid with Planetary Perturbations

Propagation of a fictitious near-Earth asteroid orbiting the Sun, perturbed by
solar radiation pressure and the gravitational pull of all eight planets.

This example requires the [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl)
package and a JPL DE440 SPK kernel (`de440.bsp`). The SPK file can be downloaded from
[NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp).

```julia
using AstroForceModels
using ComponentArrays
using FrameTransformations
using Tempo
using StaticArrays
using Ephemerides

# ── Epoch ────────────────────────────────────────────────────────────────────
JD = date_to_jd(2030, 7, 1, 0, 0, 0.0)

# ── Planetary GM values [km³/s²] (DE440 values) ─────────────────────────────
const GM_MERCURY = 2.2032080486417923e4
const GM_VENUS   = 3.2485859200000006e5
const GM_EARTH   = 3.986004415e5
const GM_MARS    = 4.282837566395650e4
const GM_JUPITER = 1.266865349115908e8
const GM_SATURN  = 3.793120623436167e7
const GM_URANUS  = 5.793951256527211e6
const GM_NEPTUNE = 6.836527100580532e6

# ── Frame System with SPK ephemeris ──────────────────────────────────────────
frames = FrameSystem{2, Float64}()
add_axes_icrf!(frames)

# Load JPL DE440 ephemeris
eph = EphemerisProvider("de440.bsp")

# Solar System Barycenter as root point
add_point!(frames, :SSB, 0, :ICRF)

# Sun and all 8 planet-system barycenters as children of SSB
add_point_ephemeris!(frames, eph, :Sun,      10, 0)
add_point_ephemeris!(frames, eph, :MercuryB,  1, 0)
add_point_ephemeris!(frames, eph, :VenusB,    2, 0)
add_point_ephemeris!(frames, eph, :EMB,       3, 0)   # Earth-Moon Barycenter
add_point_ephemeris!(frames, eph, :MarsB,     4, 0)
add_point_ephemeris!(frames, eph, :JupiterB,  5, 0)
add_point_ephemeris!(frames, eph, :SaturnB,   6, 0)
add_point_ephemeris!(frames, eph, :UranusB,   7, 0)
add_point_ephemeris!(frames, eph, :NeptuneB,  8, 0)

# ── Parameters (heliocentric propagation) ────────────────────────────────────
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = FrameAwareParams(frames, epoch, :ICRF)

# ── Central gravity: Sun (Keplerian point mass) ─────────────────────────────
sun_gravity = KeplerianGravityAstroModel(μ=AstroForceModels.μ_SUN)

# ── Helper: build a ThirdBodyModel for each planet ───────────────────────────
# center_point = 10 (Sun) because we propagate heliocentric.
# target_point = planet barycenter NAIF ID.
function planet_model(name::Symbol, naif_id::Int, μ::Float64, R::Float64)
    body = CelestialBody(name, :Sun, naif_id, μ, R)
    return ThirdBodyModel(
        body       = body,
        ephem_type = FrameEphemeris(center_point=10, target_point=naif_id, axes=:ICRF),
        frames     = frames,
    )
end

mercury = planet_model(:Mercury, 1, GM_MERCURY, 2439.7)
venus   = planet_model(:Venus,   2, GM_VENUS,   6051.8)
emb     = planet_model(:EMB,     3, GM_EARTH,   6371.0)
mars    = planet_model(:Mars,    4, GM_MARS,    3389.5)
jupiter = planet_model(:Jupiter, 5, GM_JUPITER, 69911.0)
saturn  = planet_model(:Saturn,  6, GM_SATURN,  58232.0)
uranus  = planet_model(:Uranus,  7, GM_URANUS,  25362.0)
neptune = planet_model(:Neptune, 8, GM_NEPTUNE, 24622.0)

# ── Solar radiation pressure ─────────────────────────────────────────────────
# For a heliocentric orbit the spacecraft position IS the Sun→spacecraft vector,
# so we set center=target=Sun to place the Sun at the origin.
sun_at_origin = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=10, target_point=10, axes=:ICRF),
)

srp = SRPAstroModel(;
    satellite_srp_model = CannonballFixedSRP(0.01),   # A/m · Cr  [m²/kg]
    sun_data     = sun_at_origin,
    shadow_model = NoShadow(),   # no occulting body in heliocentric space
    R_Occulting  = 0.0,
)

# ── Combined Dynamics ────────────────────────────────────────────────────────
dynamics = CentralBodyDynamicsModel(
    sun_gravity,
    (mercury, venus, emb, mars, jupiter, saturn, uranus, neptune, srp),
)

# ── Initial state: fictitious NEA at ~1.2 AU perihelion, a ≈ 1.5 AU, i ≈ 10° ──
AU_km = AstroForceModels.ASTRONOMICAL_UNIT / 1e3
r0 = 1.2 * AU_km
v0 = sqrt(AstroForceModels.μ_SUN * (2 / r0 - 1 / (1.5 * AU_km)))   # vis-viva

incl = deg2rad(10.0)
u0 = [r0, 0.0, 0.0, 0.0, v0 * cos(incl), v0 * sin(incl)]   # km, km/s

total_accel = build_dynamics_model(u0, p, 0.0, dynamics)
println("Total acceleration: ", total_accel, " km/s²")
```
