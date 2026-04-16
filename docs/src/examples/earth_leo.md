# Earth LEO Satellite

High-fidelity force model for a low-Earth orbit satellite with spherical harmonics
gravity, atmospheric drag, SRP, third-body Sun and Moon, solid-body tides, relativistic
effects, albedo, and thermal emission.

This example is fully self-contained and can be run as-is (required data files are
fetched automatically on first use).

```julia
using AstroForceModels
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxGravityModels
using ComponentArrays
using FrameTransformations
using Tempo
using StaticArrays
using SpaceIndices

# Initialize space weather indices (required for JB2008 drag model)
SpaceIndices.init()

# ── Epoch & EOP ──────────────────────────────────────────────────────────────
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
eop_data = fetch_iers_eop()

# ── Frame System ─────────────────────────────────────────────────────────────
# setup_earth_propagation_frames builds ICRF + ITRF + Earth (399) + Sun (10) + Moon (301).
epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
p = setup_earth_propagation_frames(epoch, eop_data)
frames = p.frames

# ── Force Models ─────────────────────────────────────────────────────────────

# 1. Spherical harmonics gravity (36×36 EGM96, auto-downloaded)
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
gravity = GravityHarmonicsAstroModel(;
    gravity_model  = grav_coeffs,
    body_fixed_frame = :ITRF,
    propagation_frame = :ICRF,
    order = 36, degree = 36,
    P  = MMatrix{37,37,Float64}(zeros(37, 37)),
    dP = MMatrix{37,37,Float64}(zeros(37, 37)),
    frames = frames,
)

# 2. Atmospheric drag (JB2008)
drag = DragAstroModel(;
    satellite_drag_model = CannonballFixedDrag(0.2),
    atmosphere_model     = JB2008(),
    eop_data             = eop_data,
)

# 3. Third-body Sun & Moon
sun_tb = ThirdBodyModel(
    body = SunBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
    frames = frames,
)
moon_tb = ThirdBodyModel(
    body = MoonBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
    frames = frames,
)

# 4. Solar radiation pressure (conical shadow)
srp = SRPAstroModel(;
    satellite_srp_model = CannonballFixedSRP(0.03),
    sun_data     = sun_tb,
    shadow_model = Conical(),
    R_Occulting  = AstroForceModels.R_EARTH,
)

# 5. Relativistic effects (Schwarzschild + Lense-Thirring + de Sitter)
earth_body = ThirdBodyModel(
    body = EarthBody(),
    ephem_type = FrameEphemeris(center_point=399, target_point=399, axes=:ICRF),
)
relativity = RelativityModel(;
    central_body = earth_body,
    sun_body     = sun_tb,
    J = SVector{3}(0.0, 0.0, AstroForceModels.EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS),
)

# 6. Solid-body tides (degree 2+3, Sun+Moon)
tides = SolidBodyTidesModel(;
    tide_raising_bodies = (sun_tb, moon_tb),
    R_e = AstroForceModels.R_EARTH,
)

# 7. Earth albedo radiation pressure (Lebedev quadrature)
albedo = AlbedoAstroModel(;
    satellite_shape_model = CannonballFixedSRP(0.03),
    sun_data          = sun_tb,
    body_albedo_model = UniformAlbedoModel(0.3, 0.7),
    body_fixed_frame  = :ITRF,
    propagation_frame = :ICRF,
    frames            = frames,
)

# 8. Thermal emission
thermal = ThermalEmissionAstroModel(;
    satellite_thermal_model = FixedThermalEmission(0.01),
    sun_data     = sun_tb,
    shadow_model = Conical(),
    R_Occulting  = AstroForceModels.R_EARTH,
)

# ── Combined Dynamics ────────────────────────────────────────────────────────
dynamics = CentralBodyDynamicsModel(
    gravity,
    (drag, srp, sun_tb, moon_tb, relativity, tides, albedo, thermal),
)

# ── Compute acceleration at ISS-like state ───────────────────────────────────
u0 = [6378.137 + 408.0, 0.0, 0.0, 0.0, 7.660, 0.0]  # km, km/s
total_accel = build_dynamics_model(u0, p, 0.0, dynamics)
println("Total acceleration: ", total_accel, " km/s²")
```
