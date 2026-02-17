using AstroForceModels
using BenchmarkTools
using ComponentArrays
using LinearAlgebra
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using StaticArraysCore

const SUITE = BenchmarkGroup()

SUITE["gravity"] = BenchmarkGroup(["acceleration"])
SUITE["drag"] = BenchmarkGroup(["acceleration"])
SUITE["srp"] = BenchmarkGroup(["acceleration"])
SUITE["albedo"] = BenchmarkGroup(["acceleration"])
SUITE["third_body"] = BenchmarkGroup(["acceleration"])
SUITE["relativity"] = BenchmarkGroup(["acceleration"])
SUITE["low_thrust"] = BenchmarkGroup(["acceleration"])
SUITE["shadow_models"] = BenchmarkGroup(["shadow"])
SUITE["dynamics_builder"] = BenchmarkGroup(["combined"])

# ---------------------
# Common test state
# ---------------------
const _JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
const _p = ComponentVector(; JD=_JD)
const _t = 0.0

const _state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

const _eop_data = fetch_iers_eop()

# ---------------------
# Gravity models
# ---------------------
const _grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

const _keplerian_model = KeplerianGravityAstroModel()

const _harmonics_models = [
    (
        "4x4",
        GravityHarmonicsAstroModel(;
            gravity_model=_grav_coeffs, eop_data=_eop_data, order=4, degree=4
        ),
    ),
    (
        "20x20",
        GravityHarmonicsAstroModel(;
            gravity_model=_grav_coeffs, eop_data=_eop_data, order=20, degree=20
        ),
    ),
    (
        "36x36",
        GravityHarmonicsAstroModel(;
            gravity_model=_grav_coeffs, eop_data=_eop_data, order=36, degree=36
        ),
    ),
]

SUITE["gravity"]["Keplerian"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_keplerian_model
)

for (label, model) in _harmonics_models
    SUITE["gravity"]["Harmonics $label"] = @benchmarkable acceleration(
        $_state, $_p, $_t, $model
    )
end

# ---------------------
# Drag models
# ---------------------
const _satellite_drag_model = CannonballFixedDrag(0.2)

const _ATMOSPHERE_MODELS = [
    ("JB2008", JB2008()),
    ("JR1971", JR1971()),
    ("MSIS2000", MSIS2000()),
    ("ExpAtmo", ExpAtmo()),
]

for (label, atmo) in _ATMOSPHERE_MODELS
    drag_model = DragAstroModel(;
        satellite_drag_model=_satellite_drag_model,
        atmosphere_model=atmo,
        eop_data=_eop_data,
    )
    SUITE["drag"][label] = @benchmarkable acceleration($_state, $_p, $_t, $drag_model)
end

# ---------------------
# SRP models
# ---------------------
const _satellite_srp_model = CannonballFixedSRP(0.2)
const _sun_model = ThirdBodyModel(; body=SunBody(), eop_data=_eop_data)

const _SHADOW_MODELS = [
    ("Conical", Conical()),
    ("Cylindrical", Cylindrical()),
    ("SmoothedConical", SmoothedConical()),
    ("NoShadow", NoShadow()),
]

for (label, shadow) in _SHADOW_MODELS
    srp_model = SRPAstroModel(;
        satellite_srp_model=_satellite_srp_model,
        sun_data=_sun_model,
        eop_data=_eop_data,
        shadow_model=shadow,
    )
    SUITE["srp"][label] = @benchmarkable acceleration($_state, $_p, $_t, $srp_model)
end

# Shadow model micro-benchmarks
const _sun_pos = _sun_model(_JD, Position()) ./ 1E3
const _sat_pos = SVector{3}(_state[1], _state[2], _state[3])

for (label, shadow) in _SHADOW_MODELS
    SUITE["shadow_models"][label] = @benchmarkable shadow_model(
        $_sat_pos, $_sun_pos, $shadow
    )
end

# ---------------------
# Albedo models
# ---------------------
const _uniform_albedo = UniformAlbedoModel(; visible_albedo=0.3, infrared_emissivity=0.7)

const _ALBEDO_ORDERS = [
    ("order_21", 21),
    ("order_59", 59),
    ("order_125", 125),
]

for (label, order) in _ALBEDO_ORDERS
    albedo_model = AlbedoAstroModel(;
        satellite_shape_model=_satellite_srp_model,
        sun_data=_sun_model,
        body_albedo_model=_uniform_albedo,
        eop_data=_eop_data,
        lebedev_order=order,
    )
    SUITE["albedo"][label] = @benchmarkable acceleration(
        $_state, $_p, $_t, $albedo_model
    )
end

# ---------------------
# Third body models
# ---------------------
const _moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=_eop_data)

SUITE["third_body"]["Sun"] = @benchmarkable acceleration($_state, $_p, $_t, $_sun_model)
SUITE["third_body"]["Moon"] = @benchmarkable acceleration($_state, $_p, $_t, $_moon_model)

# ---------------------
# Relativity model
# ---------------------
const _relativity_model = RelativityModel(_eop_data)

SUITE["relativity"]["Full"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_relativity_model
)
SUITE["relativity"]["Schwarzschild only"] = @benchmarkable acceleration(
    $_state,
    $_p,
    $_t,
    $(RelativityModel(_eop_data; lense_thirring_effect=false, de_Sitter_effect=false)),
)

# ---------------------
# Low thrust models
# ---------------------
const _lt_cartesian = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(1e-7, 2e-7, 3e-7)
)
const _lt_tangential = LowThrustAstroModel(; thrust_model=ConstantTangentialThrust(1e-7))
const _lt_rtn = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(0.0, 1e-7, 0.0), frame=RTNFrame()
)
const _lt_vnb = LowThrustAstroModel(;
    thrust_model=ConstantCartesianThrust(1e-7, 0.0, 0.0), frame=VNBFrame()
)
const _lt_piecewise = LowThrustAstroModel(;
    thrust_model=PiecewiseConstantThrust(
        [0.0, 3600.0, 7200.0],
        [
            SVector{3}(1e-7, 0.0, 0.0),
            SVector{3}(0.0, 1e-7, 0.0),
            SVector{3}(-1e-7, 0.0, 0.0),
        ],
    ),
    frame=RTNFrame(),
)

SUITE["low_thrust"]["Cartesian (Inertial)"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_lt_cartesian
)
SUITE["low_thrust"]["Tangential (Inertial)"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_lt_tangential
)
SUITE["low_thrust"]["Cartesian (RTN)"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_lt_rtn
)
SUITE["low_thrust"]["Cartesian (VNB)"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_lt_vnb
)
SUITE["low_thrust"]["Piecewise (RTN, 3 arcs)"] = @benchmarkable acceleration(
    $_state, $_p, $_t, $_lt_piecewise
)

# ---------------------
# Dynamics builder
# ---------------------
const _grav_model_36 = GravityHarmonicsAstroModel(;
    gravity_model=_grav_coeffs, eop_data=_eop_data, order=36, degree=36
)
const _drag_model = DragAstroModel(;
    satellite_drag_model=_satellite_drag_model,
    atmosphere_model=JB2008(),
    eop_data=_eop_data,
)
const _srp_model = SRPAstroModel(;
    satellite_srp_model=_satellite_srp_model,
    sun_data=_sun_model,
    eop_data=_eop_data,
    shadow_model=Conical(),
)

const _dynamics_keplerian = CentralBodyDynamicsModel(_keplerian_model)
const _dynamics_harmonics = CentralBodyDynamicsModel(_grav_model_36)
const _dynamics_full = CentralBodyDynamicsModel(
    _grav_model_36,
    (_sun_model, _moon_model, _srp_model, _drag_model, _relativity_model, _lt_tangential),
)

SUITE["dynamics_builder"]["Keplerian only"] = @benchmarkable build_dynamics_model(
    $_state, $_p, $_t, $_dynamics_keplerian
)
SUITE["dynamics_builder"]["Harmonics 36x36 only"] = @benchmarkable build_dynamics_model(
    $_state, $_p, $_t, $_dynamics_harmonics
)
SUITE["dynamics_builder"]["Full (gravity+drag+srp+3body+rel+lt)"] = @benchmarkable build_dynamics_model(
    $_state, $_p, $_t, $_dynamics_full
)

# ---------------------
# Tune and cache
# ---------------------
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, BenchmarkTools.params(SUITE))
end
