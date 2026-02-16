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
] #km, km/s

const _eop_data = fetch_iers_eop()
const _grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

const _grav_model = GravityHarmonicsAstroModel(;
    gravity_model=_grav_coeffs, eop_data=_eop_data, order=36, degree=36
)

const _ATMOSPHERE_MODELS = (
    ("JB2008", JB2008()),
    ("JR1971", JR1971()),
    ("MSIS2000", MSIS2000()),
    ("ExpAtmo", ExpAtmo()),
    ("NoAtmosphere", NoAtmosphere()),
)

const _satellite_drag_model = CannonballFixedDrag(0.2)
const _drag_model = DragAstroModel(;
    satellite_drag_model=_satellite_drag_model,
    atmosphere_model=JB2008(),
    eop_data=_eop_data,
)
const _BC = 0.2

const _ENZYME_RUNTIME_ACTIVITY = ["MSIS2000"]

const _relativity_model = RelativityModel()

const _satellite_srp_model = CannonballFixedSRP(0.2)
const _sun_model = ThirdBodyModel(; body=SunBody(), eop_data=_eop_data)

const _SHADOW_MODELS = (
    ("Conical", Conical()),
    ("Cylindrical", Cylindrical()),
    ("SmoothedConical", SmoothedConical()),
    ("No_Shadow", No_Shadow()),
)

const _srp_model = SRPAstroModel(;
    satellite_srp_model=_satellite_srp_model,
    sun_data=_sun_model,
    eop_data=_eop_data,
    shadow_model=Conical(),
)
const _RC = 0.2

const _moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=_eop_data)

const _model_list = CentralBodyDynamicsModel(
    _grav_model, (_sun_model, _moon_model, _srp_model, _drag_model)
)
