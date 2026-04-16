# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from the geomagnetic Lorentz force on a charged spacecraft
#
#   Earth-only. Requires IGRF model and ITRF frame in FrameSystem.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Peck, M. A. (2005). "Prospects and Challenges for Lorentz-Augmented Orbits."
#       AIAA Guidance, Navigation, and Control Conference. AIAA 2005-5995.
#
#   [2] Streetman, B. & Peck, M. A. (2007). "New Synchronous Orbits Using the
#       Geomagnetic Lorentz Force." Journal of Guidance, Control, and Dynamics, 30(6),
#       1677-1690. https://doi.org/10.2514/1.29080
#
#   [3] Khalil, K. I. & Abdel-Aziz, Y. A. (2014). "Electromagnetic effects on the
#       orbital motion of a charged spacecraft." Research in Astronomy and Astrophysics,
#       14(5), 589. https://doi.org/10.1088/1674-4527/14/5/008
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export MagneticFieldAstroModel, magnetic_field_accel
export GeomagneticFieldType, IGRFField, DipoleMagneticField

"""
    GeomagneticFieldType

Abstract type for geomagnetic field model selection.
"""
abstract type GeomagneticFieldType end

"""
    IGRFField <: GeomagneticFieldType

International Geomagnetic Reference Field (IGRF) v14.
"""
struct IGRFField <: GeomagneticFieldType end

"""
    DipoleMagneticField <: GeomagneticFieldType

Simplified geomagnetic dipole model.
"""
struct DipoleMagneticField <: GeomagneticFieldType end

"""
    MagneticFieldAstroModel

Geomagnetic Lorentz force model for charged spacecraft.

Earth-only. Uses the FrameSystem rotation for ECI↔ECEF coordinate transformations
required by the IGRF and dipole magnetic field models.

# Constructor
    MagneticFieldAstroModel(; spacecraft_charge_model, geomagnetic_field_model=DipoleMagneticField(),
        body_fixed_frame=:ITRF, propagation_frame=:ICRF, frames=nothing,
        max_degree=13, P=nothing, dP=nothing)

# Arguments
- `spacecraft_charge_model`: Model providing charge-to-mass ratio q/m [C/kg].
- `geomagnetic_field_model`: `IGRFField()` or `DipoleMagneticField()`.
- `body_fixed_frame::Symbol`: Body-fixed frame for field computation (default: `:ITRF`).
- `propagation_frame::Symbol`: Frame in which the state vector is expressed (default: `:ICRF`).
- `frames`: Optional `FrameSystem`. When provided, pre-compiles the rotation between
  `propagation_frame` and `body_fixed_frame` for allocation-free evaluation.
- `max_degree::Int`: Maximum spherical harmonic degree for IGRF (default: 13).
- `P`: Pre-allocated Legendre polynomial matrix.
- `dP`: Pre-allocated Legendre derivative matrix.
"""
struct MagneticFieldAstroModel{
    CT<:AbstractSpacecraftChargeModel,
    GT<:GeomagneticFieldType,
    PT<:Union{Nothing,AbstractMatrix},
    DPT<:Union{Nothing,AbstractMatrix},
    CR3,
} <: AbstractNonPotentialBasedForce
    spacecraft_charge_model::CT
    geomagnetic_field_model::GT
    body_fixed_frame::Symbol
    propagation_frame::Symbol
    max_degree::Int
    P::PT
    dP::DPT
    compiled_rotation3::CR3
end

function MagneticFieldAstroModel(;
    spacecraft_charge_model::CT,
    geomagnetic_field_model::GT=DipoleMagneticField(),
    body_fixed_frame::Symbol=:ITRF,
    propagation_frame::Symbol=:ICRF,
    frames=nothing,
    max_degree::Int=13,
    P::PT=nothing,
    dP::DPT=nothing,
    # Legacy: accept eop_data but ignore it
    eop_data=nothing,
) where {CT<:AbstractSpacecraftChargeModel,GT<:GeomagneticFieldType,PT,DPT}
    cr3 = nothing
    if !isnothing(frames) && propagation_frame != body_fixed_frame
        cr3 = compile_rotation(frames, propagation_frame, body_fixed_frame, Val(1))
    end
    return MagneticFieldAstroModel(
        spacecraft_charge_model,
        geomagnetic_field_model,
        body_fixed_frame,
        propagation_frame,
        max_degree,
        P,
        dP,
        cr3,
    )
end

"""
    acceleration(u::AbstractVector, p::FrameAwareParams, t::Number, model::MagneticFieldAstroModel)

Compute the geomagnetic Lorentz force acceleration on a charged spacecraft.

Earth-only.

# Arguments
- `u::AbstractVector`: Current state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `model::MagneticFieldAstroModel`: Magnetic field model.

# Returns
- `SVector{3}`: Lorentz force acceleration [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, model::MagneticFieldAstroModel
)
    JD = current_jd(p, t)
    t_ft = ft_time(p, t)

    q_over_m = charge_mass_ratio(u, p, t, model.spacecraft_charge_model)

    R_rot = if isnothing(model.compiled_rotation3)
        rotation3(p.frames, model.propagation_frame, model.body_fixed_frame, t_ft)
    else
        model.compiled_rotation3(t_ft)
    end
    R_prop2bf = R_rot.m[1]

    r_eci = SVector{3}(u[1], u[2], u[3])
    r_ecef = R_prop2bf * r_eci

    B_ecef = compute_magnetic_field(r_ecef, JD, model)

    R_bf2prop = R_prop2bf'
    B_eci = R_bf2prop * B_ecef

    return magnetic_field_accel(u, q_over_m, B_eci)
end

@inline function compute_magnetic_field(
    r_ecef::SVector{3}, JD::Number, model::MagneticFieldAstroModel{CT,IGRFField}
) where {CT}
    r_m = r_ecef * 1E3

    x, y, z = r_m[1], r_m[2], r_m[3]
    r_norm = norm(r_m)
    ρ = sqrt(x^2 + y^2)

    λ = atan(z, ρ)
    Ω = atan(y, x)

    decimal_year = 2000.0 + (JD - 2451545.0) / 365.25

    B_ned = igrf(
        decimal_year,
        r_norm,
        λ,
        Ω;
        max_degree=model.max_degree,
        P=model.P,
        dP=model.dP,
        show_warnings=false,
    )

    sinλ, cosλ = sincos(λ)
    sinΩ, cosΩ = sincos(Ω)

    B_N, B_E, B_D = B_ned[1], B_ned[2], B_ned[3]

    B_ecef_x = -sinλ * cosΩ * B_N - sinΩ * B_E - cosλ * cosΩ * B_D
    B_ecef_y = -sinλ * sinΩ * B_N + cosΩ * B_E - cosλ * sinΩ * B_D
    B_ecef_z = cosλ * B_N - sinλ * B_D

    return SVector{3}(B_ecef_x, B_ecef_y, B_ecef_z) * 1E-9
end

@inline function compute_magnetic_field(
    r_ecef::SVector{3}, JD::Number, model::MagneticFieldAstroModel{CT,DipoleMagneticField}
) where {CT}
    r_m = r_ecef * 1E3

    decimal_year = 2000.0 + (JD - 2451545.0) / 365.25

    B_ecef_nT = geomagnetic_dipole_field(r_m, decimal_year)

    return SVector{3}(B_ecef_nT[1], B_ecef_nT[2], B_ecef_nT[3]) * 1E-9
end

"""
    magnetic_field_accel(u, q_over_m, B_eci) -> SVector{3}

Compute the geomagnetic Lorentz force acceleration.

# Returns
- `SVector{3}`: Lorentz force acceleration [km/s²].
"""
@inline function magnetic_field_accel(
    u::AbstractVector{UT}, q_over_m::QT, B_eci::SVector{3,BT}
) where {UT,QT,BT}
    RT = promote_type(UT, QT, BT)

    r = SVector{3,UT}(u[1], u[2], u[3])
    v = SVector{3,UT}(u[4], u[5], u[6])

    ω_vec = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)
    v_rel = v - cross(ω_vec, r)

    v_rel_m = v_rel * 1E3
    a_m_s2 = q_over_m * cross(v_rel_m, B_eci)

    return SVector{3,RT}(a_m_s2[1] / 1E3, a_m_s2[2] / 1E3, a_m_s2[3] / 1E3)
end
