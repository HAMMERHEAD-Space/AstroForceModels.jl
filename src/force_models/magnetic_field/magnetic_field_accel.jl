# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from the geomagnetic Lorentz force on a charged spacecraft
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

Uses spherical harmonic expansion up to degree 13 for high-fidelity geomagnetic
field computation. Valid for dates between 1900 and 2035.

See SatelliteToolboxGeomagneticField.jl for details.
"""
struct IGRFField <: GeomagneticFieldType end

"""
    DipoleMagneticField <: GeomagneticFieldType

Simplified geomagnetic dipole model.

Assumes Earth's magnetic field is a perfect dipole. Less accurate than IGRF but
sufficient for preliminary analysis where uncertainties are high.
"""
struct DipoleMagneticField <: GeomagneticFieldType end

"""
    MagneticFieldAstroModel{CT,GT,EoT,PT,DPT} <: AbstractNonPotentialBasedForce

Geomagnetic Lorentz force model for charged spacecraft.

A charged spacecraft moving through Earth's magnetic field experiences a perturbative
Lorentz force. The magnetic field co-rotates with Earth, so the force depends on the
spacecraft's velocity relative to the rotating field [1,2]:

    𝐚 = (q/m) × (𝐯_rel × 𝐁)

where q/m is the charge-to-mass ratio, 𝐯_rel = 𝐯 - 𝛚×𝐫 is the velocity relative
to the co-rotating magnetic field, and 𝐁 is the geomagnetic field vector.

The force is always perpendicular to both the relative velocity and the magnetic
field, making it a non-dissipative perturbation that can modify orbital elements
without adding or removing energy from the orbit (in the rotating frame) [3].

# Type Parameters
- `CT <: AbstractSpacecraftChargeModel`: Charge-to-mass ratio model
- `GT <: GeomagneticFieldType`: Geomagnetic field model type
- `EoT <: Union{EopIau1980,EopIau2000A}`: Earth Orientation Parameters type
- `PT <: Union{Nothing,AbstractMatrix}`: Legendre polynomial buffer type
- `DPT <: Union{Nothing,AbstractMatrix}`: Legendre derivative buffer type

# Fields
- `spacecraft_charge_model::CT`: Model providing the charge-to-mass ratio q/m [C/kg]
- `geomagnetic_field_model::GT`: Geomagnetic field model — `IGRFField()` or `DipoleMagneticField()`
- `eop_data::EoT`: Earth Orientation Parameters for J2000 ↔ ITRF transformations
- `max_degree::Int`: Maximum spherical harmonic degree for IGRF (ignored for dipole). Default: 13
- `P::PT`: Pre-allocated Legendre polynomial matrix for IGRF (reduces allocations). Default: `nothing`
- `dP::DPT`: Pre-allocated Legendre derivative matrix for IGRF (reduces allocations). Default: `nothing`
"""
Base.@kwdef struct MagneticFieldAstroModel{
    CT<:AbstractSpacecraftChargeModel,
    GT<:GeomagneticFieldType,
    EoT<:Union{EopIau1980,EopIau2000A},
    PT<:Union{Nothing,AbstractMatrix},
    DPT<:Union{Nothing,AbstractMatrix},
} <: AbstractNonPotentialBasedForce
    spacecraft_charge_model::CT
    geomagnetic_field_model::GT = DipoleMagneticField()
    eop_data::EoT
    max_degree::Int = 13
    P::PT = nothing
    dP::DPT = nothing
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, model::MagneticFieldAstroModel)

Compute the geomagnetic Lorentz force acceleration on a charged spacecraft.

# Arguments
- `u::AbstractVector`: Current state [rx, ry, rz, vx, vy, vz] in km and km/s (J2000).
- `p::ComponentVector`: Parameters including `JD` (Julian Date at epoch).
- `t::Number`: Elapsed time since epoch [s].
- `model::MagneticFieldAstroModel`: Magnetic field force model configuration.

# Returns
- `SVector{3}`: Lorentz force acceleration [km/s²] in J2000 frame.
"""
@inline function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, model::MagneticFieldAstroModel
)
    JD = current_jd(p, t)

    q_over_m = charge_mass_ratio(u, p, t, model.spacecraft_charge_model)

    R_J2000_to_ITRF::SatelliteToolboxTransformations.DCM{Float64} = r_eci_to_ecef(
        DCM, J2000(), ITRF(), JD, model.eop_data
    )

    r_eci = SVector{3}(u[1], u[2], u[3])
    r_ecef = R_J2000_to_ITRF * r_eci

    B_ecef = compute_magnetic_field(r_ecef, JD, model)

    R_ITRF_to_J2000 = R_J2000_to_ITRF'
    B_eci = R_ITRF_to_J2000 * B_ecef

    return magnetic_field_accel(u, q_over_m, B_eci)
end

"""
    compute_magnetic_field(r_ecef::SVector{3}, JD::Number, model::MagneticFieldAstroModel{CT,IGRFField}) where CT

Compute the geomagnetic field vector in ECEF [T] using the IGRF model.

Converts ECEF Cartesian coordinates to geocentric spherical, calls the IGRF model
to get B in the geocentric NED frame, then rotates back to ECEF.
"""
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

"""
    compute_magnetic_field(r_ecef::SVector{3}, JD::Number, model::MagneticFieldAstroModel{CT,DipoleMagneticField}) where CT

Compute the geomagnetic field vector in ECEF [T] using the simplified dipole model.
"""
@inline function compute_magnetic_field(
    r_ecef::SVector{3}, JD::Number, model::MagneticFieldAstroModel{CT,DipoleMagneticField}
) where {CT}
    r_m = r_ecef * 1E3

    decimal_year = 2000.0 + (JD - 2451545.0) / 365.25

    B_ecef_nT = geomagnetic_dipole_field(r_m, decimal_year)

    return SVector{3}(B_ecef_nT[1], B_ecef_nT[2], B_ecef_nT[3]) * 1E-9
end

"""
    magnetic_field_accel(u::AbstractVector, q_over_m::Number, B_eci::SVector{3}) -> SVector{3}

Compute the geomagnetic Lorentz force acceleration given the state, charge-to-mass
ratio, and magnetic field vector in the inertial frame.

The Lorentz force on a charged body moving through a magnetic field is [1,2]:

    𝐅 = q(𝐯_rel × 𝐁)

where 𝐯_rel = 𝐯 - 𝛚×𝐫 is the velocity relative to the co-rotating magnetic field.
The acceleration is:

    𝐚 = (q/m)(𝐯_rel × 𝐁)

The magnetic field co-rotates with the Earth, so the relative velocity accounts for
Earth's rotation via the transport theorem, identical to the apparent velocity used
in atmospheric drag computations.

# Arguments
- `u::AbstractVector`: Spacecraft state [rx, ry, rz, vx, vy, vz] in km and km/s (J2000).
- `q_over_m::Number`: Charge-to-mass ratio [C/kg].
- `B_eci::SVector{3}`: Geomagnetic field vector in J2000 frame [T].

# Returns
- `SVector{3}`: Lorentz force acceleration [km/s²] in J2000 frame.
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
