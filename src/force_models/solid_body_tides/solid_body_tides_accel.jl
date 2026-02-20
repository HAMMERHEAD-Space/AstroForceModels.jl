# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Solid Body Tides Acceleration Model
#
#   Computes the gravitational acceleration perturbation on a satellite due to the
#   redistribution of a central body's mass from tidal deformation, following
#   IERS Conventions (2010), Section 6.2 (Step 1).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note 36,
#       Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 2010.
#       Section 6.2: Effect of Solid Earth Tides (Equations 6.6, 6.7).
#
#   [2] Montenbruck, O. and Gill, E., Satellite Orbits: Models, Methods, and Applications,
#       Springer, 2000. Section 3.2.5: Solid Earth Tides.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TODO: SPIN THIS OUT INTO ITS OWN PACKAGE WITH FULL OCEAN TIDES, POLE TIDES, ECT.
export SolidBodyTidesModel, solid_body_tides_accel

"""
    SolidBodyTidesModel{N,TBT,K2T,K3T,KP2T,RET,BT} <: AbstractNonPotentialBasedForce

Solid body tides acceleration model following IERS Conventions (2010), Section 6.2.

Models the perturbation to the geopotential caused by the tidal deformation of a central
body due to the gravitational attraction of tide-raising bodies. Implements Step 1 of the
IERS procedure using frequency-independent Love numbers and the Legendre polynomial
addition theorem for efficient computation directly in the inertial frame.

The formulation is general and works for any central body (Earth, Mars, etc.) with
any number of tide-raising bodies. For Earth, the dominant tide raisers are the Moon
and Sun; for Mars they would be the Sun, Phobos, and Deimos.

The tidal perturbation potential at satellite position ``\\mathbf{r}`` from a tide-raising
body ``j`` at degree ``n`` is:

```math
\\Delta V_j^{(n)} = k_n \\frac{GM_j R_e^{2n+1}}{r_j^{n+1} r^{n+1}} P_n(\\cos\\gamma_j)
```

where ``k_n`` is the Love number, ``R_e`` is the central body's equatorial radius, ``r_j``
is the distance to body ``j``, and ``\\gamma_j`` is the geocentric angle between the
satellite and body ``j``.

# Fields
- `tide_raising_bodies::NTuple{N,ThirdBodyModel}`: Tuple of `ThirdBodyModel`s for each
  tide-raising body. Gravitational parameters are obtained from each body's `CelestialBody`.
- `k2::Number`: Degree-2 Love number (default: 0.30190, IERS 2010 Table 6.3, anelastic).
- `k3::Number`: Degree-3 Love number (default: 0.093, IERS 2010 Table 6.3).
- `k2_plus::Number`: Degree-2 → degree-4 coupling Love number ``k_2^{(+)}``
  (default: -0.00089, IERS 2010 Table 6.3, anelastic m=0). Set to 0 to disable.
- `R_e::Number`: Central body equatorial radius [km] (default: `R_EARTH`).
- `include_degree_3::Bool`: Include degree-3 tidal contribution (default: true).
"""
Base.@kwdef struct SolidBodyTidesModel{
    N,
    TBT<:NTuple{N,ThirdBodyModel},
    K2T<:Number,
    K3T<:Number,
    KP2T<:Number,
    RET<:Number,
    BT<:Bool,
} <: AbstractNonPotentialBasedForce
    tide_raising_bodies::TBT
    k2::K2T = 0.30190
    k3::K3T = 0.093
    k2_plus::KP2T = 0.0
    R_e::RET = R_EARTH
    include_degree_3::BT = true
end

"""
    SolidBodyTidesModel(eop_data; kwargs...)

Convenience constructor for Earth that creates Sun and Moon tide-raising body models
sharing a single EOP dataset.

# Arguments
- `eop_data::Union{EopIau1980,EopIau2000A}`: Earth Orientation Parameter data.

# Keyword Arguments
All fields of [`SolidBodyTidesModel`](@ref) may be overridden.
"""
function SolidBodyTidesModel(eop_data::Union{EopIau1980,EopIau2000A}; kwargs...)
    defaults = (
        tide_raising_bodies=(
            ThirdBodyModel(; body=SunBody(), eop_data=eop_data),
            ThirdBodyModel(; body=MoonBody(), eop_data=eop_data),
        ),
    )
    return SolidBodyTidesModel(; merge(defaults, kwargs)...)
end

"""
    acceleration(u, p, t, model::SolidBodyTidesModel) -> SVector{3}

Compute the solid body tides acceleration perturbation on a satellite.

# Arguments
- `u::AbstractVector`: Spacecraft state [rx, ry, rz, vx, vy, vz] in km and km/s (J2000 ECI).
- `p::ComponentVector`: Simulation parameters (must contain `JD`).
- `t::Number`: Elapsed time since epoch [s].
- `model::SolidBodyTidesModel`: Solid body tides model configuration.

# Returns
- `SVector{3}`: Acceleration perturbation [km/s²] in J2000 ECI.
"""
@inline function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, model::SolidBodyTidesModel
)
    jd = current_jd(p, t)
    bodies = _collect_body_data(jd, model.tide_raising_bodies)
    return solid_body_tides_accel(
        u,
        bodies;
        k2=model.k2,
        k3=model.k3,
        R_e=model.R_e,
        include_degree_3=model.include_degree_3,
    )
end

@inline function _collect_body_data(jd, bodies::Tuple)
    body = first(bodies)
    r = body(jd, Position()) ./ 1E3
    μ = body.body.μ
    return ((r, μ), _collect_body_data(jd, Base.tail(bodies))...)
end

@inline _collect_body_data(jd, ::Tuple{}) = ()

"""
    solid_body_tides_accel(
        u, bodies;
        k2=0.30190, k3=0.093,
        R_e=R_EARTH, include_degree_3=true,
    ) -> SVector{3}

Compute the solid body tides acceleration on a satellite from tide-raising body data.

Uses the gradient of the tidal perturbation potential derived from IERS Conventions (2010)
Equation 6.6 with the Legendre polynomial addition theorem.

**Degree 2** acceleration from body ``j``:
```math
\\mathbf{a}_j^{(2)} = \\frac{3 k_2 GM_j R_e^5}{r_j^3 r^4}
\\left[\\frac{1 - 5\\xi^2}{2}\\hat{\\mathbf{r}} + \\xi\\hat{\\mathbf{r}}_j\\right]
```

**Degree 3** acceleration from body ``j``:
```math
\\mathbf{a}_j^{(3)} = \\frac{k_3 GM_j R_e^7}{2 r_j^4 r^5}
\\left[(15\\xi - 35\\xi^3)\\hat{\\mathbf{r}} + (15\\xi^2 - 3)\\hat{\\mathbf{r}}_j\\right]
```

where ``\\xi = \\hat{\\mathbf{r}} \\cdot \\hat{\\mathbf{r}}_j``.

# Arguments
- `u::AbstractVector`: Spacecraft state vector [rx, ry, rz, vx, vy, vz] in km, km/s.
- `bodies::Tuple`: Tuple of `(r_body, μ_body)` pairs, where `r_body` is the body position
  [km] and `μ_body` is the gravitational parameter [km³/s²].

# Keyword Arguments
- `k2::Number`: Degree-2 Love number (default: 0.30190).
- `k3::Number`: Degree-3 Love number (default: 0.093).
- `R_e::Number`: Central body equatorial radius [km] (default: `R_EARTH`).
- `include_degree_3::Bool`: Include degree-3 contribution (default: true).

# Returns
- `SVector{3}`: Tidal acceleration [km/s²].
"""
@inline function solid_body_tides_accel(
    u::AbstractVector{UT},
    bodies::Tuple;
    k2::K2T=0.30190,
    k3::K3T=0.093,
    R_e::RET=R_EARTH,
    include_degree_3::Bool=true,
) where {UT,K2T,K3T,RET}
    r_sat = SVector{3,UT}(u[1], u[2], u[3])
    r_norm = norm(r_sat)
    r_hat = r_sat / r_norm

    accel = _sum_body_degree2(r_hat, r_norm, k2, R_e, bodies)

    if include_degree_3
        accel = accel + _sum_body_degree3(r_hat, r_norm, k3, R_e, bodies)
    end

    return accel
end

# -- Compile-time recursion over bodies for degree-2 contributions --

@inline function _sum_body_degree2(r_hat, r_norm, k2, R_e, bodies::Tuple)
    r_body, μ_body = first(bodies)
    current = _tidal_degree2(r_hat, r_norm, r_body, μ_body, k2, R_e)
    rest = _sum_body_degree2(r_hat, r_norm, k2, R_e, Base.tail(bodies))
    return SVector{3}(current[1] + rest[1], current[2] + rest[2], current[3] + rest[3])
end

@inline function _sum_body_degree2(r_hat, r_norm, k2, R_e, bodies::Tuple{<:Any})
    r_body, μ_body = first(bodies)
    return _tidal_degree2(r_hat, r_norm, r_body, μ_body, k2, R_e)
end

# -- Compile-time recursion over bodies for degree-3 contributions --

@inline function _sum_body_degree3(r_hat, r_norm, k3, R_e, bodies::Tuple)
    r_body, μ_body = first(bodies)
    current = _tidal_degree3(r_hat, r_norm, r_body, μ_body, k3, R_e)
    rest = _sum_body_degree3(r_hat, r_norm, k3, R_e, Base.tail(bodies))
    return SVector{3}(current[1] + rest[1], current[2] + rest[2], current[3] + rest[3])
end

@inline function _sum_body_degree3(r_hat, r_norm, k3, R_e, bodies::Tuple{<:Any})
    r_body, μ_body = first(bodies)
    return _tidal_degree3(r_hat, r_norm, r_body, μ_body, k3, R_e)
end

# -- Single-body tidal acceleration computations --

"""
Degree-2 tidal acceleration from a single perturbing body.

Gradient of the degree-2 tidal perturbation potential:
  ``\\Delta V^{(2)} = k_2 GM_j R_e^5 / (r_j^3 r^3) \\cdot P_2(\\cos\\gamma)``
"""
@inline function _tidal_degree2(
    r_hat::SVector{3},
    r_norm::Number,
    r_body::AbstractVector,
    μ_body::Number,
    k2::Number,
    Re::Number,
)
    r_body_norm = norm(r_body)
    r_body_hat = r_body / r_body_norm

    ξ = dot(r_hat, r_body_hat)

    coeff = 3.0 * k2 * μ_body * Re^5 / (r_body_norm^3 * r_norm^4)
    radial = (1.0 - 5.0 * ξ^2) / 2.0
    body_dir = ξ

    return coeff * (radial * r_hat + body_dir * r_body_hat)
end

"""
Degree-3 tidal acceleration from a single perturbing body.

Gradient of the degree-3 tidal perturbation potential:
  ``\\Delta V^{(3)} = k_3 GM_j R_e^7 / (r_j^4 r^4) \\cdot P_3(\\cos\\gamma)``
"""
@inline function _tidal_degree3(
    r_hat::SVector{3},
    r_norm::Number,
    r_body::AbstractVector,
    μ_body::Number,
    k3::Number,
    Re::Number,
)
    r_body_norm = norm(r_body)
    r_body_hat = r_body / r_body_norm

    ξ = dot(r_hat, r_body_hat)

    coeff = k3 * μ_body * Re^7 / (2.0 * r_body_norm^4 * r_norm^5)
    radial = 15.0 * ξ - 35.0 * ξ^3
    body_dir = 15.0 * ξ^2 - 3.0

    return coeff * (radial * r_hat + body_dir * r_body_hat)
end
