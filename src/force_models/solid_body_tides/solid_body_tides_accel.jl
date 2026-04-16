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
export SolidBodyTidesModel, solid_body_tides_accel

"""
    SolidBodyTidesModel

Solid body tides acceleration model following IERS Conventions (2010), Section 6.2.

# Fields
- `tide_raising_bodies::NTuple{N,ThirdBodyModel}`: Tide-raising body models (use FrameEphemeris).
- `k2::Number`: Degree-2 Love number (default: 0.30190).
- `k3::Number`: Degree-3 Love number (default: 0.093).
- `k2_plus::Number`: Degree-2 → degree-4 coupling Love number (default: 0.0).
- `R_e::Number`: Central body equatorial radius [km]. No default — must be specified.
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
    R_e::RET
    include_degree_3::BT = true
end

"""
    acceleration(u, p, t, model::SolidBodyTidesModel) -> SVector{3}

Compute the solid body tides acceleration perturbation.

# Arguments
- `u::AbstractVector`: Spacecraft state [km, km/s].
- `p::FrameAwareParams`: Parameters with frame system.
- `t::Number`: Elapsed time since epoch [s].
- `model::SolidBodyTidesModel`: Solid body tides model.

# Returns
- `SVector{3}`: Tidal acceleration perturbation [km/s²].
"""
@inline function acceleration(
    u::AbstractVector, p::FrameAwareParams, t::Number, model::SolidBodyTidesModel
)
    t_ft = ft_time(p, t)
    bodies = _collect_body_data(p.frames, t_ft, model.tide_raising_bodies)
    return solid_body_tides_accel(
        u,
        bodies;
        k2=model.k2,
        k3=model.k3,
        R_e=model.R_e,
        include_degree_3=model.include_degree_3,
    )
end

@inline function _collect_body_data(frames, t_ft, bodies::Tuple)
    body = first(bodies)
    r = get_position(body.ephem_type, body.body, frames, t_ft, body.compiled_vector3)
    μ = body.body.μ
    return ((r, μ), _collect_body_data(frames, t_ft, Base.tail(bodies))...)
end

@inline _collect_body_data(frames, t_ft, ::Tuple{}) = ()

"""
    solid_body_tides_accel(u, bodies; k2, k3, R_e, include_degree_3) -> SVector{3}

Compute the solid body tides acceleration from tide-raising body data.

# Arguments
- `u::AbstractVector`: Spacecraft state [km, km/s].
- `bodies::Tuple`: Tuple of `(r_body, μ_body)` pairs [km, km³/s²].

# Keyword Arguments
- `k2::Number`: Degree-2 Love number. Default: `0.30190`.
- `k3::Number`: Degree-3 Love number. Default: `0.093`.
- `R_e::Number`: Central body equatorial radius [km]. **Required — no default**
  (prevents silent misuse for non-Earth missions).
- `include_degree_3::Bool`: Include degree-3 contribution. Default: `true`.

# Returns
- `SVector{3}`: Tidal acceleration [km/s²].
"""
@inline function solid_body_tides_accel(
    u::AbstractVector{UT},
    bodies::Tuple;
    k2::K2T=0.30190,
    k3::K3T=0.093,
    R_e::RET,
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
