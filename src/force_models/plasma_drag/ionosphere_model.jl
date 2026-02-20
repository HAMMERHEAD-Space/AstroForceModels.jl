# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Ionospheric ion density models for plasma drag computation.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Lafleur, T. (2023). "Charged aerodynamics: Ionospheric plasma drag on objects in
#       low-Earth orbit." Acta Astronautica, 212, 370-386.
#
#   [2] Rishbeth, H. & Garriott, O.K. (1969). "Introduction to Ionospheric Physics."
#       Academic Press.
#
#   [3] Bilitza, D. (2018). "IRI the International Standard for the Ionosphere."
#       Advances in Radio Science, 16, 1-11.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export AbstractIonosphereModel, ChapmanIonosphere, ConstantIonosphere, NoIonosphere
export compute_ion_density

"""
    AbstractIonosphereModel

Abstract type for ionospheric density models used to compute ion mass density for 
plasma drag calculations.

# Implementations
- [`ChapmanIonosphere`](@ref): Chapman layer model for the F2 region
- [`ConstantIonosphere`](@ref): Fixed ion density for testing and parametric studies
- [`NoIonosphere`](@ref): Zero ion density (disables plasma drag)
"""
abstract type AbstractIonosphereModel end

"""
    ChapmanIonosphere{NmT,HmT,HsT,MiT} <: AbstractIonosphereModel

Chapman layer model for the F2 region of the ionosphere.

The ion number density follows the Chapman production function profile:

    nᵢ(h) = Nₘₐₓ ⋅ exp(½(1 - z - exp(-z)))

where z = (h - hₘₐₓ) / Hₛ, and the ion mass density is ρᵢ = nᵢ ⋅ mᵢ.

Default parameters represent typical mid-latitude, moderate solar activity conditions. 
The dominant ion in the F2 region (200-500 km) is O⁺.

# Fields
- `Nmax::NmT`: Peak ion number density [m⁻³]. Typical range: 10¹⁰ (solar min) to 10¹² (solar max).
- `hmax::HmT`: Altitude of the F2 peak [km]. Typical range: 250-450 km.
- `Hs::HsT`: Scale height of the F2 layer [km]. Typical range: 50-80 km.
- `mi::MiT`: Mean ion mass [kg]. Default is O⁺ mass (2.6567628e-26 kg).

# References
- Rishbeth & Garriott (1969), "Introduction to Ionospheric Physics"
- Lafleur (2023), "Charged aerodynamics", Acta Astronautica 212, 370-386
"""
Base.@kwdef struct ChapmanIonosphere{NmT<:Number,HmT<:Number,HsT<:Number,MiT<:Number} <:
                   AbstractIonosphereModel
    Nmax::NmT = 3e11
    hmax::HmT = 350.0
    Hs::HsT = 60.0
    mi::MiT = MASS_O_PLUS
end

"""
    ConstantIonosphere{DT} <: AbstractIonosphereModel

Ionosphere model with a fixed ion mass density, useful for testing and parametric studies.

# Fields
- `rho_i::DT`: Fixed ion mass density [kg/m³].

# Example
```julia
iono = ConstantIonosphere(rho_i=1e-17)
```
"""
Base.@kwdef struct ConstantIonosphere{DT<:Number} <: AbstractIonosphereModel
    rho_i::DT
end

"""
    NoIonosphere <: AbstractIonosphereModel

Null ionosphere model that returns zero density. Use to disable plasma drag while 
keeping the model in the perturbation tuple for interface consistency.
"""
struct NoIonosphere <: AbstractIonosphereModel end

"""
    compute_ion_density(JD, u, eop_data, model::ChapmanIonosphere)

Compute the ion mass density at the spacecraft position using a Chapman layer profile.

The geodetic altitude is computed from the J2000 state vector via ECEF transformation,
and the Chapman function is evaluated for the F2 layer. Density is zero above 1500 km
where the ionosphere becomes negligible.

# Arguments
- `JD::Number`: Current Julian Date.
- `u::AbstractVector`: Spacecraft state vector [r; v] in J2000 ECI [km; km/s].
- `eop_data`: Earth Orientation Parameters for coordinate transformations.
- `model::ChapmanIonosphere`: Chapman layer model parameters.

# Returns
- `rho_i::Number`: Ion mass density [kg/m³].
"""
@inline function compute_ion_density(
    JD::Number, u::AbstractVector, eop_data::EopIau1980, model::ChapmanIonosphere
)
    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), JD, eop_data)
    ecef_pos = R_J20002ITRF * SVector{3}(u[1], u[2], u[3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    h_km = geodetic_pos[3] / 1E3

    z = (h_km - model.hmax) / model.Hs
    n_i = model.Nmax * exp(0.5 * (1.0 - z - exp(-z)))

    # Zero above 1500 km where ionosphere is negligible
    return (h_km < 1500.0) * n_i * model.mi
end

@inline function compute_ion_density(
    JD::Number, u::AbstractVector, eop_data::EopIau1980, model::ConstantIonosphere
)
    return model.rho_i
end

@inline function compute_ion_density(
    JD::Number, u::AbstractVector, eop_data::EopIau1980, model::NoIonosphere
)
    return zero(eltype(u))
end
