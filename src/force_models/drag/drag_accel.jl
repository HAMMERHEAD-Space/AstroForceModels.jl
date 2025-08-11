# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Drag
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#TODO: REFERENCE
#   [1] 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export DragAstroModel, drag_accel

"""
    DragAstroModel{ST,AT,EoT,RT,PT} <: AbstractNonPotentialBasedForce

Atmospheric drag force model for spacecraft orbital dynamics.

This model computes the acceleration due to atmospheric drag acting on a spacecraft in 
Earth's atmosphere. The drag force is proportional to the atmospheric density and the 
square of the relative velocity between the spacecraft and the atmosphere.

# Type Parameters
- `ST <: AbstractSatelliteDragModel`: Type of the satellite drag model
- `AT <: AtmosphericModelType`: Type of the atmospheric model
- `EoT <: Union{EopIau1980,EopIau2000A}`: Type of Earth Orientation Parameters
- `RT <: Union{Nothing,AbstractVector}`: Type for RTS (optional)
- `PT <: Union{Nothing,AbstractMatrix}`: Type for P matrix (optional)

# Fields
- `satellite_drag_model::ST`: The satellite drag model for computing ballistic coefficient (Cd*A/m)
- `atmosphere_model::AT`: The atmospheric model for computing density (e.g., JR1971, JB2008, NRLMSISE00)
- `eop_data::EoT`: Earth Orientation Parameters for coordinate transformations and atmospheric modeling
- `rts::RT`: Optional RTS parameter for atmospheric model (defaults to nothing)
- `P::PT`: Optional P matrix parameter for atmospheric model (defaults to nothing)

# Example
```julia
# Create satellite drag model
sat_drag = CannonballDragModel(
    area = 10.0,        # [m²]
    drag_coeff = 2.2,   # dimensionless
    mass = 1000.0       # [kg]
)

# Create drag force model
drag_model = DragAstroModel(
    satellite_drag_model = sat_drag,
    atmosphere_model = JB2008(),
    eop_data = eop_data
)
```

# See Also
- [`acceleration`](@ref): Compute drag acceleration
- [`AbstractSatelliteDragModel`](@ref): Satellite drag model interface
- [`density_calculator`](@ref): Atmospheric density computation
"""
@with_kw struct DragAstroModel{ST,AT,EoT,RT,PT} <: AbstractNonPotentialBasedForce where {
    ST<:AbstractSatelliteDragModel,
    AT<:AtmosphericModelType,
    EoT<:Union{EopIau1980,EopIau2000A},
    RT<:Union{Nothing,AbstractVector},
    PT<:Union{Nothing,AbstractMatrix},
}
    satellite_drag_model::ST
    atmosphere_model::AT
    eop_data::EoT

    rts::RT = nothing
    P::PT = nothing
end

"""
    acceleration(u::AbstractVector, p::ComponentVector, t::Number, drag_model::DragAstroModel)

Computes the drag acceleration acting on a spacecraft given a drag model and current state and 
parameters of an object.

# Arguments
- `u::AbstractVector`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `drag_model::DragAstroModel`: Drag model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional drag acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractVector, p::ComponentVector, t::Number, drag_model::DragAstroModel
)
    # Compute density at the satellite's current position
    rho = compute_density(
        p.JD + t / 86400.0,
        u,
        drag_model.eop_data,
        drag_model.atmosphere_model;
        roots_container=drag_model.rts,
        P=drag_model.P,
    )

    #TODO: OFFER OPTION TO COMPUTE FROM EOP or SPICE EPHEMERIS 
    ω_vec = SVector{3,Float64}(0.0, 0.0, EARTH_ANGULAR_SPEED)

    # Compute the ballistic coefficient
    BC = ballistic_coefficient(u, p, t, drag_model.satellite_drag_model)

    # Return the 3-Dimensional Drag Force
    return drag_accel(u, rho, BC, ω_vec)
end

"""
    drag_accel(u::AbstractVector, rho::Number, BC::Number, ω_vec::AbstractVector, t::Number, [DragModel]) -> SVector{3}{Number}

Compute the Acceleration Atmospheric Drag

The atmosphere is treated as a solid revolving with the Earth and the apparent velocity of the satellite is computed
using the transport theorem

                𝐯_app = 𝐯 - 𝛚 x 𝐫

The acceleration from drag is then computed with a cannonball model as
                
                𝐚 = 1/2 * ρ * BC * |𝐯_app|₂^2 * v̂


!!! note
    Currently only fixed cannonball state based ballistic coefficients are supported, custom models can be created for
    higher fidelity.

# Arguments

- `u::AbstractVector`: The current state of the spacecraft in the central body's inertial frame.
- `rho::Number`: Atmospheric density at (t, u) [kg/m^3].
- `BC::Number`: The ballistic coefficient of the satellite -- (area/mass) * drag coefficient [kg/m^2].
- `ω_vec::AbstractVector`: The angular velocity vector of Earth. Typically approximated as [0.0; 0.0; ω_Earth]
- `t::Number`: Current time of the simulation.

# Returns

- `SVector{3}{Number}`: Inertial acceleration from drag
"""
@inline function drag_accel(
    u::AbstractVector{UT}, rho::RT, BC::BT, ω_vec::AbstractVector{WT}
) where {UT,RT,BT,WT}
    AT = promote_type(UT, RT, BT, WT)

    # Compute Apparent Velocity w.r.t the Atmosphere using the Transport Theorem
    apparent_vel = SVector{3}(u[4], u[5], u[6]) - cross(ω_vec, SVector{3}(u[1], u[2], u[3]))

    # Scaled by 1E3 to convert to km/s
    drag_force = -0.5 * BC * rho * norm(apparent_vel) * 1E3
    # TODO: HANDLE UNITS BETTER
    accel = SVector{3,AT}(
        drag_force * apparent_vel[1],
        drag_force * apparent_vel[2],
        drag_force * apparent_vel[3],
    )

    return accel
end
