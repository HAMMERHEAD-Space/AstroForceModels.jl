using AllocCheck
using Aqua
using AstroForceModels
using ComponentArrays
using JET
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using StaticArraysCore
using Test

using DifferentiationInterface
using FiniteDiff, ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote

@testset "AstroForceModels.jl" begin
    # Drag Tests
    include("drag/test_satellite_shape_model.jl")
    include("drag/test_density_calculator.jl")
    include("drag/test_drag_accel.jl")

    # SRP Tests
    include("solar_radiation_pressure/test_satellite_shape_models.jl")
    include("solar_radiation_pressure/test_shadow_models.jl")
    include("solar_radiation_pressure/test_srp_accel.jl")

    # Third Body Tests
    include("third_body/test_celestial_body.jl")
    include("third_body/test_third_body_model.jl")
    include("third_body/test_third_body_accel.jl")

    # Zonal Harmonics Tests
    include("gravity/test_grav_accel.jl")

    # Relativity Tests
    include("relativity/test_relativity.jl")

    # Dynamics Builder
    include("test_dynamics_builder.jl")
end

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

@testset "Differentiability" begin
    include("differentiability/test_model_parameters.jl")
    include("differentiability/test_drag.jl")
    include("differentiability/test_srp.jl")
    include("differentiability/test_gravity.jl")
    include("differentiability/test_relativity.jl")
    include("differentiability/test_third_body.jl")
    include("differentiability/test_dynamics_builder.jl")
end

@testset "Performance" begin
    # Force Model Allocation Check
    include("test_allocations.jl")
    include("test_JET.jl")
    
    Aqua.test_all(AstroForceModels; ambiguities=(recursive = false))
end
