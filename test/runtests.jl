using AstroForceModels
using ComponentArrays
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using StaticArraysCore
using Test

using AllocCheck
using Aqua
using JET

@testset "AstroForceModels.jl" begin
    # Drag Tests
    include("drag/test_satellite_shape_model.jl")
    include("drag/test_density_calculator.jl")
    include("drag/test_drag_accel.jl")

    # SRP Tests
    include("solar_radiation_pressure/test_satellite_shape_models.jl")
    include("solar_radiation_pressure/test_shadow_models.jl")
    include("solar_radiation_pressure/test_srp_accel.jl")
    include("solar_radiation_pressure/test_albedo_accel.jl")

    # Third Body Tests
    include("third_body/test_celestial_body.jl")
    include("third_body/test_third_body_model.jl")
    include("third_body/test_third_body_accel.jl")

    # Zonal Harmonics Tests
    include("gravity/test_grav_accel.jl")

    # Relativity Tests
    include("relativity/test_relativity.jl")

    # Solid Body Tides Tests
    include("solid_body_tides/test_solid_body_tides.jl")

    # Thermal Emission Tests
    include("thermal_emission/test_thermal_emission_accel.jl")

    # Low Thrust Tests
    include("low_thrust/test_thrust_model.jl")
    include("low_thrust/test_low_thrust_accel.jl")

    # Dynamics Builder
    include("test_dynamics_builder.jl")
end

# Differentiability tests are gated behind an environment variable to keep the default
# test suite fast. Set ASTROFORCEMODELS_TEST_DIFF to run them:
#   "true"  / "all"       → all 5 backends (ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote)
#   "ForwardDiff"         → ForwardDiff only (fast smoke-test for AD compatibility)
#   unset   / "false"     → skip differentiability tests entirely
const _DIFF_ENV = get(ENV, "ASTROFORCEMODELS_TEST_DIFF", "false")

if _DIFF_ENV ∉ ("false", "")
    using DifferentiationInterface
    using FiniteDiff

    _run_all = _DIFF_ENV ∈ ("true", "all")
    _requested = _run_all ? Set{String}() : Set(strip.(split(_DIFF_ENV, ",")))
    _need(name) = _run_all || name ∈ _requested

    # Only load and instantiate backends that are actually requested.
    # This avoids compiling Enzyme/Mooncake/Zygote when only ForwardDiff is needed.
    _backend_list = Tuple{String,Any}[]

    if _need("ForwardDiff")
        using ForwardDiff
        push!(_backend_list, ("ForwardDiff", AutoForwardDiff()))
    end
    if _need("Enzyme")
        using Enzyme
        push!(
            _backend_list,
            ("Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))),
        )
    end
    if _need("Mooncake")
        using Mooncake
        push!(_backend_list, ("Mooncake", AutoMooncake(; config=nothing)))
    end
    if _need("PolyesterForwardDiff")
        using PolyesterForwardDiff
        push!(_backend_list, ("PolyesterForwardDiff", AutoPolyesterForwardDiff()))
    end
    if _need("Zygote")
        using Zygote
        push!(_backend_list, ("Zygote", AutoZygote()))
    end

    if isempty(_backend_list)
        error(
            "ASTROFORCEMODELS_TEST_DIFF=\"$_DIFF_ENV\" did not match any backend. " *
            "Valid names: ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote",
        )
    end

    const _BACKENDS = Tuple(_backend_list)

    @info "Running differentiability tests with backends: $(join([b[1] for b in _BACKENDS], ", "))"

    @testset "Differentiability" begin
        include("differentiability/test_model_parameters.jl")
        include("differentiability/test_drag.jl")
        include("differentiability/test_srp.jl")
        include("differentiability/test_gravity.jl")
        include("differentiability/test_relativity.jl")
        include("differentiability/test_third_body.jl")
        include("differentiability/test_low_thrust.jl")
        include("differentiability/test_albedo.jl")
        include("differentiability/test_solid_body_tides.jl")
        include("differentiability/test_thermal_emission.jl")
        include("differentiability/test_dynamics_builder.jl")
    end
else
    @info "Skipping differentiability tests (set ASTROFORCEMODELS_TEST_DIFF to enable)"
end

@testset "Performance" begin
    # Force Model Allocation Check
    include("test_allocations.jl")
    include("test_JET.jl")

    @testset "Aqua.jl" begin
        Aqua.test_all(
            AstroForceModels;
            ambiguities=(recursive = false),
            deps_compat=(check_extras = false),
        )
    end
end
