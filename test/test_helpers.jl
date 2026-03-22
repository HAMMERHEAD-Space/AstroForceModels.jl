# Test helper: creates a FrameSystem with analytical Vallado ephemeris for Sun/Moon
# and wraps parameters in FrameAwareParams for use in all tests.

using SatelliteToolboxCelestialBodies
using FrameTransformations
using Tempo
using ComponentArrays

"""
    create_test_frames(eop_data)

Create a FrameSystem with ICRF axes, ITRF axes, and Earth/Sun/Moon points
using Vallado analytical ephemeris wrapped as dynamical points.
"""
function create_test_frames(eop_data)
    epoch = Epoch(0.0, TDB)  # dummy epoch, only need frames
    p = setup_earth_propagation_frames(epoch, eop_data)
    return p.frames
end

"""
    create_test_params(; JD, eop_data, extra_params...)

Create FrameAwareParams for testing with the standard test FrameSystem.
"""
function create_test_params(; JD, eop_data, extra_params...)
    frames = create_test_frames(eop_data)
    epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
    params = isempty(extra_params) ? ComponentVector() : ComponentVector(; extra_params...)
    return FrameAwareParams(params, frames, epoch, :ICRF)
end

# Standard third-body models for tests (Earth-centric, ICRF)
function test_sun_model(; frames=nothing)
    ThirdBodyModel(;
        body=SunBody(),
        ephem_type=FrameEphemeris(; center_point=399, target_point=10, axes=:ICRF),
        frames=frames,
    )
end

function test_moon_model(; frames=nothing)
    ThirdBodyModel(;
        body=MoonBody(),
        ephem_type=FrameEphemeris(; center_point=399, target_point=301, axes=:ICRF),
        frames=frames,
    )
end
