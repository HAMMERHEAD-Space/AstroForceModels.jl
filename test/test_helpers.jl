# Test helper: creates a FrameSystem with analytical Vallado ephemeris for Sun/Moon
# and wraps parameters in FrameAwareParams for use in all tests.

import SatelliteToolboxTransformations as STT
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
    frames = FrameSystem{2, Float64}()
    add_axes_icrf!(frames)

    # Add ITRF rotating frame using EOP data
    add_axes_rotating!(frames, :ITRF, 2, 1,
        t -> STT.r_eci_to_ecef(STT.DCM, STT.J2000(), STT.ITRF(), 2451545.0 + t / 86400.0, eop_data),
        t -> begin
            jd = 2451545.0 + t / 86400.0
            dt = 0.01
            R1 = STT.r_eci_to_ecef(STT.DCM, STT.J2000(), STT.ITRF(), jd, eop_data)
            R2 = STT.r_eci_to_ecef(STT.DCM, STT.J2000(), STT.ITRF(), jd + dt/86400.0, eop_data)
            return (R2 - R1) / dt
        end,
    )

    # Earth as root point
    add_point!(frames, :Earth, 399, :ICRF)

    # Sun relative to Earth via Vallado analytical ephemeris
    function sun_state(t)
        jd = 2451545.0 + t / 86400.0
        pos_mod = sun_position_mod(jd)
        vel_mod = sun_velocity_mod(jd)
        R = STT.r_eci_to_eci(STT.MOD(), STT.J2000(), jd)
        pos_j2000 = R * pos_mod ./ 1e3  # m -> km
        vel_j2000 = R * vel_mod ./ 1e3  # m/s -> km/s
        return vcat(pos_j2000, vel_j2000)
    end
    add_point_dynamical!(frames, :Sun, 10, 399, :ICRF, sun_state)

    # Moon relative to Earth via Vallado analytical ephemeris
    function moon_state(t)
        jd = 2451545.0 + t / 86400.0
        pos_mod = moon_position_mod(jd)
        R = STT.r_eci_to_eci(STT.MOD(), STT.J2000(), jd)
        pos_j2000 = R * pos_mod ./ 1e3
        # Finite-difference velocity (Vallado doesn't provide moon velocity)
        dt = 1.0
        jd2 = jd + dt / 86400.0
        pos_mod2 = moon_position_mod(jd2)
        R2 = STT.r_eci_to_eci(STT.MOD(), STT.J2000(), jd2)
        pos2 = R2 * pos_mod2 ./ 1e3
        vel_j2000 = (pos2 - pos_j2000) / dt
        return vcat(pos_j2000, vel_j2000)
    end
    add_point_dynamical!(frames, :Moon, 301, 399, :ICRF, moon_state)

    return frames
end

"""
    create_test_params(; JD, eop_data, extra_params...)

Create FrameAwareParams for testing with the standard test FrameSystem.
Returns (params, frames) tuple so frames can be reused for compiled transforms.
"""
function create_test_params(; JD, eop_data, extra_params...)
    frames = create_test_frames(eop_data)
    base = ComponentVector(; JD=JD, extra_params...)
    epoch = Epoch((JD - 2451545.0) * 86400.0, TDB)
    return FrameAwareParams(base, frames, epoch, :ICRF)
end

# Standard third-body models for tests (Earth-centric, ICRF)
function test_sun_model(; frames=nothing)
    ThirdBodyModel(
        body=SunBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=10, axes=:ICRF),
        frames=frames,
    )
end

function test_moon_model(; frames=nothing)
    ThirdBodyModel(
        body=MoonBody(),
        ephem_type=FrameEphemeris(center_point=399, target_point=301, axes=:ICRF),
        frames=frames,
    )
end
