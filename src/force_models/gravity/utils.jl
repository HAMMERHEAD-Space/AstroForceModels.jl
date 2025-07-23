"""
    err_iau2006(JD_TT::Number) -> Float64

Compute the Earth Rotation Rate (time derivative of the Greenwich Mean Sidereal Time)
according to the IAU 2006 Conventions.

# Args

- `JD_TT`: Terrestrial Time (TT) Julian date.

# Returns

- `Float64`: The Earth Rotation Rate [rad/s].
"""
function err_iau2006(JD_TT::Number)
    # Coefficients for the rate of change of GMST.
    c0 = 4612.156534
    c1 = 2.7831634
    c2 = - 0.00000132
    c3 = - 0.000119824
    c4 = - 0.000000184

    # Arcseconds to radians.
    as2rad = 4.848136811095359935899141E-6

    # Seconds to centuries.
    sec2cy = 3155760000.0

    # The Earth Rotation Angle (ERA) is a linear function of UT1, thus its derivative
    # is constant.
    dotERA = 7.292115146706980494985261831431E-5 # [rad/s]

    # Following SOFA, we compute the TT Julian centuries since J2000 as
    T = (JD_TT - 2451545.0) / 36525.0

    dotDelta = @evalpoly(T, c0, c1, c2, c3, c4) # [as/cy]
    dotDelta = dotDelta * as2rad / sec2cy       # [rad/s]

    return dotERA + dotDelta
end
