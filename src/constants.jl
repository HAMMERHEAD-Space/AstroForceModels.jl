# Body constants (R_SUN, R_EARTH, R_MOON, μ_SUN, μ_EARTH, μ_MOON, etc.)
# are provided by CelestialBodies.jl — do not redefine here.

export SPEED_OF_LIGHT,
    SOLAR_FLUX,
    SOLAR_IRRADIANCE,
    EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS,
    MASS_O_PLUS

# Speed of Light [km/s]
const SPEED_OF_LIGHT::Float64 = 2.99792458E5
# Solar Irradiance at 1 AU [W/m^2]
const SOLAR_IRRADIANCE::Float64 = 1360.8
# Solar Flux [N/m^2] (solar irradiance / speed of light in m/s)
const SOLAR_FLUX::Float64 = SOLAR_IRRADIANCE / (SPEED_OF_LIGHT * 1E3)
# Earth's Specific Angular Momentum (angular momentum per unit mass) [km^2/s]
const EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS::Float64 = 0.4 * R_EARTH^2 * EARTH_ANGULAR_SPEED
# Mass of O⁺ ion [kg] — dominant ionospheric ion in the F2 region (200-600 km)
const MASS_O_PLUS::Float64 = 15.999 * 1.66053906660e-27
