export R_SUN,
    R_EARTH,
    R_MOON,
    SPEED_OF_LIGHT,
    SOLAR_FLUX,
    SOLAR_IRRADIANCE,
    μ_MOON,
    μ_SUN,
    μ_EARTH,
    EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS

# Radius of the Sun [km]
const R_SUN::Float64 = 6.955E5
# Radius of the Earth [km]
const R_EARTH::Float64 = 6378.1363
# Radius of the Moon [km]
const R_MOON::Float64 = 1738.1
# Speed of Light [km/s]
const SPEED_OF_LIGHT::Float64 = 2.99792458E5
# Solar Irradiance at 1 AU [W/m^2]
const SOLAR_IRRADIANCE::Float64 = 1360.8
# Solar Flux [N/m^2] (solar irradiance / speed of light in m/s)
const SOLAR_FLUX::Float64 = SOLAR_IRRADIANCE / (SPEED_OF_LIGHT * 1E3)
# Gravitational Parameter of the Moon [km^3/s^2]
const μ_MOON::Float64 = 4.902800118457551E3
# Gravitational Parameter of the Sun [km^3/s^2]
const μ_SUN::Float64 = 1.3271244004127946E11
# Gravitational Parameter of the Earth [km^3/s^2]
const μ_EARTH::Float64 = 3.986004415e5
# Earth's Specific Angular Momentum (angular momentum per unit mass) [km^2/s]
const EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS::Float64 = 0.4 * R_EARTH^2 * EARTH_ANGULAR_SPEED
