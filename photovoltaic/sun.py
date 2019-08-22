from .core import q, k, h, c, Stefan_Boltzmann, cosd, sind, arcsind, arccosd, Wien
from scipy import integrate
import numpy as np
import os

def am_intensity(airmass):
    """Return radiation intensity (W/m**2) given airmass (units) """
    It = 1.353 * (0.7 ** (airmass ** 0.678))
    Ig = It * 1.1
    return It, Ig


def air_mass(angle):
    """Return air mass (units) where *angle* is the zenith angle (degrees) """
    return 1 / cosd(angle)


def am_shadow(s, h):
    """Return the air mass (units) where h is the height of a pole and s is the length of its shadow.
     s and h are the same length units, i.e. both in m, ft, cm etc."""
    am = np.sqrt(1 + (s / h) ** 2)
    return am


def blackbody_peak(T):
    """Return the peak wavelength (nm) of a black body spectrum where T is the temperature (K)"""
    return 1e9 * Wien / T


def blackbody_integrated(T):
    """Return integrated irradiance (W/m2/steradian) from a blackbody where T is the temperature (K)"""
    return Stefan_Boltzmann * T ** 4


def blackbody_spectrum(wavelength, T=6000):
    """Return the blackbody irradaiance (W/m2/nm) at a given wavelength (nm) and temperature, T (K). """
    wavelength = wavelength * 1e-9
    F = 2 * π * h * c ** 2 / ((wavelength ** 5) * (np.exp(h * c / (wavelength * T * k)) - 1))
    return F * 1e-9  # convert to W/m3 to W/m2/nm


def equal_spacing(x, y, x_min, x_max, x_step):
    """Returns spectra with equal spacking and truncation (W/m2) (NOT W/m2/nm)
    given a spectrum (W/m2/nm) as a function of wavelength (nm)
    wavelength minimum (nm) and wavlength maximum (nm)
    Note: would actually work for any x y data"""
    y_integrated = integrate.cumtrapz(y, x, initial=0)
    x_midpoint = np.arange(x_min, x_max + x_step, x_step)
    x_extents = np.arange(x_min - x_step / 2, x_max + 3 / 2 * x_step, x_step)
    y_spaced = np.diff(np.interp(x_extents, x, y_integrated))
    return x_midpoint, y_spaced


def space_solar_power(x):
    """Return the radiant power density (W/m²) where x is the distance from  the sun (m)"""
    return 2.8942e25 / (x ** 2)


def solar_spectra(fname=None):
    """Return wavelength (nm) and AM0, AM15G, AM15D (W/m²/nm)
    of the standard spectrum.
    reference XXX, DOI XXX"""
    if fname is None:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'ASTMG173.txt')
    wavelength, AM0, AM15G, AM15D = np.genfromtxt(fname, skip_header=2, unpack=True)
    return wavelength, AM0, AM15G, AM15D


def etr_earth(day_no):
    """Return extraterrestrial radiation at earth (W/m**2) where day_no is the day of the year (day).
     January 1 has a day_no of 1. There is no correction for leap year."""
    return (1 + .033 * (cosd((day_no - 2.0) * (360.0 / 365.0))) * 1353)


def declination(day_no):
    """Return declination angle of sun (degrees) where day_no is the day of the year (day).
    For Jan 1 day_no = 1, Dec 31 dayno = 365. There is no correction for leap years"""
    return 23.45 * sind((day_no - 81) * (360 / 365))


def equation_of_time(day_no):
    """Return the equation of time (minutes) where day_no is the day of the year (day). """
    B = 360.0 / 365.0 * (day_no - 81.0)
    EoT = 9.87 * sind(2 * B) - 7.53 * cosd(B) - 1.5 * sind(B)
    # print('EoT', EoT)
    return EoT


def time_correction(EoT, longitude, GMTOffset):
    """Return the time correction (minutes) where EoT is the equation of time,
    and given location longitude (degrees) and the GMT offset (hours)"""
    LSTM = 15.0 * GMTOffset
    TimeCorrection = 4.0 * (longitude - LSTM) + EoT
    return TimeCorrection


def elevation(declination, latitude, local_solar_time):
    """Return the elevation angle of the sun (degrees)
    given declination (degrees), latitude (degrees) and local_solar_time (hours) """
    hra = 15.0 * (local_solar_time - 12.0)
    return arcsind(sind(declination) * sind(latitude) + cosd(declination) * cosd(latitude) * cosd(hra))


def sun_rise_set(latitude, declination, time_correction):
    """Return the sunrise and sunset times in hours
    given the latitude (degrees) and the declination (degrees)
    """
    A = -1 * (sind(latitude) * sind(declination)) / (cosd(latitude) * cosd(declination))
    local_solar_time = arccosd(A) / 15.0
    sunrise = 12.0 - local_solar_time - (time_correction / 60.0)
    sunset = 12 + local_solar_time - (time_correction / 60.0)
    return sunrise, sunset


def elev_azi(declination, latitude, local_solar_time):
    """Return the elevation (degrees) and azimuth (degrees)"""
    hour_angle = 15.0 * (local_solar_time - 12.0)
    elevation = arcsind(sind(declination) * sind(latitude) + cosd(declination) * cosd(latitude) * cosd(hour_angle))
    azimuth = arccosd(
        (cosd(latitude) * sind(declination) - cosd(declination) * sind(latitude) * cosd(hour_angle)) / cosd(elevation))
    # the multiplication by 1.0 causes a single value return for single inputs, otherwise it returns an array of one element
    azimuth = np.where(hour_angle > 0, 360.0 - azimuth, azimuth) * 1.0
    return elevation, azimuth


def module_direct(azimuth, elevation, module_azimuth, module_tilt):
    """Returns the faction of light on a arbtrarily tilted surface
     given sun azimuth (degrees) where north is zero and elevation
     module_azimuth and module_tilt, where """
    fraction = cosd(elevation) * sind(module_tilt) * cosd(azimuth - module_azimuth) + sind(elevation) * cosd(
        module_tilt)
    return fraction


def sun_position(dayNo, latitude, longitude, GMTOffset, H, M):
    """Return the position of the sun as a elevation and azimuth given
    latitude, logitude and the GMTOffset, """
    EoT = equation_of_time(dayNo)
    TimeCorrection = time_correction(EoT, longitude, GMTOffset)
    local_solar_time = H + (TimeCorrection + M) / 60.0
    elevation, azimuth = elev_azi(declination(dayNo), latitude, local_solar_time)
    return elevation, azimuth



