#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Basic photovoltaic functions
Provides the formulas and equations typically used in an introductory photocoltaics textbook.

Typical solar units are used, NOT SI units. The units are denoted in parenthesis on the comment lines.
wavelength (nm)
Energy of  photon (eV)
semiconductor dimensions (cm)
degrees instead of radians.
Temperature of 298.15 K (25 degC) not 300 K

The first line on all input files is ignored to allow for column headers
# denotes a comment in input files and is ignored.

Contributions by: sgbowden, richter, heraimenka, jhul etc

Variables and acronyms
ARC - anti-reflection coating
"""
__version__ = '0.1.3'

import numpy as np
import os
"""from scipy import integrate"""

# define constants
q = 1.60217662e-19  # (coulombs) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
k_eV = 8.6173303e-05  # (eV K^-1)
Wien = 2.898e-3  # (m K)
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
π = np.pi  # yes, I use unicode
pi = np.pi  # compatibility with class
h = 6.62607004e-34  # (J.s)
hbar = 6.62607004e-34 / (2 * π)  # usable
c = 299792458.0  # (m s^-1)
hc_q = h * c / q  # 1.2398419745831506e-06


# ******** Section: helpers *********
def sind(angle):
    """Return the sine of the angle(degrees)
    Example:
    >>>sind(0)
    0
    """
    return np.sin(np.radians(angle))


def cosd(angle):
    """Return the cosine of the angle(degrees)"""
    return np.cos(np.radians(angle))


def tand(angle):
    """Return the tangent of the angle(degrees)"""
    return np.tan(np.radians(angle))


def arcsind(x):
    """Return the arcsin (degrees)"""
    return np.degrees(np.arcsin(x))


def arccosd(x):
    """Return the arccos (degrees)"""
    return np.degrees(np.arccos(x))


def module_path():
    # no longer used as I put it in the individual functions
    import os
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    return dir_path


# ******** Section: Solar Radiation *********
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
    """Return the blackbody irradaiance (W/nm) at a given wavelength (nm) and temperature, T (K). """
    wavelength = wavelength * 1e-9
    F = 2 * π * h * c ** 2 / ((wavelength ** 5) * (np.exp(h * c / (wavelength * T * k)) - 1))
    return F * 1e-9  # convert to W/m to W/nm


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


def spectrum_spacing(x, y, x_min, x_max, x_step):
    """Returns spectra in W/m2 (NOT W/m2/nm) with equal spacing
    given a spectrum (W/m2/nm) as a function of wavelength (nm)
    wavelength minimum (nm) and wavlength maximum (nm) """
    y_integrated = integrate.cumtrapz(y, x, initial=0)
    x_midpoint = np.arange(x_min, x_max + x_step, x_step)
    x_extents = np.arange(x_min - x_step / 2, x_max + 3 / 2 * x_step, x_step)
    y_spaced = np.diff(np.interp(x_extents, x, y_integrated))
    return x_midpoint, y_spaced


# ******** Section: Light *********

def nm2eV(x):
    """ Given wavelength (nm) of a photon return the energy (eV) """
    return hc_q * 1e9 / x


def eV2nm(x):
    """ Given energy (eV) of a photon return the wavelength (nm) """
    return hc_q * 1e9 / x


def nm2joule(x):
    """ Given wavelength (nm) of a photon return the energy (eV) """
    return h * c * 1e9 / x


def photon_flux(power, wavelength):
    """Return the photon flux (/s) given the power of light (watts) and wavelength (nm)
    If power is in W/m2 then flux is in m-2s-1"""
    return power / nm2joule(wavelength)


def k2alpha(kd, wavelength):
    """Quick convert of extinction coefficient (units) to absorption coefficient(cm-1) given the wavelength (nm)"""
    return 1e7 * 4 * π * kd / wavelength


def refraction(n1, n2, θ1):
    """Return the refracted angle of light (degrees) where n1 is the refractive index of incident medium (units),
    n2 is the refractive index of the transmission medium (units) and θ1 is the incident angle to the normal"""
    θ2 = arcsind(n1 / n2 * sind(θ1))
    return θ2


def absorption_coeff(kd, wavelength):
    """absorption coefficient (cm-1) from extinction coefficient (units) and wavelength (nm)
    wavelength """
    return 1e7 * 4 * π * kd / wavelength


def transmittance(abs_coeff, thickness):
    """Return the fraction of light transmitted (units) where abs_coeff is the absorption coefficient (cm-1)
    and 'thickness' is the depth in the material (cm)
    """
    return np.exp(-abs_coeff * thickness)


def ARC_thick(wavelength, n1):
    """Return optimal anti-reflection coating thickness (nm) at a given wavelength (nm)
    where n1 is the refractive index of the antireflection coating.
    The returned unit is in the same as the wavelength unit."""
    return wavelength / (4 * n1)


def ARC_opt_n(n0, n2):
    """Return the optimal refractive index, typically denoted as n1, for an antireflection coating (units)
    where n0 is the refractive index of the incident medium and n2 is the refractive index of the object (units)"""
    return np.sqrt(n0 * n2)


def ARC_refl(wavelength, n0, n1, nSemi, thickness):
    """Return the reflectivity from a object (units) that has an anti-reflection coating.
    Where:
    n0 - the ambient refractive index units),
    n1 - refractive index of the dielectric layer (units)
    nSemi - refractive index of the semiconductor (units)
    The reflectivity is in the range of 0 to 1.
    """
    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - nSemi) / (n1 + nSemi)
    θ = (2 * π * n1 * thickness) / wavelength
    reflectivity = 100 * (r1 * r1 + r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ)) / (
        1 + r1 * r1 * r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ))
    return reflectivity


def DLARC_refl(wavelength, n0, n1, n2, nSemi, thickness1, thickness2):
    """Return the reflectivity from a object (units) that has a double layer anti-reflection coating, DLARC.
    Where:
    n0 - refractive index of the ambient (units)
    n1 - refractive index of the dielectric layer 1 (units)
    n2 - refractive index of the dielectric layer 2 (units)
    nSemi - refractive index of the semiconductor (units)
    wavelength, thickness1, thickness 2 all in same units (m) or (nm) etc.
    """
    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - n2) / (n1 + n2)
    r3 = (n2 - nSemi) / (n2 + nSemi)
    θ1 = (2 * π * n1 * thickness1) / wavelength
    θ2 = (2 * π * n2 * thickness2) / wavelength

    numerator = r1 * r1 + r2 * r2 + r3 * r3 + r1 * r1 * r2 * r2 * r3 * r3 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))
    denominator = 1 + r1 * r1 * r2 * r2 + r1 * r1 * r3 * r3 + r3 * r3 * r2 * r2 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))

    return numerator / denominator


# ******** Section: Semiconductors *********


def probability_fermi_dirac(E, Ef, T):
    """Return the fermi dirac function (units) where E is the energy (), Ef is the fermi
    given the energies in electron volts  """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) + 1.0)


def probability_maxwell_boltzmann(E, Ef, T):
    """ Given the energies in electron volts return the fermi dirac function """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT))


def probability_bose_einstein(E, Ef, T):
    """ Given the energies in electron volts return the fermi dirac function """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) - 1.0)


def lifetime_f_length(L, D):
    """
    Return the lifetime (s)
    given the  diffusion length (cm) and diffusivity (cm2/s)
    """
    return np.sqrt(L * D)


def diffusion_length(lifetime, diffusivity):
    """Return carrier diffusion length (cm)
    given carrier lifetime(s) and diffusivity (units)
    """
    return np.sqrt(lifetime * diffusivity)


def tau_b__tau_eff(tau_eff, S, thickness):
    """Return the bulk lifetime (s)
    Given tau_eff (s)
    surface recombination (cm/s)
    thickness (cm)
    """
    return tau_eff - thickness / (2 * S)


def Vt(T=298.15):
    """Return thermal voltage (volts) at given temperature, T(Kelvin).
    The default temperature is 298.15 K, which is equal to 25 °C"""
    return k * T / q


def diffusivity(mobility, T=298.15):
    """Return the diffusivity (cm²/s) given the mobility (cm²/Vs)
    This is also known as the Einstein relation"""
    return mobility * k * T / q


def mobility(D, T=298.15):
    """Return the mobility of carriers (cm²/Vs) where D is the diffusivity (cm²/s).
    This is also known as the Einstein relation"""
    return D * q / (k * T)


def equilibrium_carrier(N, ni=8.6e9):
    """Return the majority and minority carrier concentrations (cm-3) of a semiconductor at equilibrium
    where N is the doping level (cm-3) and ni is the intrinsic carrier concentratoin (cm-3)
    Strictly N and ni just have to be in the same units but (cm-3 is almost always used."""
    majority = N
    minority = N / (ni ** 2)
    return majority, minority


def conductivity(n, p, ue, uh):
    """Return the conductivity of a material(siemens)
    Where:
    n - concentration of electrons (cm-3)
    p - concentration of holes (cm-3)
    ue - electron mobility (cm²/Vs)
    uh - hole mobility (cm²/Vs)"""
    return q * ue * n + q * uh * p


def resistivity_Si_n(Ndonor):
    """Return the resistivity of n-type silicon (ohm cm)
    given the doping of donors(cm-3)"""
    n_minority = ni_Si() ** 2 / Ndonor
    return 1 / ((q * mob_thurber(Ndonor, False) * Ndonor) + (q * mob_thurber(n_minority, False, False) * n_minority))


def resistivity_Si_p(Nacceptor):
    """Return the resistivity of p-type silicon (ohm cm)
    given the doping of acceptors(cm-3)"""
    n_minority = ni_Si() ** 2 / Nacceptor
    return 1 / ((q * mob_thurber(Nacceptor) * Nacceptor) + (q * mob_thurber(n_minority, True, False) * n_minority))


def mob_thurber(N, p_type=True, majority=True):
    """Return the mobility of carriers in silicon according to the model of Thurbur as a function of doping
    Where:
    N - doping level (cm-3)
    p_type is True or 1 for p doped material and False or 0 for n-type.
    majority is True or 1 for majority carriers and False or 0 for minority carriers.
    https://archive.org/details/relationshipbetw4006thur"""
    i = 2 * p_type + majority
    # n-type minority, n-type majority, p-type minority, p-type majority
    umax = [1417, 1417, 470, 470][i]
    umin = [160, 60, 155, 37.4][i]
    Nref = [5.6e16, 9.64E16, 1e17, 2.82E17][i]
    a = [0.647, 0.664, 0.9, 0.642][i]
    return umin + (umax - umin) / (1 + ((N / Nref) ** a))


def mob_masetti_phos(N):
    """Return the mobility of carriers (cm²/Vs) in phosphorus doped silicon accoording to the model of Masetti1983.
    Where N is the doping (cm-3)
    http://dx.doi.org/10.1109%2FT-ED.1983.21207"""
    µmax = 1414
    µmin = 68.5
    u1 = 56.1
    Cr = 9.20e16
    Cs = 3.41e20
    a = 0.711
    b = 1.98
    return µmin + (µmax - µmin) / (1 + ((N / Cr) ** a)) - u1 / (1 + ((Cs / N) ** b))


def mob_klassen(Nd, Na, Δn=1, T=298.16):
    """Return the mobility (cm2/Vs)
    given the doping etc."""
    s1 = 0.89233
    s2 = 0.41372
    s3 = 0.19778
    s4 = 0.28227
    s5 = 0.005978
    s6 = 1.80618
    s7 = 0.72169
    r1 = 0.7643
    r2 = 2.2999
    r3 = 6.5502
    r4 = 2.367
    r5 = -0.01552
    r6 = 0.6478
    fCW = 2.459
    fBH = 3.828
    mh_me = 1.258
    me_m0 = 1

    T = 298.16
    n0, p0 = equilibrium_carrier(Nd)

    n_i = 8.31E+09

    cA = 0.5
    cD = 0.21
    Nref_A = 7.20E+20
    Nref_D = 4.00E+20

    p = p0 + Δn
    n = n0 + Δn
    cc = p + n

    Za_Na = 1 + 1 / (cA + (Nref_A / Na) ** 2)
    Zd_Nd = 1 + 1 / (cD + (Nref_D / Nd) ** 2)

    Na_h = Za_Na * Na
    Nd_h = Zd_Nd * Nd

    phosphorus = [1414, 68.5, 9.20E16, 0.711, 2.285]
    boron = [470.5, 44.9, 2.23E+17, 0.719, 2.247]

    boron_µmax = 470.5
    boron_µmin = 44.9
    boron_Nref_1 = 2.23E+17
    boron_α = 0.719
    boron_θ = 2.247

    phosphorus_µmax = 1414
    phosphorus_µmin = 68.5
    phosphorus_Nref_1 = 9.20E16
    phosphorus_α = 0.711
    phosphorus_θ = 2.285

    µ_eN = phosphorus_µmax ** 2 / (phosphorus_µmax - phosphorus_µmin) * (T / 300) ** (3 * phosphorus_α - 1.5)
    µ_hN = boron_µmax ** 2 / (boron_µmax - boron_µmin) * (T / 300) ** (3 * boron_α - 1.5)

    µ_ec = phosphorus_µmax * phosphorus_µmin / (phosphorus_µmax - phosphorus_µmin) * (300 / T) ** 0.5
    µ_hc = boron_µmax * boron_µmin / (boron_µmax - boron_µmin) * (300 / T) ** 0.5

    µ_eD = µ_eN * (phosphorus_Nref_1 / Nd) ** phosphorus_α + µ_ec * (cc / Nd)
    µ_hA = µ_hN * (boron_Nref_1 / Na) ** boron_α + µ_hc * (cc / Na)

    Ne_sc = Na_h + Nd_h + p
    Nh_sc = Na_h + Nd_h + n

    PBHe = 1.36e+20 / cc * me_m0 * (T / 300) ** 2
    PBHh = 1.36e+20 / cc * mh_me * (T / 300) ** 2

    PCWe = 3.97e+13 * (1 / (Zd_Nd ** 3 * (Nd_h + Na_h + p)) * ((T / 300) ** 3)) ** (2 / 3)
    PCWh = 3.97e+13 * (1 / (Za_Na ** 3 * (Nd_h + Na_h + n)) * ((T / 300) ** 3)) ** (2 / 3)

    Pe = 1 / (fCW / PCWe + fBH / PBHe)
    Ph = 1 / (fCW / PCWh + fBH / PBHh)

    G_Pe = 1 - s1 / ((s2 + (1 / me_m0 * 300 / T) ** s4 * Pe) ** s3) + s5 / (((me_m0 * 300 / T) ** s7 * Pe) ** s6)
    G_Ph = 1 - s1 / ((s2 + (1 / (me_m0 * mh_me) * T / 300) ** s4 * Ph) ** s3) + s5 / (
        ((me_m0 * mh_me * 300 / T) ** s7 * Ph) ** s6)

    F_Pe = (r1 * Pe ** r6 + r2 + r3 / mh_me) / (Pe ** r6 + r4 + r5 / mh_me)
    F_Ph = (r1 * Ph ** r6 + r2 + r3 * mh_me) / (Ph ** r6 + r4 + r5 * mh_me)

    Ne_sc_eff = Nd_h + G_Pe * Na_h + p / F_Pe
    Nh_sc_eff = Na_h + G_Ph * Nd_h + n / F_Ph

    # Lattice Scattering
    µ_eL = phosphorus_µmax * (300 / T) ** phosphorus_θ
    µ_hL = boron_µmax * (300 / T) ** boron_θ

    µe_Dah = µ_eN * Ne_sc / Ne_sc_eff * (phosphorus_Nref_1 / Ne_sc) ** phosphorus_α + µ_ec * ((p + n) / Ne_sc_eff)
    µh_Dae = µ_hN * Nh_sc / Nh_sc_eff * (boron_Nref_1 / Nh_sc) ** boron_α + µ_hc * ((p + n) / Nh_sc_eff)

    µe = 1 / (1 / µ_eL + 1 / µe_Dah)
    µh = 1 / (1 / µ_hL + 1 / µh_Dae)

    return µe, µh


def Eg0_paessler(T=298.15):
    """Return the bandgap of silicon (eV) according to Paessler2002, where T is the temperature (K).
    Code adapted from Richter Fraunhofer ISE
    https://doi.org/10.1103/PhysRevB.66.085201
    """
    # constants from Table I on page 085201-7
    α = 3.23 * 0.0001  # (eV/K)
    Θ = 446  # (K)
    Δ = 0.51
    Eg0_T0 = 1.17  # eV     band gap of Si at 0 K

    Tdelta = 2 * T / Θ
    wurzel = (1 + π ** 2 / (3 * (1 + Δ ** 2)) * Tdelta ** 2 + (
        3 * Δ ** 2 - 1) / 4 * Tdelta ** 3 + 8 / 3 * Tdelta ** 4 + Tdelta ** 6) ** (1 / 6)
    Eg0 = Eg0_T0 - α * Θ * ((1 - 3 * Δ ** 2) / (np.exp(Θ / T) - 1) + 3 / 2 * Δ ** 2 * (wurzel - 1))
    return Eg0


def ni_Si(T=298.15):
    """Return the intrinsic carrier concentration of silicon (cm**-3) according to Sproul94, where T is the temperature (K)
    http://dx.doi.org/10.1063/1.357521"""
    return 9.38e19 * (T / 300) * (T / 300) * np.exp(-6884 / T)


def ni_misiakos(T=298.15):
    """
    Return the intrinsic carrier concentration (cm-3) without band gap narrowing according to Misiakos,
    where T is the temperature (K).
    DOI http://dx.doi.org/10.1063/1.354551
    """
    return 5.29E+19 * (T / 300) ** 2.54 * np.exp(-6726 / T)


def n_ieff(N_D, N_A, Δn, T=298.15):
    """Return effective ni (cm-3)
    given
    donor concentration N_D=n0 (1/cm³)      only one dopant type possible
    acceptor concentration N_A=p0 (1/cm³)    only one dopant type possible
    excess carrier density (1/cm³)
    temperature (K)
    calculation of the effective intrinsic concentration n_ieff including BGN
    according to Altermatt JAP 2003
    """
    # calculation of fundamental band gap according to Pässler2002
    Eg0 = Eg0_paessler(T)

    # n_i without BGN according to Misiakos93, parameterization fits very well
    # to value of Altermatt2003 at 300K
    ni0 = ni_misiakos(T)

    ni = ni0  # ni0 as starting value for n_ieff for calculation of n0 & p0

    n0 = np.where(N_D > N_A, N_D, N_A / ni ** 2)
    p0 = np.where(N_D > N_A, N_D / ni ** 2, N_A)

    # self-conistent iterative calculation of n_ieff

    for i in range(5):  # lazy programmer as it converges pretty fast anyway
        n = n0 + Δn
        p = p0 + Δn
        dEc, dEv = bandgap_schenk(n, p, N_A, N_D, Δn, T)
        ni = ni0 * np.exp(q * (dEc + dEv) / (2 * k * T))  # there is something wrong here as the units don't match up.
        n0 = np.where(N_D > N_A, N_D, N_A / ni ** 2)
        p0 = np.where(N_D > N_A, N_D / ni ** 2, N_A)

    # print('iterations',ni)
    return ni


def bandgap_schenk(n_e, n_h, N_D, N_A, Δn, T=298.15):
    """
    returns the band gap narowing in silicon
    delta conduction band, delta valence band in eV
    given:

    n_e => total electron density with Δn (1/cm³)
    n_h => total hole density with Δn (1/cm³)
    N_A => acceptor concentration (1/cm³)
    N_D => donor concentration (1/cm³)
    Δn  => excess carrier density (1/cm³)
    T   => temperature (K)

    Band-gap narrowing after Schenk 1998, JAP 84(3689))
    model descriped very well in K. McIntosh IEEE PVSC 2010
    model confirmed by Glunz2001 & Altermatt2003
    nomenclatur and formula no. according to McIntosh2010, table no. according to Schenk1998
    ==========================================================================
    Input parameters:

    ==========================================================================
    Code adapted from Richter at Fraunhofer ISE
    http://dx.doi.org/10.1063%2F1.368545
    """

    # Silicon material parameters (table 1)
    m_e_ = 0.321  # m_e/m_0 -> relative electron mass
    m_h_ = 0.346  # m_h/m_0 -> relative hole mass
    g_e = 12  # degeneracy factor for electrons
    g_h = 4  # degeneracy factor for holes
    m_ = 0.1665  # µ*/m_0 -> reduced effective mass / m_0
    alfa_e = 0.5187  # µ*/m_e
    alfa_h = 0.4813  # µ*/m_h
    Ry_ex = 0.01655  # eV    excitonic Rydberg constant
    alfa_ex = 0.0000003719  # cm     excitonic Bohr radius
    epsilon_s = 11.7  # static dielectric constant

    # Parameters for Pade-Approximation (tab. 2 & 3)
    b_e = 8
    b_h = 1
    c_e = 1.3346
    c_h = 1.2365
    d_e = 0.893
    d_h = 1.153
    p_e = 7 / 30
    p_h = 7 / 30
    h_e = 3.91
    h_h = 4.2
    j_e = 2.8585
    j_h = 2.9307
    k_e = 0.012
    k_h = 0.19
    q_e = 3 / 4
    q_h = 1 / 4

    # ==========================================================================
    # pre-calculations:
    F = (k * T / q) / Ry_ex  # eq. 29
    a3 = alfa_ex ** 3

    # Normalizing of the densities
    n_e *= a3
    n_h *= a3
    N_D *= a3
    N_A *= a3

    # for eq. 33 (normalized)
    n_sum_xc = n_e + n_h
    n_p_xc = alfa_e * n_e + alfa_h * n_h

    # for eq. 37 (normalized)
    n_sum_i = N_D + N_A  # eq.39 bzw. eq. 29
    n_p_i = alfa_e * N_D + alfa_h * N_A  # eq.39 bzw. eq. 29

    Ui = n_sum_i ** 2 / F ** 2  # eq. 38
    n_ionic = n_sum_i  # McIntosh2010

    # exchange quasi-partical shift Eq33:
    delta_xc_h = -(
        (4 * π) ** 3 * n_sum_xc ** 2 * ((48 * n_h / (π * g_h)) ** (1 / 3) + c_h * np.log(1 + d_h * n_p_xc ** p_h)) + (
            8 * π * alfa_h / g_h) * n_h * F ** 2 + np.sqrt(8 * π * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * π) ** 3 * n_sum_xc ** 2 + F ** 3 + b_h * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)
    delta_xc_e = -(
        (4 * π) ** 3 * n_sum_xc ** 2 * ((48 * n_e / (π * g_e)) ** (1 / 3) + c_e * np.log(1 + d_e * n_p_xc ** p_e)) + (
            8 * π * alfa_e / g_e) * n_e * F ** 2 + np.sqrt(8 * π * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * π) ** 3 * n_sum_xc ** 2 + F ** 3 + b_e * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)

    # ionic quasi-partical shift Eq37:
    delta_i_h = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / π) * (1 + h_h * np.log(1 + np.sqrt(n_sum_i) / F)) + j_h * Ui * n_p_i ** 0.75 * (
            1 + k_h * n_p_i ** q_h))
    delta_i_e = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / π) * (1 + h_e * np.log(1 + np.sqrt(n_sum_i) / F)) + j_e * Ui * n_p_i ** 0.75 * (
            1 + k_e * n_p_i ** q_e))

    # rescale BGN
    dE_gap_h = -Ry_ex * (delta_xc_h + delta_i_h)
    dE_gap_e = -Ry_ex * (delta_xc_e + delta_i_e)
    return dE_gap_e, dE_gap_h


# ******** SubSection: Bulk Recombination *********
def U_radiative(n, p):
    B_rad = 4.73e-15
    U_radiative = n * p * B_rad
    return U_radiative


def U_radiative_alt(n0, p0, Δn, T=298.15):
    n_p = n0 + p0 + 2 * Δn
    n = n0 + Δn
    p = p0 + Δn
    B_low = 4.73e-15
    b_min = 0.2 + (0 - 0.2) / (1 + (T / 320) ** 2.5)
    b1 = 1.5E+18 + (10000000 - 1.5E+18) / (1 + (T / 550) ** 3)
    b3 = 4E+18 + (1000000000 - 4E+18) / (1 + (T / 365) ** 3.54)
    B_rel = b_min + (1 - b_min) / (1 + (0.5 * n_p / b1) ** 0.54 + (0.5 * n_p / b3) ** 1.25)
    B_rad = B_low * B_rel
    U_radiative_alt = n * p * B_rad
    return U_radiative_alt


def U_SRH(n, p, Et, τ_n, τ_p, ni_eff=8.5e9, T=298.15):
    """Return the shockley read hall recombination cm-3
    given Et (eV) trap level from intrinsic"""
    n1 = ni_eff * np.exp(q * Et / k / T)
    p1 = ni_eff * np.exp(-q * Et / k / T)
    U_SRH = (n * p - ni_eff ** 2) / (τ_p * (n + n1) + τ_n * (p + p1))
    return U_SRH


def U_auger_richter(n0, p0, Δn, ni_eff):
    """Return the auger recombination
    18 and 19
    https://doi.org/10.1016/j.egypro.2012.07.034"""
    B_n0 = 2.5E-31
    C_n0 = 13
    D_n0 = 3.3E+17
    exp_n0 = 0.66
    B_p0 = 8.5E-32
    C_p0 = 7.5
    D_p0 = 7E+17
    exp_p0 = 0.63
    C_dn = 3E-29
    D_dn = 0.92
    g_eeh = (1 + C_n0 * (1 - np.tanh((n0 / D_n0) ** exp_n0)))
    g_ehh = (1 + C_p0 * (1 - np.tanh((p0 / D_p0) ** exp_p0)))
    np_ni2 = (n0 + Δn) * (p0 + Δn) - ni_eff ** 2
    U = np_ni2 * (B_n0 * n0 * g_eeh + B_p0 * p0 * g_ehh + C_dn * Δn ** D_dn)
    return U


def U_low_doping(n0, p0, Δn):
    """recombination due to Auger and radiative
    equation 21 in DOI: 10.1103/PhysRevB.86.165202"""
    B_low = 4.73e-15
    n = n0 + Δn
    p = p0 + Δn
    U = Δn / (n * p * (8.7e-29 * n0 ** 0.91 + 6.0e-30 * p0 ** 0.94 + 3.0e-29 * Δn ** 0.92 + B_low))
    return U


def lifetime(U, Δn):
    """Return the lifetime (seconds) where U is the recombination  and Δn is the excess minority carrier density.
        This is the definition of lifetime"""
    return Δn / U


# not sure if I should keep these
def lifetime_auger(Δn, Ca=1.66e-30):
    """Returns the Auger lifetime (s) at high level injection
    given the injection level (cm-3)"""
    return 1 / (Ca * Δn ** 2)


def lifetime_SRH(N, Nt, Et, σ_n, σ_p, Δn, T=298.15):
    Nv = 31000000000000000000 * (T / 300) ** 1.85
    Nc = 28600000000000000000 * (T / 300) ** 1.58
    Eg = 1.1246
    vth = 11000000 * (T / 300) ** 0.5
    p0 = N
    n0 = (ni_Si(T) ** 2) / N
    τ_n0 = 1 / (Nt * σ_n * vth)
    τ_p0 = 1 / (Nt * σ_p * vth)
    n1 = Nc * np.exp(-Et / Vt())
    p1 = Nv * np.exp((-Et - Eg) / Vt())
    k_ratio = σ_n / σ_p
    τ_SRH = (τ_p0 * (n0 + n1 + Δn) + τ_n0 * (p0 + p1 + Δn)) / (n0 + p0 + Δn)
    return τ_SRH


# surface recombination
def U_surface(n, p, Sn, Sp, n1=8.3e9, p1=8.3e9, ni=8.3e9):
    """Return the carrier recombination (/s) at a surface.
    Where.
    Sn, Sp: surface recombination for electrons and holes
    n1, p1 XXX
    ni - intrinsice carrier concentratoin (cm-3)"""
    U_surface = Sn * Sp * (n * p - ni ** 2) / (Sn * (n + n1) + Sp * (p + p1))
    return U_surface


# ******** Section: Solar Cells *********


def IQE_emitter(ab, We, Le, De, Se):
    """Return the internal quantum efficiency of a solar cell emitter
    Where:
    ab - absorption coefficient (/cm)
    We - thickness of the emitter (cm)
    De - diffusivty of carriers in the emitter (cm²/s)
    Se - recombination at the front surface (cm/s)
    Hovel, I think."""
    GF = ((Se * Le / De) + ab * Le - np.exp(-ab * We) * ((Se * Le / De) * np.cosh(We / Le) + np.sinh(We / Le))) / (
        (Se * Le / De) * np.sinh(We / Le) + np.cosh(We / Le)) - Le * ab * np.exp(-ab * We)
    QEE = (Le * ab / (ab * ab * Le * Le - 1)) * GF
    return QEE


def IQE_base(ab, We_Wd, Wb, Lb, Db, Sb):
    """Return quantum efficiency of the base of a solar cell
    where:
    ab -  absorption coefficient (cm)
    We_Wd - junction depth (cm)
    Sb - surface recombination velocity (cm/s)
    Lb - diffusion length of minority carrier in the base (cm)
    Db - diffusivity of minority carriers in the base (cm²/Vs)
    """
    GF = (ab * Lb - (
        (Sb * Lb / Db) * (np.cosh(Wb / Lb) - np.exp(-ab * Wb)) + np.sinh(Wb / Lb) + Lb * ab * np.exp(-ab * Wb)) / (
              (Sb * Lb / Db) * np.sinh(Wb / Lb) + np.cosh(Wb / Lb)))
    QEB = (np.exp(-ab * We_Wd)) * (Lb * ab / (ab ** 2 * Lb ** 2 - 1)) * GF
    return QEB


def IQE_depletion(ab, We, Wd):
    QED = np.exp(-ab * We) * (1 - np.exp(-ab * Wd))
    return QED


def IQE(ab, Wd, Se, Le, De, We, Sb, Wb, Lb, Db):
    """ We is the thickness of emitter and start of the junction
    """
    QEE = IQE_emitter(ab, We, Le, De, Se)
    QEB = IQE_base(ab, We + Wd, Wb, Lb, Db, Sb)
    QED = IQE_depletion(ab, We, Wd)
    IQEt = QEE + QEB + QED
    return QEE, QEB, QED, IQEt


def QE2SR(wavelength, QE):
    """'converts a QE in units to spectral response
    given the wavelength (nm)"""
    spectral_response = QE * wavelength / 1239.8
    return spectral_response


def SR2QE(wavelength, spectral_response):
    """convert SR (A/W) to QE (unit 0 to 1)
    assumes that the wavelegth is in  nm"""
    QE = spectral_response * wavelength / 1239.8
    return QE


def impliedV(Δn, N, T=298.15):
    """Return voltage (V) where Δn is the excess carrier concentration (cm-3), N is the doping (cm-3) and
    T is the temperature (K). Implied voltage is often used to convert the carrier concentration in a lifetime
    tester to voltage."""
    return Vt(T) * np.log((Δn + N) * Δn / ni_Si(T) ** 2)


def implied_carrier(V, N, ni=8.6e9, T=298.15):
    """Return excess carrier concentration (cm-3)
    Given voltage and doping determine """
    Δn = (-N + np.sqrt(N ** 2 + 4 * ni ** 2 * np.exp(V / Vt(T)))) / 2
    return Δn


def J0_layer(W, N, D, L, S, ni=8.6e9):
    """Return the saturation current density (A/cm2) for the narrow case.
    Whete:
    W - layer thickness (cm)
    N - doping (cm-3)
    L - diffusion length (cm)
    S - surface recombination velocity (cm/s)
    Optional:
    ni - intrinsic carrier concentration (cm-3)
    """
    F = (S * np.cosh(W / L) + D / L * np.sinh(W * L)) / (D / L * np.cosh(W * L) + S * np.sinh(W / L))
    return q * ni ** 2 * F * D / (L * N)


# def J0(ni, We, Ne, De, Le, Se, Nb, Wb, Db, Lb, Sb):
#    '''determines J0, the dark saturation current, under the narrow base diode
# condition where L > W.'''
#    Fe = (Se*np.cosh(We/Le)+De/Le*np.sinh(We*Le))/(De/Le*np.cosh(We*Le)+Se*np.sinh(We/Le))
#    Fb = (Sb*np.cosh(Wb/Lb)+Db/Lb*np.sinh(Wb*Lb))/(Db/Lb*np.cosh(Wb*Lb)+Sb*np.sinh(Wb/Lb))
#    J0 = q*ni**2*(Fe*De/(Le*Ne)+ Fb*Db/(Lb*Nb))
#    return J0

def efficiency(Voc, Isc, FF, A=1):
    """Return the efficiency of a solar cell (units not percentage)given Voc (volts), Isc in (amps) and  FF (units).
    also works for Jsc since area of 1 is assumed
    """
    return 1000 * Voc * Isc * FF / A


def current2gen(I):
    """Return generation (eh pairs/s)
    given current (amps)
    """
    return I / q


def I_diode(V, I0, T=298.15):
    """Return the current (A) in an ideal diode where I0 is the saturation current (A),
    V is the voltage across the junction (volts), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return I0 * np.exp(V / Vt(T) - 1)


def I_cell(V, IL, I0, T=298.15):
    """Return current (amps) of a solar cell
    given voltage, light generated current, I0
    also works for J0
    """
    return IL - I0 * np.exp(V / Vt(T))


def I_cell_Rshunt(V, IL, I0, Rshunt, T=298.15):
    """Return current (A) of a solar cell from   """
    return IL - I0 * np.exp(V / Vt(T)) - V / Rshunt


def V_Rseries(voltage, I, Rs):
    """Returns the voltage of a solar cells under the effect of series resistance"""
    return voltage - I * Rs


def Voc(IL, I0, n=1, T=298.15):
    """Return the open circuit voltage, Voc, (volts) from IL(A) and I0(A).
    IL and Io must be in the same units, Eg, (A), (mA) etc
    Using (mA/cm**2) uses J0 and JL instead.
    """
    return n * Vt(T) * np.log(IL / I0 + 1)


def V_cell(I, IL, I0, T=298.15):
    """Return the voltage (V) in an ideal solar cell where I0 is the saturation current (A),
    I is the current (A), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return Vt(T) * np.log((IL - I) / I0 + 1)


def cell_params(V, I):
    """Return key parameters of a solar cell IV curve where V is a voltage array and
    I is a current array, both with type numpy.array.
    Voc (V), Isc (A), FF, Vmp(V), Imp(A) given voltage vector in (volts)
    current vector in (amps) or (A/cm²)
    If I is in (A/cm²) then Isc will be Jsc and Imp will be Jmp. No attempt is made to fit the fill factor.
    """
    Voc = np.interp(0, -I, V)
    Isc = np.interp(0, V, I)
    idx = np.argmax(V * I)
    Vmp = V[idx]
    Imp = I[idx]
    FF = Vmp * Imp / (Voc * Isc)
    return Voc, Isc, FF, Vmp, Imp


def finger_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp):
    """Return the fractional resistivity power loss in a finger (0 to 1)
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    return (L ** 2 * Jmp * Sf * resistivity) / (3 * wf * df * Vmp) * 100.0


def finger_shading(wf, Sf):
    """Return the fractional power loss due to finger shading (0 to 1) where wf is the wideth of the finger and Sf is the finger spacing."""
    return (wf / Sf) * 100.0


def finger_sheet(Sf, Jmp, Rsheet, Vmp):
    return (Sf ** 2 * Jmp * Rsheet) / (12 * Vmp) * 100.0


def finger_total_loss(L, Jmp, Sf, resistivity, Rsheet, wf, df, Vmp):
    """Return the fractional power loss in a finger
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    Presistivity = finger_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp)
    Pshading = finger_shading(wf, Sf)
    Psheet = finger_sheet(Sf, Jmp, Rsheet, Vmp)
    return Presistivity + Pshading + Psheet, Presistivity, Pshading, Psheet


def FF(Vmp, Imp, Voc, Isc):
    """Return FFv the fill factor of a solar cell.
    given Voc - open circuit voltage (volts)"""
    return (Vmp * Imp) / (Voc * Isc)


def FF_ideal(Voc, ideality=1, T=298.15):
    """Return the FF (units)
    given Voc - open circuit voltage (volts), ideality factor, defaults to 1 (units) """
    voc = normalised_Voc(Voc, ideality, T)
    FF0 = (voc - np.log(voc + 0.72)) / (voc + 1)
    return FF0


def normalised_Voc(Voc, ideality, T=298.15):
    """Return the normalised voc of a solar cell where Voc is the open-circuit voltage, 'ideality' is the ideality factor
     and T is the temperature (K)"""
    return Voc / (ideality * Vt(T))


def FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
        Isc - short circuit current (amps)
        Rseries - series resistance (ohms)
        ideality factor (units)
        T - temperature (K)
    """
    # voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rs = Rseries / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - 1.1 * rs) + (rs ** 2 / 5.4)
    return FF


def FF_Rsh(Voc, Isc, Rshunt, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
    """
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - ((voc + 0.7) * FF0) / (voc * rsh))
    return FF


def FF_RsRsh(Voc, Isc, Rseries, Rshunt, ideality=1, T=298.15):
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FFRs = FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15)
    FFRsRsh = FFRs * (1 - ((voc + 0.7) * FFRs) / (voc * rsh))
    return FFRsRsh


# ******** Section: silicon material properties *********
# silicon material properties

def optical_properties(fname=None):
    """Returns an array with the optical properties of a material
    column 0 - wavelngth (nm)
    column 1 - absorption coefficient (/cm)
    column 2 - real refractive index
    column 3 - imaginary refractive index
    if so file is given then silicon is used
    Eg: wavelength, abs_coeff, n, k = optical_properties()
    """
    if fname is None:
        package_path = os.path.dirname(os.path.abspath(__file__))
        fname = os.path.join(package_path, 'silicon_optical_properties.txt')
    wavelength, abs_coeff, nd, kd = np.loadtxt(fname, skiprows=1, unpack=True)
    return wavelength, abs_coeff, nd, kd


# processing
def phos_active(T):
    """Return the active limit of phosphorous in silicon
    given temperature (K)"""
    return 1.3e22 * np.exp(-0.37 * eV / (k * T))


def phos_solubility(T):
    """Return the solubility limit of phosphorous in silicon
     given the temperature (K)"""
    return 2.45e23 * np.exp(-0.62 * eV / (k * T))

# modules(Pedro)
def read_cell_info(selected):
    package_path = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(package_path, 'cell_info.txt')

    with open(fname,"r") as f:
        for line in f:
            col1, col2, col3, col4 = line.split()
            if (col1==selected):
                semicondutor=col1
                J_SC=float(col2)
                V_OC=float(col3)
                J_0=float(col4)
    return semicondutor,J_SC,V_OC,J_0
def module_current(M,N,T,material):
    semicondutor,J_SC,V_OC,J_0 = read_cell_info(material)
    I_0=J_0*15.6*15.6
    I_L=J_SC*15.6*15.6
    V_T=V_OC*N
    I_total = M*(I_L-I_0*np.exp((q*V_T/N)/(k_eV*T)-1))
    return I_total
