from .core import arcsind, sind, joulefnm, pi, π, k, q
import numpy as np

def photon_flux(power, wavelength):
    """Return the photon flux (/s) given the power of light (watts) and wavelength (nm)
    If power is in W/m2 then flux is in m-2s-1"""
    return power / joulefnm(wavelength)


def absorptionfextinction(kd, wavelength):
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


