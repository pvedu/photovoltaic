import numpy as np

# define constants
q = 1.60217662e-19  # (coulombs) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
k_eV = 8.6173303e-05  # (eV K^-1)
Wien = 2.898e-3  # (m K)
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
pi = np.pi  #
h = 6.62607004e-34  # (J.s)
hbar = 6.62607004e-34 / (2 * pi)  # usable
c = 299792458.0  # (m s^-1)
hc_q = h * c / q  # 1.2398419745831506e-06




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


def Vt(T=298.15):
    """Return thermal voltage (volts) at given temperature, T(Kelvin).
    The default temperature is 298.15 K, which is equal to 25 °C"""
    return k * T / q

def eVfnm(x):
    """ Given wavelength (nm) of a photon return the energy (eV) """
    return hc_q * 1e9 / x


def nmfeV(x):
    """ Given energy (eV) of a photon return the wavelength (nm) """
    return hc_q * 1e9 / x


def joulefnm(x):
    """ Given wavelength (nm) of a photon return the energy (eV) """
    return h * c * 1e9 / x

