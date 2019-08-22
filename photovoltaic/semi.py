""""Basic semiconductor equations. The cell module contains code related to solar cells and devices.
"""

from .core import k, q
import numpy as np

def lifetime(U, Δn):
    """Return the lifetime (seconds) where U is the recombination  and Δn is the excess minority carrier density.
        This is the definition of lifetime"""
    return Δn / U


def probability_fermi_dirac(E, Ef, T=298.15):
    """Return the fermi dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """

    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) + 1.0)


def probability_maxwell_boltzmann(E, Ef, T=298.15):
    """Return the Maxwell Bolzman dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT))


def probability_bose_einstein(E, Ef, T=298.15):
    """Return the Maxwell Bolzman dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) - 1.0)


def length2lifetime(length, diffusivity):
    """
    Convert the carrier diffusion length, L (cm), to minority lifetime (s), given the D, the diffusivity (cm2/s).
    """
    return np.sqrt(length * diffusivity)


def lengthflifetime(lifetime, diffusivity):
    """Convert the carrier diffusion length, L (cm), to minority lifetime (s), given the D, the diffusivity (cm2/s).
    """
    return np.sqrt(lifetime * diffusivity)


def bulkfeffective(tau_eff, S, thickness):
    """Return the bulk lifetime (s) where taus
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
    """Return the mobility of carriers (cm²/Vs) given D is the diffusivity (cm²/s).
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
	
def U_SRH(n, p, Et, τ_n, τ_p, ni_eff=8.5e9, T=298.15):
    """Return the shockley read hall recombination cm-3
    given Et (eV) trap level from intrinsic"""
    n1 = ni_eff * np.exp(q * Et / k / T)
    p1 = ni_eff * np.exp(-q * Et / k / T)
    U_SRH = (n * p - ni_eff ** 2) / (τ_p * (n + n1) + τ_n * (p + p1))
    return U_SRH
