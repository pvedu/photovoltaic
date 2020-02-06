""""Basic semiconductor equations. The cell module contains code related to solar cells and devices.
"""

from .core import k, q
import numpy as np

def Ei(Eg, Nc, Nv, T=298.15):
    """Return Ei (eV) the intrinsic Fermi level
    Given:
        Nc (cm-3) density of states in the conduction band
        Nv (cm-3) density of states in the valence band
        Eg (eV) band gap"""
    return Eg/2 + 0.5*k*T*np.log(Nv/Nc)


def U_SRH(n, p, τ_n0, τ_p0, Et, Ei = 0.561, ni_eff=8.5e9, T=298.15):
    """Return the Shockley Read Hall recombination cm-3
    Given:
        n (cm-3) is the total concentration of electrons,
        p (cm-3) is the total concentration of holes,
        Et (eV) is energy level of the trap from VALENCE band.
        Ei (eV) energy level of the intrinsic level,
        τ_n0 (s) lifetime of electrons,
        τ_p0 (s) lifetime of holes,
        ni_eff the intrinsic carrier concentration,
        T (K)) temperature.
    """
    n1 = ni_eff * np.exp(q * (Et - Ei) / (k*T))
    p1 = ni_eff * np.exp(-q * (Et - Ei) / (k * T))
    U_SRH = (n * p - ni_eff ** 2) / (τ_p0 * (n + n1) + τ_n0 * (p + p1))
    return U_SRH


def Vt(T=298.15):
    """Return thermal voltage (volts) at given temperature, T(Kelvin).
    The default temperature is 298.15 K, which is equal to 25 °C"""
    return k * T / q


def bulkfeffective(tau_eff, S, thickness):
    """Return the bulk lifetime (s) where taus
    Given:
        tau_eff (s)
        surface recombination (cm/s)
        thickness (cm)
    """
    return tau_eff - thickness / (2 * S)


def conductivity(n, p, ue, uh):
    """Return the conductivity of a material(siemens)
    Given:
        n (cm-3) concentration of electrons
        p (cm-3) concentration of holes
        ue (cm²/Vs) electron mobility
        uh (cm²/Vs) hole mobility """
    return q * ue * n + q * uh * p

def diffusivity(mobility, T=298.15):
    """Return the diffusivity (cm²/s) given the mobility (cm²/Vs)
    This is also known as the Einstein relation"""
    return mobility * k * T / q


def equilibrium_carrier(N, ni=8.6e9):
    """Return the majority and minority carrier concentrations (cm-3) of a semiconductor at equilibrium
    where N is the doping level (cm-3) and ni is the intrinsic carrier concentratoin (cm-3)
    Strictly N and ni just have to be in the same units but (cm-3) is almost always used."""
    majority = N
    minority = N / (ni ** 2)
    return majority, minority


def lengthflifetime(lifetime, diffusivity):
    """Return carrier diffusion length, L (cm)
    Given:
        minority lifetime (s),
        diffusivity (cm2/s).
    """
    return np.sqrt(lifetime * diffusivity)


def lifetime(U, Δn):
    """Return the lifetime (seconds) where U is the recombination  and Δn is the excess minority carrier density.
        This is the definition of lifetime"""
    return Δn / U


def lifetime0(Nt, sigma, vth=2e7):
    """Return the SRH fundemetal lifetime (s) for a trap level
    Given: Nt (/cm3) trap density
    signma - capture cross section
    vth (cm/s) thermal velocity of carriers )"""
    return 1 / (Nt * sigma * vth)


def lifetimeflength(length, diffusivity):
    """
    Convert the carrier diffusion length, L (cm), to minority lifetime (s), given the D, the diffusivity (cm2/s).
    """
    return np.sqrt(length * diffusivity)


def mobility(D, T=298.15):
    """Return the mobility of carriers (cm²/Vs) given D is the diffusivity (cm²/s).
    This is also known as the Einstein relation"""
    return D * q / (k * T)


def probability_bose_einstein(E, Ef, T=298.15):
    """Return the Maxwell Bolzman dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) - 1.0)


def probability_fermi_dirac(E, Ef, T=298.15):
    """Return the fermi dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """

    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) + 1.0)


def probability_maxwell_boltzmann(E, Ef, T=298.15):
    """Return the Maxwell Bolzman dirac function (units) where E is the energy (eV), Ef is the fermi energy (eV) """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT))


