import numpy as np

from .core import pi, Vt, k, q, eV
import os

def resistivity_Si_n(Ndonor):
    """Return the resistivity of n-type silicon (ohm cm)
    given the doping of donors(cm-3)"""
    n_minority = ni_sproul() ** 2 / Ndonor
    return 1 / ((q * mob_thurber(Ndonor, False) * Ndonor) + (q * mob_thurber(n_minority, False, False) * n_minority))


def resistivity_Si_p(Nacceptor):
    """Return the resistivity of p-type silicon (ohm cm)
    given the doping of acceptors(cm-3)"""
    n_minority = ni_sproul() ** 2 / Nacceptor
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
    n0 = Nd
    p0 = Nd / (1e10 ** 2)


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
    wurzel = (1 + pi ** 2 / (3 * (1 + Δ ** 2)) * Tdelta ** 2 + (
        3 * Δ ** 2 - 1) / 4 * Tdelta ** 3 + 8 / 3 * Tdelta ** 4 + Tdelta ** 6) ** (1 / 6)
    Eg0 = Eg0_T0 - α * Θ * ((1 - 3 * Δ ** 2) / (np.exp(Θ / T) - 1) + 3 / 2 * Δ ** 2 * (wurzel - 1))
    return Eg0


def ni_sproul(T=298.15):
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
        (4 * pi) ** 3 * n_sum_xc ** 2 * ((48 * n_h / (pi * g_h)) ** (1 / 3) + c_h * np.log(1 + d_h * n_p_xc ** p_h)) + (
            8 * pi * alfa_h / g_h) * n_h * F ** 2 + np.sqrt(8 * pi * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * pi) ** 3 * n_sum_xc ** 2 + F ** 3 + b_h * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)
    delta_xc_e = -(
        (4 * pi) ** 3 * n_sum_xc ** 2 * ((48 * n_e / (pi * g_e)) ** (1 / 3) + c_e * np.log(1 + d_e * n_p_xc ** p_e)) + (
            8 * pi * alfa_e / g_e) * n_e * F ** 2 + np.sqrt(8 * pi * n_sum_xc) * F ** (5 / 2)) / (
                     (4 * pi) ** 3 * n_sum_xc ** 2 + F ** 3 + b_e * np.sqrt(n_sum_xc) * F ** 2 + 40 * n_sum_xc ** (
                         3 / 2) * F)

    # ionic quasi-partical shift Eq37:
    delta_i_h = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / pi) * (1 + h_h * np.log(1 + np.sqrt(n_sum_i) / F)) + j_h * Ui * n_p_i ** 0.75 * (
            1 + k_h * n_p_i ** q_h))
    delta_i_e = -n_ionic * (1 + Ui) / (
        np.sqrt(0.5 * F * n_sum_i / pi) * (1 + h_e * np.log(1 + np.sqrt(n_sum_i) / F)) + j_e * Ui * n_p_i ** 0.75 * (
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
    n0 = (ni_sproul(T) ** 2) / N
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

# processing
def phos_active(T):
    """Return the active limit of phosphorous in silicon
    given temperature (K)"""
    return 1.3e22 * np.exp(-0.37 * eV / (k * T))


def phos_solubility(T):
    """Return the solubility limit of phosphorous in silicon
     given the temperature (K)"""
    return 2.45e23 * np.exp(-0.62 * eV / (k * T))

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