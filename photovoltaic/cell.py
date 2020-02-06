from .core import q, pi, Vt
from .si import ni_sproul
import numpy as np


def FF(Vmp, Imp, Voc, Isc):
    """Return FFv the fill factor of a solar cell.
    given Voc - open circuit voltage (volts)"""
    return (Vmp * Imp) / (Voc * Isc)


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


def FF_RsRsh(Voc, Isc, Rseries, Rshunt, ideality=1, T=298.15):
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FFRs = FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15)
    FFRsRsh = FFRs * (1 - ((voc + 0.7) * FFRs) / (voc * rsh))
    return FFRsRsh


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


def FF_ideal(Voc, ideality=1, T=298.15):
    """Return the FF (units)
    given Voc - open circuit voltage (volts), ideality factor, defaults to 1 (units) """
    voc = normalised_Voc(Voc, ideality, T)
    FF0 = (voc - np.log(voc + 0.72)) / (voc + 1)
    return FF0


def iqe(ab, Wd, Se, Le, De, We, Sb, Wb, Lb, Db):
    """ We is the thickness of emitter and start of the junction
    """
    QEE = iqe_emitter(ab, We, Le, De, Se)
    QEB = iqe_base(ab, We + Wd, Wb, Lb, Db, Sb)
    QED = iqe_depletion(ab, We, Wd)
    IQEt = QEE + QEB + QED
    return QEE, QEB, QED, IQEt


def I_cell(V, IL, I0, T=298.15):
    """Return current (amps) of a solar cell
    given voltage, light generated current, I0
    also works for J0
    """
    return IL - I0 * np.exp(V / Vt(T))


def I_cell_Rshunt(V, IL, I0, Rshunt, T=298.15):
    """Return current (A) of a solar cell from   """
    return IL - I0 * np.exp(V / Vt(T)) - V / Rshunt


def I_diode(V, I0, T=298.15):
    """Return the current (A) in an ideal diode where I0 is the saturation current (A),
    V is the voltage across the junction (volts), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return I0 * np.exp(V / Vt(T) - 1)


def J0_layer(W, N, D, L, S, ni=8.6e9):
    """Return the saturation current density (A/cm2) for the narrow case.
    Where:
    W - layer thickness (cm)
    N - doping (cm-3)
    L - diffusion length (cm)
    S - surface recombination velocity (cm/s)
    Optional:
    ni - intrinsic carrier concentration (cm-3)
    """
    F = (S * np.cosh(W / L) + D / L * np.sinh(W / L)) / (D / L * np.cosh(W / L) + S * np.sinh(W / L))
    return q * ni ** 2 * F * D / (L * N)


def V_Rseries(voltage, I, Rs):
    """Returns the voltage of a solar cells under the effect of series resistance"""
    return voltage - I * Rs


def V_cell(I, IL, I0, T=298.15):
    """Return the voltage (V) in an ideal solar cell where I0 is the saturation current (A),
    I is the current (A), T is the temperature (K) and n is the ideallity factor (units).
    For current density. I0 is in A/cm² and current density is returned"""
    return Vt(T) * np.log((IL - I) / I0 + 1)


def Voc(IL, I0, n=1, T=298.15):
    """Return the open circuit voltage, Voc, (volts) from IL(A) and I0(A).
    IL and Io must be in the same units, Eg, (A), (mA) etc
    Using (mA/cm²) uses J0 and JL instead.
    """
    return n * Vt(T) * np.log(IL / I0 + 1)


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


def collection(x, thickness, S, L, D):
    """Returns the collection probability (unit 0 to 1) at a distance x (cm) from the junction.
    Where: thickness is the layer thickness (cm),
    S is the surface recombination velocity (cm/s),
    L is the minority carrier diffusion length (cm) and
    D is the diffusivity (cm²/Vs)
    """
    hypSin_xL = np.sinh(-x / L)
    hypCos_xL = np.cosh(-x / L)
    hypSin_WL = np.sinh(-thickness / L)
    hypCos_WL = np.cosh(-thickness / L)
    num = S * L / D * hypCos_WL + hypSin_WL
    denom = S * L / D * hypSin_WL + hypCos_WL
    return (hypCos_xL - num / denom * hypSin_xL)


def genfcurrent(current):
    """Return generation (eh pairs/s) given current (amps)
    """
    return current / q


def efficiency(Voc, Isc, FF, A=1):
    """Return the efficiency of a solar cell (units not percentage)given Voc (volts), Isc in (amps) and  FF (units).
    also works for Jsc since area of 1 is assumed
    """
    return 1000 * Voc * Isc * FF / A


def finger_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp):
    """Return the fractional resistivity power loss in a finger (0 to 1)
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    return (L ** 2 * Jmp * Sf * resistivity) / (3 * wf * df * Vmp) * 100.0


def finger_shading(wf, Sf):
    """Return the fractional power loss due to finger shading (0 to 1)
    where wf is the width of the finger (cm) and Sf is the finger spacing (cm)."""
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


def impliedV(Δn, N, T=298.15):
    """Return voltage (V) where Δn is the excess carrier concentration (cm-3), N is the doping (cm-3) and
    T is the temperature (K). Implied voltage is often used to convert the carrier concentration in a lifetime
    tester to voltage."""
    return Vt(T) * np.log((Δn + N) * Δn / ni_sproul(T) ** 2)


def implied_carrier(V, N, ni=8.6e9, T=298.15):
    """Return excess carrier concentration (cm-3)
    Given voltage and doping determine """
    Δn = (-N + np.sqrt(N ** 2 + 4 * ni ** 2 * np.exp(V / Vt(T)))) / 2
    return Δn


def iqe_base(ab, We_Wd, Wb, Lb, Db, Sb):
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


def iqe_depletion(ab, We, Wd):
    QED = np.exp(-ab * We) * (1 - np.exp(-ab * Wd))
    return QED


def iqe_emitter(ab, We, Le, De, Se):
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


def normalised_Voc(Voc, ideality, T=298.15):
    """Return the normalised voc of a solar cell where Voc is the open-circuit voltage, 'ideality' is the ideality factor
     and T is the temperature (K)"""
    return Voc / (ideality * Vt(T))


def qefsr(spectral_response, wavelength):
    """Return the QE (unit 0 to 1) given the spectral response (A/W) and the wavelength (nm) """
    QE = spectral_response * wavelength / 1239.8
    return QE


def spectral_mismatch(spectrum, lamp, ref_SR, device_SR):
    '''returns the spectral mismatch. Where spectrum is the reference spectrum in W/m2 or similar,
    lamp is the spectrum of the lamp in W/m2 or similar,
    ref_SR is the spectral response of the reference cell used for calibration in A/W and
    device_SR is the spectral response of the device under test in A/W. '''
    M = (np.sum(lamp * device_SR) / np.sum(lamp * ref_SR)) * (np.sum(spectrum * ref_SR) / np.sum(spectrum * device_SR))
    return M


def srfqe(qe, wavelength):
    """'Return the spectral response (A/W) given the quantum efficiency (0 yo 1 units) and wavelength (nm)"""
    spectral_response = qe * wavelength / 1239.8
    return spectral_response
