import photovoltaic as pv
import numpy as np

print('Thermal voltage 25 degC (V):', pv.Vt())  # default is 25degC
print('Thermal voltage 300 K (V):', pv.Vt(300))
print('Silicon ni at 25 degC {:.3e} cm-3'.format(pv.ni_misiakos()))
print('Silicon ni at 300 K {:.3e} cm-3'.format(pv.ni_misiakos(300)))

print('n-type cell with doping level of 2e15')
ni = pv.ni_misiakos()
n0, p0 = pv.equilibrium_carrier(2e15, ni)  # n-type dioping at 1e15
print('Majority n {:.3e}, Minority p: {:.3e}'.format(n0, p0))
dn = 1e15
dEc, dEv = pv.bandgap_schenk(n0 + dn, p0 + dn, n0, p0, dn)
print('BGN Ec: {:.2e} eV, Ev {:.2e} eV'.format(dEc, dEv))

print('nieff {:.3e}'.format(pv.n_ieff(n0, p0, 1e15)))
B = pv.U_radiative_alt(n0, p0, dn)
print('radiative: {:.3}'.format(B))
print('Mobility of electrons as minority carriers: ', pv.u_Si_e_min(1e15))
