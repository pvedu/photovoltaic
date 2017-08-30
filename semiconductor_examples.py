import photovoltaic as pv
import numpy as np
import matplotlib.pyplot as plt

# from scipy.constants import *

print('Thermal voltage 25 degC (V):', pv.Vt())  # default is 25degC
print('Thermal voltage 300 K (V):', pv.Vt(300))
print('Silicon ni at 25 degC {:.3e} cm-3'.format(pv.ni0_Misiakos()))
print('Silicon ni at 300 K {:.3e} cm-3'.format(pv.ni0_Misiakos(300)))

print('n-type cell with doping level of 2e15')
ni = pv.ni0_Misiakos()
n0, p0 = pv.carrier_conc(2e15, ni)  # n-type dioping at 1e15
print('Majority n {:.3e}, Minority p: {:.3e}'.format(n0, p0))
dn = 1e15
dEc, dEv = pv.BGN_Schenk(n0 + dn, p0 + dn, n0, p0, dn)
print('BGN Ec: {:.2e} eV, Ev {:.2e} eV'.format(dEc, dEv))

print('nieff {:.3e}'.format(pv.n_ieff(n0, p0, 1e15)))
B = pv.B_altermatt(n0, p0, dn)
print('radiative: {:.3}'.format(B))


print('Mobility of electrons as minority carriers: ', pv.u_Si_e_min(1e15))  # should I use cm-3 or m?

N = np.logspace(12, 21, 200)  # sweep the doping in the base (cm-3)

# plot the mobilities
plt.plot(N, pv.u_Si_e_maj(N), label='electron majority')
plt.plot(N, pv.u_Si_e_min(N), label='electron minority')
plt.plot(N, pv.u_Si_h_maj(N), label='hole majority')
plt.plot(N, pv.u_Si_h_min(N), label='hole minority')
plt.semilogx()
plt.legend(loc='upper right')
plt.title('carrier mobilities in silicon')
plt.xlabel('doping (cm**{-3})')  # add axis labels and plot title
plt.ylabel('mobility (cm**2Vs**{-1})µ²')

# plot the resisitivity
##figure(2)
plt.figure()
plt.plot(N, pv.resistivity_Si_n(N), label='n-type')
plt.plot(N, pv.resistivity_Si_p(N), label='p-type')
plt.loglog()
plt.legend()
plt.title('resisitivity in silicon as a function of doping')
plt.xlabel('doping (cm^{-3})')
plt.ylabel('resistivity (\u03A9cm)')
plt.xlim(1e12, 1e21)
plt.ylim(1e-4, 1e4)
# plt.show()
delta_n = pv.carrier_V(0.7, 1e15)
print(pv.impliedV(1e15, 6.7e15))
G = pv.current2gen(0.038)
tau = delta_n / G

print(delta_n, G, tau)

# plot the mobility models
N = np.logspace(14, 22)
mmobility = pv.masetti_mobility(N)
PC1Dmobility = pv.u_Si_e_maj(N)
plt.plot(N, mmobility, label='Masetti')
plt.plot(N, PC1Dmobility, label='PC1D')
plt.xlabel('doping (/cm³)')
plt.ylabel('Mobility (cm²/V·s)')
plt.loglog()
plt.ylim(10, 3000)
plt.xlim(1e14, 1e22)
plt.legend(loc=0)
plt.title('electron mobility in silicon')
