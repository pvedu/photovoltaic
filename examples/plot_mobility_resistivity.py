import photovoltaic as pv
import numpy as np
import matplotlib.pyplot as plt


N_D = np.logspace(12, 21, 200)  # sweep the doping in the base (cm-3)

# plot the mobilities
plt.plot(N_D, pv.u_Si_e_maj(N_D), label='electron majority')
plt.plot(N_D, pv.u_Si_e_min(N_D), label='electron minority')
plt.plot(N_D, pv.u_Si_h_maj(N_D), label='hole majority')
plt.plot(N_D, pv.u_Si_h_min(N_D), label='hole minority')
plt.semilogx()
plt.legend(loc='upper right')
plt.title('carrier mobilities in silicon')
plt.xlabel('doping (cm**{-3})')  # add axis labels and plot title
plt.ylabel('mobility (cm**2Vs**{-1})µ²')

# plot the resisitivity
plt.figure()
plt.plot(N_D, pv.resistivity_Si_n(N_D), label='n-type')
plt.plot(N_D, pv.resistivity_Si_p(N_D), label='p-type')
plt.loglog()
plt.legend()
plt.title('resisitivity in silicon as a function of doping')
plt.xlabel('doping (cm^{-3})')
plt.ylabel('resistivity (ohm cm)')
plt.xlim(1e12, 1e21)
plt.ylim(1e-4, 1e4)
# plt.show()
delta_n = pv.implied_carrier(0.7, 1e15)
print(pv.impliedV(1e15, 6.7e15))
G = pv.current2gen(0.038)
tau = delta_n / G
print(delta_n, G, tau)

# plot the mobility models
N_D = np.logspace(14, 22)
mmobility = pv.masetti_mobility(N_D)
PC1Dmobility = pv.u_Si_e_maj(N_D)
plt.plot(N_D, mmobility, label='Masetti')
plt.plot(N_D, PC1Dmobility, label='PC1D')
plt.xlabel('doping (/cm³)')
plt.ylabel('Mobility (cm²/V·s)')
plt.loglog()
plt.ylim(10, 3000)
plt.xlim(1e14, 1e22)
plt.legend(loc=0)
plt.title('electron mobility in silicon')
plt.show()