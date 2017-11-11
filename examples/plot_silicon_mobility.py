# -*- coding: utf-8 -*-
"""
mobility as a function of doping with different models

@author: Stuart Bowden
"""
import numpy as np
import matplotlib.pyplot as plt
import photovoltaic as pv

# plot the mobility models
N_D = np.logspace(14, 22) # log space from 1e14 to 1e22
mobility_masetti = pv.mob_masetti_phos(N_D)
mobility_thurber = pv.mob_thurber(N_D, False)
plt.plot(N_D, mobility_masetti, label='Masetti')
plt.plot(N_D, mobility_thurber, label='Thurber')
plt.xlabel('doping (/cm³)')
plt.ylabel('Mobility (cm²/V·s)')
plt.loglog()
plt.ylim(10, 3000)
plt.xlim(1e14, 1e22)
plt.legend(loc=0)
plt.title('electron majority mobility in silicon')
plt.savefig('plot_silicon_mobility.png')
plt.show()
