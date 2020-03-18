"""Simple slider to demonstrate the effect of an antir-reflection coating thickness and refractie index on the reflectivity"""

import matplotlib.pyplot as plt
import numpy as np
import photovoltaic as pv
from matplotlib.widgets import Slider

n0 = 1

nSemi = 3.5
Eg = 1.12 # (eV)
Wg = pv.eVfnm(Eg) # (nm) the wavelength corresponding to Eg

# read standard specta
wavelength, AMO, AM15G, AM15D = pv.sun.solar_spectra()
# truncate to wavelengths below bandgap
wavelength_Eg = wavelength[wavelength < Wg] # wavelengths with photon greater than the bandgap
AM15G_Eg = AM15G[wavelength < Wg] # corresponding spectrum

# current to band gap
photons = AM15G_Eg / (pv.eVfnm(wavelength_Eg))
Jsc_total = np.trapz(photons, wavelength_Eg)
print('Jsc for photons above bandgap: {:.1f} A/m2'.format(Jsc_total))

# setup the plot
fig, ax = plt.subplots(figsize=(7, 5))
plt.subplots_adjust(left=0.15, bottom=0.3, top=0.97, right=0.97)
plt.axis([250, 1200, 0, 0.40])
plt.xlabel('wavelength (nm)')
plt.ylabel('reflectivity')
txt = plt.text(300, 0.35, 'refl', size=12)
# define the plot with dummy values for Y
l, = plt.plot(wavelength_Eg, wavelength_Eg, lw=2, color='red')

# setup the slider. In matplot lib its a hacked plot
ax_s1 = plt.axes([0.25, 0.15, 0.6, 0.03], fc='white')
s1 = Slider(ax_s1, 'thickness', 0.0, 200, valinit=100, valfmt='%0.0f nm')  # old style format
s1_xticks = np.arange(0.0, 200, 50)
ax_s1.xaxis.set_visible(True)
ax_s1.set_xticks(s1_xticks)

ax_s2 = plt.axes([0.25, 0.06, 0.6, 0.03], fc='white')
s2 = Slider(ax_s2, 'index', 1.0, 3.0, valinit=1.9, valfmt='%0.1f nm')  # old style format
s2_xticks = np.arange(1.0, 3.0, 0.5)
ax_s2.xaxis.set_visible(True)
ax_s2.set_xticks(s2_xticks)


def update(val):  # why do we need a val?
    t1 = s1.val
    n1 = s2.val
    refl = pv.optic.ARC_refl(wavelength_Eg, n0, n1, nSemi, t1) / 100
    Jsc_loss = np.trapz(photons * refl, wavelength_Eg)  # A/m
    R_weighted = Jsc_loss / Jsc_total
    txt.set_text('weighted reflectance {:.3f}'.format(R_weighted))
    l.set_ydata(refl)
    fig.canvas.draw_idle()


update(0)  # initial plot

s1.on_changed(update)
s2.on_changed(update)
plt.show()
