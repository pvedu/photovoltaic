import matplotlib.pyplot as plt
import numpy as np
import photovoltaic as pv

print('reading spectrum')
wavelength, AM0, AM15G, AM15D = pv.sun.solar_spectra()

#spectral irradiance vs wavelength

plt.figure('plot_solar_reference_spectra')
plt.plot(wavelength,AM0, label="AM0")
plt.plot(wavelength,AM15G, label="AM1.5 Global")
plt.plot(wavelength, AM15D, label="AM1.5 Direct")
plt.legend(loc='upper right')
plt.ylim(0, 2.2)
plt.xlim(200, 2000)
plt.grid(True)
plt.xlabel('wavelength (nm)')       #  add axis labels and plot title
plt.ylabel('spectral irradiance (W m$^{-2}$ nm$^{-1}$)')
plt.title('Solar Reference Spectra')
plt.savefig('plot_solar_reference_spectra.png')

plt.show()
