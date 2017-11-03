'''
Blackbody emmission as a function of wavelength.
pveducation.org S.Bowden 2016
http://pveducation.org/pvcdrom/properties-of-sunlight/blackbody-radiation
'''
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import photovoltaic as pv

wavelength = 600 # (µm) Python 3 allows unicode variables. λ is u03BB
print(pv.blackbody_spectrum(wavelength)) # just a single variable, using default T = 6000 K
print('Blackbody radiation is {:.2e} at {:.2f} µm using the default 6000 K'.format(pv.blackbody_spectrum(wavelength),wavelength)) # same but nicely formatted.

wavelength = np.linspace(200,3000,200) # 100 nm to 3 um, 200 points
temperature = 4500 # (K)

plt.plot(wavelength, pv.blackbody_spectrum(wavelength,temperature))
plt.xlabel('wavelength (µm)')       #  add axis labels and plot title
plt.ylabel('spectral irradiance (W m$^{-2}$ µm$^{-1}$)') # tex formating works in Matplotlib but not for all strings.
plt.title('Blackbody Radiation at '+str(temperature)+' K')
plt.ylim(0)
plt.xlim(0)
plt.show()
