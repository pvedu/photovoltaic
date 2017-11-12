# -*- coding: utf-8 -*-
"""
Declination angle. Shows the equivalence of the three simple formula.
and the PSA algorithm

"""
import matplotlib.pyplot as plt
import numpy as np
import photovoltaic as pv


def declinationPSA(dayNo):
    '''converted from C++ code at www.psa.es/sdg/sunpos.htm Please check website for latest code'''
    dElapsedJulianDays = dayNo
    dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
    dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
    dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
    dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) + 0.00034894 * np.sin(
        2 * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * np.cos(dOmega)
    dSin_EclipticLongitude = np.sin(dEclipticLongitude)
    dDeclination = np.arcsin(np.sin(dEclipticObliquity) * dSin_EclipticLongitude)
    return (np.degrees(dDeclination))


dayNo = 263

dec1 = 23.45 * pv.sind(360 / 365 * (dayNo + 284))
dec2 = -23.45 * pv.cosd(360 / 365 * (dayNo + 10))
dec3 = 23.45 * pv.sind(360 / 365 * (dayNo - 81))
dec4 = declinationPSA(dayNo)

print('declination (degrees) at day ' + str(dayNo))
print(dec1, dec2, dec3, dec4)

dayNo = np.linspace(0, 365, 365)  # day number from 1 to 365

dec1 = 23.45 * pv.sind(360 / 365 * (dayNo + 284))
dec2 = -23.45 * pv.cosd(360 / 365 * (dayNo + 10))
dec3 = 23.45 * pv.sind(360 / 365 * (dayNo - 81))
dec4 = declinationPSA(dayNo)

plt.xlabel('day number')
plt.ylabel('declination (degrees)')
plt.xlim(1, 365)
plt.grid(True)
plt.plot(dayNo, dec1, label=r'$\delta=-23.45 \cos(\frac{360}{365}(d+10))$', color='r')
plt.plot(dayNo, dec2, label=r'$\delta=23.45 \sin(\frac{360}{365}(d+284))$')
plt.plot(dayNo, dec3, label=r'$\delta=23.45 \sin(\frac{360}{365}(d-81))$')
plt.plot(dayNo, dec4, label='PSA algorithm')
plt.legend(loc='lower center')
plt.savefig('plot_declination.png')
plt.show()
