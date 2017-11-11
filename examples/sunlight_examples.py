import photovoltaic as pv
import numpy as np

print('declination',pv.declination(1))
print('Equation of time (min): ', pv.equation_of_time(1))
print('Time Correction (min): ', pv.time_correction(-3.71,140, 10))
print('Elevation, Azimuth (degrees): ', pv.elev_azi(-23.0116367279, -34, 13.27))
print('Sun Azimuth and Elevation (degrees)', pv.sun_position(1, -34, 140, 10, 14, 0))
print('Sunrise, Sunset (hours)', pv.sun_rise_set(-34, -23.0116367279, -43.71))
print('Direct light on a module (fraction)', pv.module_direct(94.7078, 36.9758, 0, 45))

'''
Code
print('declination',pv.declination(1))
print('Equation of time (min): ', pv.equation_of_time(1))
print('Time Correction (min): ', pv.time_correction(-3.71,140, 10))
print('Elevation, Azimuth (degrees): ', pv.elev_azi(-23.0116367279, -34, 13.27))
print('Sun Azimuth and Elevation (degrees)', pv.sun_position(1, -34, 140, 10, 14, 0))
print('Sunrise, Sunset (hours)', pv.sun_rise_set(-34, -23.0116367279, -43.71))
print('Direct light on a module (fraction)', pv.module_direct(94.7078, 36.9758, 0, 45))

Results
declination -23.0116367279
Equation of time (min):  -3.7051783234
Time Correction (min):  -43.71
Elevation, Azimuth (degrees):  (70.030717757890969, 298.39779695285648)
Sun Azimuth and Elevation (degrees) (70.01342676739516, 298.35884675563125)
Sunrise, Sunset (hours) (5.6187032686031166, 19.838296731396884)
Direct light on a module (fraction) 0.378945195148
'''
