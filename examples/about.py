"Set of functions to give version, file path etc"

import photovoltaic as pv
print(dir(pv))
print('Doc String: ',pv.__doc__)
print(pv.__path__)
print(pv.__file__)
print('Version: ',pv.__version__)
print(pv.sind(45))