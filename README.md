# photovoltaic


photovoltaic is a library of python functions used in photovoltaics. The functions are simple enough to include in your code

Help Index: http://htmlpreview.github.io/?https://github.com/pvedu/photovoltaic/blob/master/html/photovoltaic.html


Typical usage:

    import photovoltaic as pv
    irradiance = pv.sun.blackbody_spectrum(800)
    print(irradiance)

This would print the blackbody irradiance at 800 nm with the default temperature of 6000 K in W/m2/nm.


## Installation


Installation is via pip from the pypi repositry. From a command propmpt:

    pip install photovoltaic

	
Tha above command should also install the latest scipy and numpy packages. They can also be installed directly with:

    pip install numpy

    pip install scipy

## Requirements

Known to work under plain vanilla Python 3.6 using the standard IDLE editor with Numpy and Scipy installed. The examples also make use of matplotlib. It should also work with the  various Python systems such as Anaconda Jupyter etc.


Anaconda includes a wealth of scientific packagkes and is available at: https://www.anaconda.com/download/ 

Standard Python is at https://www.python.org/downloads/

For the graphs, Matplotlib is needed in addition to the above numpy and scipy packages:
    pip install matplotlib



## Examples


Examples of using the library are available at:
https://github.com/pvedu/pvon

But its much more fun to go to: https://mybinder.org/v2/gh/pvedu/pvon/master
where you can run the code without installing anything.


## Other

f is used to mean from in some of the function names. For example:

nmfeV() converts the energy of a photon from electron volts to a nm.

This follows the conventions of other python functions such as strfdatetime.


The library is designed to be as simple as possible and an "algorithm that runs". While it is easier to install the whole library, it is also straighforward to cut/paste parts of the code.
