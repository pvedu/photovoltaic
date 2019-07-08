============
photovoltaic
============

photovoltaic is a library of python functions used in photovoltaics. Individual functions are in photovoltaic.html, which can be viewed at:
http://htmlpreview.github.io/?https://github.com/trautsned/photovoltaic/blob/master/photovoltaic.html


Typical usage:


    import photovoltaic as pv

    irradiance = pv.sun.blackbody_spectrum(800)

    print(irradiance)

This would print the blackbody irradiance at 800 nm with the default temperature of 6000 K in W/m2/nm.


Installation
---------------

Installation is via pip from the pypi repositry. From a command propmpt:

    pip install photovoltaic
	
Tha above command should also install the latest scipy and numpy packages. They can also be installed directly with:

	pip install numpy
	
	pip install scipy

Requirements
------------
Known to work under plain vanilla Python 3.6 using the standard IDLE editor with Numpy and Scipy installed. The examples also make use of matplotlib. It should also work with the  various Python systems such as Anaconda Jupyter etc.


Anaconda includes a wealth of scientific packagkes and is available at: https://www.anaconda.com/download/ 

Standard Python is at https://www.python.org/downloads/

For the graphs, Matplotlib is needed in addition to the above numpy and scipy packages:
	pip install matplotlib



Examples
--------

There are many more examples on github at:
https://github.com/trautsned/photovoltaic/tree/master/examples

Other
-----
f is used to mean from in some of the function names. For example:

nmfeV() converts the energy of a photon from electron volts to a nm.

This follows the conventions of other python functions such as strfdatetime.


The library is designed to be as simple as possible and an "algorithm that runs". As such, most of the library is unashamedly procedural for simplicity. The syntax of object orientated code varies from language to language.
