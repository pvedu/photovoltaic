============
photovoltaic
============

photovoltaic is a library of python functions used in photovoltaics. Individual functions are in photovoltaic.html, shich can be viewed at:
http://htmlpreview.github.io/?https://github.com/trautsned/photovoltaic/blob/master/photovoltaic.html

Typical usage:


    import photovoltaic as pv

    irradiance = pv.blackbody_spectrum(800)

    print(irradiance)

This would print the blackbody irradiance at 800 nm with the default temperature of 6000K in W/m2/nm.


Installation
---------------

Installation is via pip and pypi:

    pip install photovoltaic

Requirements
------------
Known to work under plain vanilla Python 3.6 using the standard IDLE editor with Numpy and Scipy installed. The examples also make use of matplotlib. It should also work with the  various Python systems such as Anaconda Jupyter etc.

Examples
--------

There are many more examples on github at:
https://github.com/trautsned/photovoltaic/tree/master/examples
