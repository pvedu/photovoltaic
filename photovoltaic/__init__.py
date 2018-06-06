#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Basic photovoltaic functions
Provides the formulas and equations typically used in an introductory photocoltaics textbook.

Typical solar units are used, NOT SI units. The units are denoted in parenthesis on the comment lines.
wavelength (nm)
Energy of  photon (eV)
semiconductor dimensions (cm)
degrees instead of radians.
Temperature of 298.15 K (25 degC) not 300 K

The first line on all input files is ignored to allow for column headers
# denotes a comment in input files and is ignored.

Contributions by: sgbowden, richter, heraimenka, jhul etc

Variables and acronyms
ARC - anti-reflection coating
"""
__version__ = '0.1.5'

#import the submodules

from . import optic
from . import sun
from . import si
from . import money
from . import cell

#import the basic definitions
from .core import *