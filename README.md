# CS173

<a href="https://github.com/sceccode/cs173.git"><img src="https://github.com/sceccode/cs173/wiki/images/cs173_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/cs173)

CS173

This model is the velocity model used to produce CyberShake Study 17.3 in Central California. 
The model was constructed by tiling together the CCA-06, CVM-S4.26, USGS Bay Area v8.3.0, and 
the CVM-S4.26 1D model.  The minimum Vs is 900 m/s, Vp is 1600 m/s, and density is 1600 kg/m3. 
Smoothing was applied within 10 km of tiling model interfaces. The model has a resolution of 
175m and extends down to a depth of 50.4 km'.

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.7 or higher.

## Library

The library ./lib/libcs173.a may be statically linked into any
user application. Also, if your system supports dynamic linking,
you will also have a ./lib/libcs173.so file that can be used
for dynamic linking. The header file defining the API is located
in ./include/cs173.h.


