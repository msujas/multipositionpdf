A library for processing total-scattering data with multiple detector positions on BM31.

Only tested with 1D poni interpolation so far.

Usage: Must make ponis for a handful of calibration images (5 or so). These are loaded together with detector y (and optional z) positions into the PoniData class. The ponis for the rest of the positions are then interpolated in 1 or 2 dimensions (only 1D tested so far) with the PoniList class. CBF files are loaded into the FilePoni class together with their detector y (and z) positions, and the PoniList, and a interpolated poni is calculated for it, then it is integrated.

The MultiFile class then takes a list of all the FilePoni objects. It regrids the cake arrays to the same axes and merges them together, saves a merged cake file, and 1D merged patterns.