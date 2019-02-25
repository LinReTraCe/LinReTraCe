#! /usr/bin/env python

import btp
import dft
from inout import h5output

import numpy as np
import h5py

# calc = dft.vaspcalculation('../PbTe.vasp.collinear')

# calc = dft.vaspcalculation('../PbTe.vasp.collinear')
# calc.readData()

calc = dft.w2kcalculation('../aFe-so/', optic=True)
calc.readData()

interp = btp.BTPInterpolation(calc, niter=10)
interp.interpolate()

h5output('preprocessed_data.hdf5', calc, interp)


# print(interp.energies)
# print(interp.velocities)
# print(interp.curvatures)

# print(calc.energies[0].shape)
# print(calc.energies[1].shape)
# print(calc.moments[0].shape)
# print(calc.moments[1].shape)
