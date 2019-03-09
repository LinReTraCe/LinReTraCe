#! /usr/bin/env python

import numpy as np
import os
import h5py

def h5output(outfile, dftcalc, btpinterp=None):
  print("\nInitializing output: ", outfile)
  h5out = h5py.File(outfile, 'w')

  h5out['.kmesh/kpoints']      = dftcalc.kpoints
  h5out['.kmesh/nkp']          = dftcalc.nkp
  h5out['.kmesh/weights']      = dftcalc.weights
  h5out['.kmesh/weightsum']    = dftcalc.weightsum

  h5out['.unitcell/latLength'] = dftcalc.latLength
  h5out['.unitcell/latAngle']  = dftcalc.latAngle
  h5out['.unitcell/volume']    = dftcalc.vol

  h5out['.bands/charge']             = dftcalc.charge
  h5out['.bands/energyBandMax']      = dftcalc.energyBandMax
  h5out['.bands/opticalBandMin']     = dftcalc.opticalBandMin
  h5out['.bands/opticalBandMax']     = dftcalc.opticalBandMax

  # unpolarized
  if dftcalc.spins == 1:
    for ikp in range(dftcalc.nkp):
      h5out['/kpoint/{:06}/energies'.format(ikp+1)]          = dftcalc.energies[0][ikp]

      if btpinterp is not None:
        h5out['/kpoint/{:06}/derivatives'.format(ikp+1)]     = btpinterp.velocities[0][ikp]
        h5out['/kpoint/{:06}/curvatures'.format(ikp+1)]      = btpinterp.curvatures[0][ikp]

      if (dftcalc.optic):
        h5out['/kpoint/{:06}/moments'.format(ikp+1)]         = dftcalc.moments[0][ikp]

  # spin polarized
  elif dftcalc.spins == 2:
    for ikp in range(dftcalc.nkp):
      h5out['/up/kpoint/{:06}/energies'.format(ikp+1)]       = dftcalc.energies[0][ikp]
      h5out['/dn/kpoint/{:06}/energies'.format(ikp+1)]       = dftcalc.energies[1][ikp]

      if btpinterp is not None:
        h5out['/up/kpoint/{:06}/derivatives'.format(ikp+1)]  = btpinterp.velocities[0][ikp]
        h5out['/up/kpoint/{:06}/curvatures'.format(ikp+1)]   = btpinterp.curvatures[0][ikp]
        h5out['/dn/kpoint/{:06}/derivatives'.format(ikp+1)]  = btpinterp.velocities[1][ikp]
        h5out['/dn/kpoint/{:06}/curvatures'.format(ikp+1)]   = btpinterp.curvatures[1][ikp]

      if (dftcalc.optic):
        h5out['/up/kpoint/{:06}/moments'.format(ikp+1)]      = dftcalc.moments[0][ikp]
        h5out['/dn/kpoint/{:06}/moments'.format(ikp+1)]      = dftcalc.moments[1][ikp]

  print("Output file successfully created.")
  print("File size: ", os.stat(outfile).st_size/1000000, " MB.")
