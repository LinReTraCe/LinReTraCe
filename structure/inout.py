#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import warnings
import os
import sys
import logging
logger = logging.getLogger(__name__)

import numpy as np
with warnings.catch_warnings():
  warnings.filterwarnings("ignore",category=FutureWarning)
  import h5py

from structure.dft import DftCalculation
from structure.model import Model

def h5output(outfile, escalc, velcalc=None, peierls=False):
  '''
  Output function which takes an instance of an ElectronicStructure calculation
  and writes the internal data to an hdf5 file.

  Secondary argument is optional and contains the derivatives of the bands
  These are required to calculate the magnetic field optical elements.
  derivatives must therefore be either a BoltztrapInterpolation object or a Model object.

  Third argument is also optional and functions as a switch
  to use the velocities squared as the diagonal optical elements (peierls approximation)
  '''

  logger.info("Initializing output: {}".format(outfile))
  with h5py.File(outfile, 'w') as h5out:

    # identifier
    h5out.attrs['identifier'] = 'LRTCinput'

    # output of attributes
    h5out['.kmesh/nkp']           = int(escalc.nkp)
    h5out['.kmesh/nkx']           = int(escalc.nkx)
    h5out['.kmesh/nky']           = int(escalc.nky)
    h5out['.kmesh/nkz']           = int(escalc.nkz)
    if escalc.weights.dtype != np.float64:
      raise IOError('weights must be array of np.float64')
    h5out['.kmesh/weights']       = escalc.weights
    h5out['.kmesh/weightsum']     = int(escalc.weightsum)
    if escalc.multiplicity.dtype != int:
      raise IOError('Multiplicity must be array of integers')
    h5out['.kmesh/multiplicity']  = escalc.multiplicity
    if escalc.kpoints.dtype != np.float64:
      raise IOError('kpoints must be array of np.float64')
    h5out['.kmesh/points']        = escalc.kpoints
    h5out['.kmesh/irreducible']   = escalc.irreducible

    h5out['.unitcell/volume']     = float(escalc.vol)
    h5out['.unitcell/ndim']       = int(escalc.ndim)
    if escalc.dims.dtype != bool:
      raise IOError('dims must be array of boolean')
    h5out['.unitcell/dims']       = escalc.dims

    h5out['.bands/charge']        = float(escalc.charge)
    h5out['.bands/energyBandMax'] = int(escalc.energyBandMax)
    h5out['.bands/ispin']         = int(escalc.spins)
    h5out['.bands/mu']            = float(escalc.mu)

    for ispin in range(escalc.spins):
      if escalc.spins == 1:
        prefix = '/'
      else:
        prefix = '/up/' if ispin == 0 else '/dn/'

      # band gap information
      h5out['.bands/bandgap'+prefix+'gapped']      = escalc.gapped[ispin]
      if escalc.gapped[ispin]:
        h5out['.bands/bandgap'+prefix+'gapsize']   = escalc.gap[ispin]
        h5out['.bands/bandgap'+prefix+'cband']     = escalc.cb[ispin] + 1 # python -> Fortran
        h5out['.bands/bandgap'+prefix+'vband']     = escalc.vb[ispin] + 1
        h5out['.bands/bandgap'+prefix+'ene_cband'] = escalc.ecb[ispin]
        h5out['.bands/bandgap'+prefix+'ene_vband'] = escalc.evb[ispin]


      # energies
      h5out[prefix+'energies'] = escalc.energies[ispin]
      if isinstance(escalc, DftCalculation):
        if ispin == 0: # only warn once
          if escalc.opticdiag and velcalc is not None:
            if peierls:
              logger.warning('CAREFUL: Overwriting intra-band optical elements from DFT calculation.')
            else:
              if escalc != velcalc:
                logger.warning('CAREFUL: You are mixing DFT and Peierls optical elements.')

      # optical elements
      # use the peierls approximation
      # if a) specificially asked for
      #    b) its the only thing we have
      if not escalc.opticdiag and velcalc is not None:
        peierls = True

      if (peierls and velcalc is not None):
        if escalc.opticfull:
          ''' if the full elements are present we need to truncate accordingly '''
          h5out[prefix+'momentsDiagonal'] = velcalc.opticalDiag[ispin][:,escalc.opticalBandMin:escalc.opticalBandMax,:]
          if ispin == 0:
            h5out['.bands/opticalBandMin']  = escalc.opticalBandMin + 1 # internal -> Fortran
            h5out['.bands/opticalBandMax']  = escalc.opticalBandMax
        else:
          ''' if they are not : use the full range from the interpolated data '''
          h5out[prefix+'momentsDiagonal'] = velcalc.opticalDiag[ispin]
          if ispin == 0:
            h5out['.bands/opticalBandMin']  = velcalc.opticalBandMin + 1 # internal -> Fortran
            h5out['.bands/opticalBandMax']  = velcalc.opticalBandMax
      elif escalc.opticdiag:
        h5out[prefix+'momentsDiagonal'] = escalc.opticalDiag[ispin]
        if ispin == 0:
          h5out['.bands/opticalBandMin']  = escalc.opticalBandMin + 1 # internal -> Fortran
          h5out['.bands/opticalBandMax']  = escalc.opticalBandMax
      else:
        if ispin == 0: # only warn once
          # alibi variable so no errors occur
          h5out['.bands/opticalBandMin']  = 0
          h5out['.bands/opticalBandMax']  = 0
          logger.critical('No optical elements available. Use band interpolation (--interp) or provide them (--optic).')
          logger.critical('Output file is not able to produce transport results with linretrace.\n')

      if escalc.opticfull:
        for ikp in range(escalc.nkp):
          h5out[prefix+'kPoint/{:010}/moments'.format(ikp+1)] = escalc.opticalMoments[ispin][ikp]

      # no truncation necessary, we have the b-field quantities always on the full energy range
      if velcalc is not None:
        if velcalc.bopticdiag:
          h5out[prefix+'momentsDiagonalBfield'] = velcalc.BopticalDiag[ispin]
        if velcalc.bopticfull:
          for ikp in range(escalc.nkp):
            h5out[prefix+'kPoint/{:010}/momentsBfield'.format(ikp+1)] = velcalc.BopticalMoments[ispin][ikp]
