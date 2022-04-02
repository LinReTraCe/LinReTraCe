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

from structure.dft import DFTcalculation
# from structure.model import Model
# from structure.btp   import BTPInterpolation

def h5output(outfile, escalc, btpinterp=None, peierls=False):
  '''
  Output function which takes an instance of an ElectronicStructure calculation
  and writes the internal data to an hdf5 file.
  Secondary argument is optional and contains the interpolated
  band velocities / band curvatures.
  In the case of models this can be computed explicitly.
  Third argument is also optional and functions as a switch
  to use the velocities squared as the diagonal optical elements (peierls approximation)
  '''

  # if not isinstance(escalc, ElectronicStructure):
  #   sys.exit('Provided electronic structure object for h5output is illegal.')

  # if btpinterp is not None:
  #   if not (isinstance(btpinterp, Model) or \
  #      isinstance(btpinterp, BTPInterpolation)):
  #     sys.exit('Provided interpolation object for h5output is illegal.')


  logger.info("Initializing output: {}".format(outfile))
  with h5py.File(outfile, 'w') as h5out:

    # identifier
    h5out.attrs['identifier'] = 'LRTCinput'

    # output of attributes
    h5out['.kmesh/nkp']           = escalc.nkp
    h5out['.kmesh/nkx']           = escalc.nkx
    h5out['.kmesh/nky']           = escalc.nky
    h5out['.kmesh/nkz']           = escalc.nkz
    h5out['.kmesh/weights']       = escalc.weights
    if type(escalc.weightsum) is not int:
      raise IOError('Weightsum must be integer')
    h5out['.kmesh/weightsum']     = escalc.weightsum
    if escalc.multiplicity.dtype != int:
      raise IOError('Multiplicity must be array of integers')
    h5out['.kmesh/multiplicity']  = escalc.multiplicity
    h5out['.kmesh/points']        = escalc.kpoints
    h5out['.kmesh/irreducible']   = escalc.irreducible

    h5out['.unitcell/volume']     = escalc.vol
    h5out['.unitcell/ndim']       = escalc.ndim
    h5out['.unitcell/dims']       = escalc.dims

    h5out['.bands/charge']        = escalc.charge
    h5out['.bands/energyBandMax'] = escalc.energyBandMax
    h5out['.bands/ispin']         = escalc.spins
    h5out['.bands/mu']            = escalc.mu

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


      # energies + derivatives
      h5out[prefix+'energies']          = escalc.energies[ispin]
      if btpinterp is not None:
        # we do not need to truncate here
        ''' these derivatives used to be used in LRTC: now not necessary anymore '''
        if False:
          h5out[prefix+'derivatives']     = btpinterp.velocities[ispin]
          h5out[prefix+'curvatures']      = btpinterp.curvatures[ispin]
        h5out[prefix+'momentsDiagonalBfield'] = btpinterp.BopticalDiag[ispin]


      if isinstance(escalc, DFTcalculation):
        if ispin == 0: # only warn once
          if escalc.opticdiag and btpinterp is not None:
            if peierls:
              logger.warning('CAREFUL: Overwriting intra-band optical elements from DFT calculation.')
            else:
              if escalc != btpinterp:
                logger.warning('CAREFUL: You are mixing DFT and Peierls optical elements.')

      # optical elements
      # use the peierls approximation
      # if a) specificially asked for
      #    b) its the only thing we have

      if not escalc.opticdiag and btpinterp is not None:
        peierls = True

      if (peierls and btpinterp is not None):
        if escalc.opticfull:
          ''' if the full elements are present we need to truncate accordingly '''
          h5out[prefix+'momentsDiagonal'] = btpinterp.opticalDiag[ispin][:,escalc.opticalBandMin:escalc.opticalBandMax,:]
          if ispin == 0:
            h5out['.bands/opticalBandMin']  = escalc.opticalBandMin + 1 # internal -> Fortran
            h5out['.bands/opticalBandMax']  = escalc.opticalBandMax
        else:
          ''' if they are not : use the full range from the interpolated data '''
          h5out[prefix+'momentsDiagonal'] = btpinterp.opticalDiag[ispin]
          if ispin == 0:
            h5out['.bands/opticalBandMin']  = btpinterp.opticalBandMin + 1 # internal -> Fortran
            h5out['.bands/opticalBandMax']  = btpinterp.opticalBandMax
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

      if btpinterp is not None and btpinterp.bopticfull:
        for ikp in range(escalc.nkp):
          h5out[prefix+'kPoint/{:010}/momentsBfield'.format(ikp+1)] = btpinterp.BopticalMoments[ispin][ikp]
