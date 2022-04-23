#!/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import logging
logger = logging.getLogger(__name__)

import numpy as np
try:
  import spglib
  spglib_exist = True
except ImportError:
  spglib_exist = False

# python 2 & 3 compatible usage of abstract base classes
if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

import structure.symmetries.C1
from   structure.es import ElectronicStructure

class Model(ElectronicStructure, ABC):
  '''
  Abstract Base Class for generic model Hamiltonians.
  Currently we only use this class to generate unpolarized polarizations:
  -> weightsum = 2; spins = 1
  '''

  def __init__(self, nkx, nky, nkz):
    super(Model, self).__init__()

    self.nkx            = nkx
    self.nky            = nky
    self.nkz            = nkz
    self.nkp            = self.nkx*self.nky*self.nkz
    self.weightsum      = 2  # always spin-unpolarized (paramagnetic)
    self.spins          = 1

    ''' this should always be there for models '''
    self.opticdiag      = True
    self.bopticdiag     = True

    self.velocities     = []
    self.curvatures     = []

    self._defineDimensions()

  @abc.abstractmethod
  def computeData(self):
    '''
    Abstract method which is always called when computing the
    electronic structure and optical elements of the model
    '''
    pass

  def _setupArrays(self, ortho):
    '''
    Setup energy, derivative, curvature and optical element arrays
    Currently set up to support non-orthogonal optical elements without spin-orbit coupling
    '''
    logger.info('Setting up arrays with {} band{}.'.format(self.energyBandMax, 's' if self.energyBandMax > 1 else ''))

    if ortho:
      ioptical = 3
    else:
      ioptical = 9

    ''' generic array structure for spin unpolarized calculations '''
    self.energies        = [np.zeros((self.nkp, self.energyBandMax,), dtype=np.float64)]
    self.opticalMoments  = [np.zeros((self.nkp, self.energyBandMax, self.energyBandMax, ioptical), dtype=np.float64)]
    self.opticalDiag     = [np.zeros((self.nkp, self.energyBandMax, ioptical), dtype=np.float64)]
    self.BopticalDiag    = [np.zeros((self.nkp, self.energyBandMax, 3, 3, 3), dtype=np.complex128)]
    self.BopticalMoments = [np.zeros((self.nkp, self.energyBandMax, self.energyBandMax, 3, 3, 3), dtype=np.complex128)]

    ''' velocities and curvatures depend strongly on the setup
        hence they will be initiated according to the needs of the specific class '''

  def setDiagonal(self,value):
    logger.info('Setting diagonal optical elements to value: {}'.format(value))
    self.opticalDiag[0][:,:,:] = float(value)
    self.BopticalDiag[0][:,:,:,:,:] = float(value)

  def setOffDiagonal(self,value):
    if self.energyBandMax > 1:
      logger.info('Setting offdiagonal optical elements to value: {}'.format(value))
      ''' keep the intra band elements here '''
      self.opticalMoments[0][:,:,:,:] = float(value)
      self.opticalMoments[0][:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
        = self.opticalDiag[0][:,:,:]

      self.BopticalMoments[0][:,:,:,:,:,:] = float(value)
      self.BopticalMoments[0][:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
        = self.BopticalDiag[0][:,:,:,:,:]
      logger.info('Setting offdiagonal output switch')
      self.opticfull  = True
      self.bopticfull = True
    else:
      logger.warning('Cannot set offdiagonal optical elements since there is only one band')
