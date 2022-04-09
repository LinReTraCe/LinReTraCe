#!/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np

from   structure.model import Model

class Quadratic(Model):
  '''
  Class to create a quadratic dispersion
  format compared to the DFT version.
  Here we immediately combine the DFT class and the BTP interpolation class.
  Hence for the output we simply supply the same object for both the dftcalc
  and the btpinterp argument.
  '''

  def __init__(self, spacing, nkx=1, nky=1, nkz=1):
    super(Quadratic, self).__init__(nkx,nky,nkz)
    raise NotImplementedError('Quadratic models have been deprecated.')
    self.spacing = spacing

    self._defineDimensions()
    self.vol = self.spacing**self.ndim
    self._setupKmesh()

  def computeData(self, e0, mass, charge):
    self.e0      = np.array(e0, dtype=np.float64)
    self.mass    = np.array(mass, dtype=np.float64)
    self.charge  = charge

    if self.e0.shape[0] != self.mass.shape[0]:
      raise IOError('Offset energies and mass parameters must have same dimensionality')

    self.energyBandMax  = self.e0.shape[0]
    self.opticalBandMax = self.energyBandMax

    if self.charge <= 0 or self.charge >= self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands')

    self._setupArrays()
    self._computeDispersion()
    self._calcFermiLevel()

  def _setupKmesh(self):
    '''
    Setup kmesh from [-1,1) pi/a
    '''

    logger.info('Setting up kmesh with {} kpoints.'.format(self.nkp))
    logger.info('nkx = {}'.format(self.nkx))
    logger.info('nky = {}'.format(self.nky))
    logger.info('nkz = {}'.format(self.nkz))

    if self.nkx > 1:
      self._kmeshx = np.linspace(-1,1,self.nkx,endpoint=False)
    else:
      self._kmeshx = np.array([0.0], dtype=np.float64) # this is necessary here

    if self.nky > 1:
      self._kmeshy = np.linspace(-1,1,self.nky,endpoint=False)
    else:
      self._kmeshy = np.array([0.0], dtype=np.float64) # this is necessary here

    if self.nkz > 1:
      self._kmeshz = np.linspace(-1,1,self.nkz,endpoint=False)
    else:
      self._kmeshz = np.array([0.0], dtype=np.float64) # this is necessary here

    self.kpoints = np.empty((self.nkp,3), dtype=np.float64)
    for ikx in range(self.nkx):
      for iky in range(self.nky):
        for ikz in range(self.nkz):
          self.kpoints[ikz + iky*self.nkz + ikx*self.nky*self.nkz,:] = \
            np.array([self._kmeshx[ikx], self._kmeshy[iky], self._kmeshz[ikz]])


  def _computeDispersion(self):
    '''
    create the quadradic energy dispersion
    '''

    for iband in range(self.energyBandMax):
      ek = np.zeros((self.nkx, self.nky, self.nkz),     dtype=np.float64)
      vk = np.zeros((self.nkx, self.nky, self.nkz, 3) , dtype=np.float64)
      ck = np.zeros((self.nkx, self.nky, self.nkz, 6),  dtype=np.float64)

      # dispersion
      # e(k) = k**2 / 2m   --- k = -pi / a ... pi/a
      # internally k runs from -1 ... 1 (not included)

      ek += ((self._kmeshx*np.pi/self.spacing)[:,None,None]**2) \
         +  ((self._kmeshy*np.pi/self.spacing)[None,:,None]**2) \
         +  ((self._kmeshz*np.pi/self.spacing)[None,None,:]**2)

      ek /= (2.0*self.mass[iband])

      # there is no dimensionality offset
      # due to the way we set up the meshes

      # actual suppiled offset
      ek += self.e0[iband]

      # first derivative -> velocities
      vk[...,0] += (self._kmeshx*np.pi/self.spacing)[:,None,None] / self.mass[iband]
      vk[...,1] += (self._kmeshy*np.pi/self.spacing)[None,:,None] / self.mass[iband]
      vk[...,2] += (self._kmeshz*np.pi/self.spacing)[None,None,:] / self.mass[iband]

      ck[...,0] += 1./self.mass[iband]
      ck[...,1] += 1./self.mass[iband]
      ck[...,2] += 1./self.mass[iband]
      # 3 , 4  and 6
      # xy, xz and yz is identical to 0 in an orthogonal lattice!

      # now we save all this data into the preallocated arrays
      self.energies[0][:,iband] = ek.flatten()

      self.velocities[0][:,iband,0] = vk[...,0].flatten()
      self.velocities[0][:,iband,1] = vk[...,1].flatten()
      self.velocities[0][:,iband,2] = vk[...,2].flatten()

      self.curvatures[0][:,iband,0] = ck[...,0].flatten()
      self.curvatures[0][:,iband,1] = ck[...,1].flatten()
      self.curvatures[0][:,iband,2] = ck[...,2].flatten()
      self.curvatures[0][:,iband,3] = ck[...,3].flatten()
      self.curvatures[0][:,iband,4] = ck[...,4].flatten()
      self.curvatures[0][:,iband,5] = ck[...,5].flatten()

      self.opticalMoments[0][:,iband,iband,0] = self.velocities[0][:,iband,0]**2
      self.opticalMoments[0][:,iband,iband,1] = self.velocities[0][:,iband,1]**2
      self.opticalMoments[0][:,iband,iband,2] = self.velocities[0][:,iband,2]**2

      self.opticalDiag[0][:,iband,0] = self.opticalMoments[0][:,iband,iband,0]
      self.opticalDiag[0][:,iband,1] = self.opticalMoments[0][:,iband,iband,1]
      self.opticalDiag[0][:,iband,2] = self.opticalMoments[0][:,iband,iband,2]



