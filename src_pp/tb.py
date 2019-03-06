#! /usr/bin/env python

import dft
import numpy as np
import sys

class tightbinding(object):
  ''' Class to create tight binding band structures with an identical
  format compared to the DFT version.
  Here we immediately combine the DFT class and the BTP interpolation class.
  Hence for the output we simply supply the same object for both the dftcalc
  and the btpinterp argument. '''

  def __init__(self, nbands, nkx=1, nky=1, nkz=1):
    self.optic          = True
    self.nkx            = nkx
    self.nky            = nky
    self.nkz            = nkz
    self.nkp            = nkx*nky*nkz
    self.charge         = None
    self.weights        = np.full((self.nkp,), fill_value=2./self.nkp, dtype=np.float64)
    self.weightsum      = 2.
    self.energyBandMax  = nbands
    self.opticalBandMin = 1
    self.opticalBandMax = nbands
    self.spins          = 1
    self.energies       = []
    self.derivatives    = []
    self.curvatures     = []
    self.moments        = []

    self._defineDimension()
    self._setupKmesh()
    self._setuparrays()

  def _defineDimension(self):
    self.ndim = 0
    for i in [self.nkx, self.nky, self.nkz]:
      if i < 1:
        raise ValueError("Number of kpoints in each direction have to be positive")
      if i > 1:
        self.ndim += 1

  def _setupKmesh(self):
    self._kmeshx = np.linspace(0,1,self.nkx,endpoint=False)
    self._kmeshy = np.linspace(0,1,self.nky,endpoint=False)
    self._kmeshz = np.linspace(0,1,self.nkz,endpoint=False)

  def _setuparrays(self):
    self.energies.append(np.zeros((self.nkp, self.energyBandMax,), dtype=np.float64))
    self.derivatives.append(np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64))
    self.curvatures.append(np.zeros((self.nkp, self.energyBandMax, 3, 3), dtype=np.float64))
    self.moments.append(np.zeros((self.nkp, self.energyBandMax, self.energyBandMax, 6), dtype=np.complex128))

  def computeData(self, e0, dist, hopping):
    ''' create the energy dispersion from the following arrays:
      e0[nbands]
      dist[number of nearest neighbours]
      hopping[number of nearest neighbours, nbands]

      1 nearest neighbor : t
      2 nearest neighbors: t, t'
      etc.

      e0 : 0 point energies of the bands
      dist : distances normalized to the lattice constant a
             i.e. the first entries should always be 1
      hopping : hopping parameters t in the Hubbard Hamiltonian
                one minus sign already included
      '''

    nnn = len(dist) # number of nearest neighbors
    for iband in range(self.energyBandMax):
      ek = 0 # energy
      vk = 0 # velocity
      ck = 0 # curvature
      for i in range(nnn):
        ek += -2. * hopping[i,iband] * (np.cos(self._kmeshx*2*np.pi * dist[i])[:,None,None] \
                                     +  np.cos(self._kmeshy*2*np.pi * dist[i])[None,:,None] \
                                     +  np.cos(self._kmeshz*2*np.pi * dist[i])[:,None,None])

        vk += 2. * hopping[i,iband] * dist[i] * (np.sin(self._kmeshx*2*np.pi * dist[i])[:,None,None] \
                                              +  np.sin(self._kmeshy*2*np.pi * dist[i])[None,:,None] \
                                              +  np.sin(self._kmeshz*2*np.pi * dist[i])[:,None,None])

        ck += 2. * hopping[i,iband] * dist[i]**2 * (np.cos(self._kmeshx*2*np.pi * dist[i])[:,None,None] \
                                                 +  np.cos(self._kmeshy*2*np.pi * dist[i])[None,:,None] \
                                                 +  np.cos(self._kmeshz*2*np.pi * dist[i])[:,None,None])

      temp += e0[iband]



if __name__ == '__main__':
  tb = tightbinding(2, 5, 5, 1)
