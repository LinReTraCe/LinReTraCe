#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import logging
import collections
logger = logging.getLogger(__name__)

import numpy as np

import structure.symmetries.C1

import structure.symmetries.D1x
import structure.symmetries.D1y
import structure.symmetries.D1z

import structure.symmetries.C2xy
import structure.symmetries.C2xz
import structure.symmetries.C2yz

import structure.symmetries.C4xy
import structure.symmetries.C4xz
import structure.symmetries.C4yz

import structure.symmetries.D2H

import structure.symmetries.D4Hxy
import structure.symmetries.D4Hxz
import structure.symmetries.D4Hyz

import structure.symmetries.OH

from   structure import es
from   structure.aux import levicivita

# python 2 & 3 compatible usage of abstract base classes
if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

from structure.es import ElectronicStructure

class Model(ElectronicStructure, ABC):
  '''
  Abstract Base Class for generic model Hamiltonians.
  '''

  def __init__(self, nkx, nky, nkz):
    super(Model, self).__init__()

    self.optic          = False
    self.opticdiag      = True
    self.opticfull      = False

    self.nkx            = nkx
    self.nky            = nky
    self.nkz            = nkz
    self.nkp            = self.nkx*self.nky*self.nkz
    self.charge         = None
    self.weights        = None
    self.weightsum      = 2. # always spin-unpolarized (paramagnetic)
    self.multiplicity   = None

    self.energyBandMax  = 1
    self.opticalBandMin = 0
    self.opticalBandMax = 1
    self.spins          = 1

    self.velocities     = []
    self.curvatures     = []

  @abc.abstractmethod
  def computeData(self):
    pass

  def _defineDimension(self):
    '''
    Count dimension as every k-axis with more than one k-point
    1-dimension k-axis is enforced to be kx
    2-dimension k-axis are enforced to be kx and ky
    '''
    self.ndim = 0
    self.dims = []
    for i in [self.nkx, self.nky, self.nkz]:
      if i < 1:
        raise ValueError("Number of kpoints in each direction have to be positive")
      if i > 1:
        self.ndim += 1
        self.dims.append(True)
      else:
        self.dims.append(False)
    self.dims = np.array(self.dims)
    logger.info('Detected {} dimensions.'.format(self.ndim))

  def _defineSymmetries(self):
    '''
    Define symmetry and inverse symmetry operations necessary for the construction
    of an irreducible k-grid
    In this parent class we set only the unit.
    '''
    self.nsym     = structure.symmetries.C1.nsym
    self.symop    = structure.symmetries.C1.symop
    self.invsymop = structure.symmetries.C1.invsymop

  def _setupArrays(self):
    '''
    Setup energy, derivative, curvature and optical element arrays
    Currently set up to support non-orthogonal optical elements without spin-orbit coupling
    '''
    logger.info('Setting up arrays with {} bands.'.format(self.energyBandMax))
    self.energies       = [np.zeros((self.nkp, self.energyBandMax,), dtype=np.float64)]
    self.velocities     = [np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64)]
    self.curvatures     = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalMoments = [np.zeros((self.nkp, self.energyBandMax, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalDiag    = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.BopticalDiag   = [np.zeros((self.nkp, self.energyBandMax, 3, 3, 3), dtype=np.float64)]

  def setDiagonal(self,value):
    logger.info('Setting diagonal optical elements to value: {}'.format(value))
    self.opticalDiag[0][:,:,:] = float(value)
    self.BopticalDiag[0][:,:,:,:,:] = float(value)

  def setOffDiagonal(self,value):
    if self.energyBandMax > 1:
      logger.info('Setting offdiagonal optical elements to value: {}'.format(value))
      self.opticalMoments[0][:,:,:,:] = float(value)
      self.opticalMoments[0][:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
        = self.opticalDiag[0][:,:,:]

      logger.info('Setting offdiagonal output switch')
      self.opticfull = True
    else:
      logger.warning('Cannot set offdiagonal optical elements since there is only one band')


class quadratic(Model):
  '''
  Class to create a quadratic dispersion
  format compared to the DFT version.
  Here we immediately combine the DFT class and the BTP interpolation class.
  Hence for the output we simply supply the same object for both the dftcalc
  and the btpinterp argument.
  '''

  def __init__(self, spacing, nkx=1, nky=1, nkz=1):
    super(quadratic, self).__init__(nkx,nky,nkz)
    self.spacing = spacing

    self._defineDimension()
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


class tightbinding(Model):
  '''
  Tightbinding class to create (ir-)reducible k-grids
  for orthogonal lattices in 1 / 2 / 3 dimensions
  Main class that gets interfaces with ltb script
  '''

  def __init__(self, ax=1, ay=1, az=1, nkx=1, nky=1, nkz=1, irreducible=True, kshift=False):
    super(tightbinding, self).__init__(nkx,nky,nkz)
    self.ax          = ax
    self.ay          = ay
    self.az          = az
    self.spacing     = [self.ax, self.ay, self.az]
    self.irreducible = irreducible  # generate irreducible grid instead of reducible
    self.kshift      = kshift       # shift by half a k-point to avoid Gamma point

    self._defineDimension() # method from parent class
    self._defineSymmetries()
    self._setupKmesh()      # setup reducible and reduce to irreducible

    self.vol = self.ax*self.ay*self.az

  def _setupArrays(self):
    '''
    Overwrite parent method in order to restrict to orthogonal arrays
    '''
    logger.info('Setting up arrays with {} bands.'.format(self.energyBandMax))
    self.energies       = [np.zeros((self.nkp, self.energyBandMax,), dtype=np.float64)]
    self.velocities     = [np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64)]
    self.curvatures     = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalMoments = [np.zeros((self.nkp, self.energyBandMax, self.energyBandMax, 3), dtype=np.float64)]
    self.opticalDiag    = [np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64)]
    self.BopticalDiag   = [np.zeros((self.nkp, self.energyBandMax, 3, 3, 3), dtype=np.float64)]

  def _defineSymmetries(self):
    '''
    Define symmetries for all possible orthogonal unit cells
    All these symmetry properties are hardcoded
    '''

    if self.ndim == 0:
      sym = structure.symmetries.C1
    if self.ndim == 1:
      if self.nkx != 1:
        sym = structure.symmetries.D1x
      elif self.nky != 1:
        sym = structure.symmetries.D1y
      elif self.nkz != 1:
        sym = structure.symmetries.D1z

    if self.ndim == 2:
      if self.nkz == 1:
        if self.ax == self.ay:
          sym = structure.symmetries.C4xy
        else:
          sym = structure.symmetries.C2xy
      if self.nky == 1:
        if self.ax == self.az:
          sym = structure.symmetries.C4xz
        else:
          sym = structure.symmetries.C2xz
      if self.nkx == 1:
        if self.ay == self.az:
          sym = structure.symmetries.C4yz
        else:
          sym = structure.symmetries.C2yz

    if self.ndim == 3:
      if self.ax == self.ay and self.ax == self.az:
        sym = structure.symmetries.OH
      elif self.ax == self.ay:
        sym = structure.symmetries.D4Hxy
      elif self.ax == self.az:
        sym = structure.symmetries.D4Hxz
      elif self.ay == self.az:
        sym = structure.symmetries.D4Hyz
      else:
        sym = structure.symmetries.D2H

    self.nsym     = sym.nsym
    self.symop    = sym.symop
    self.invsymop = sym.invsymop
    self.strsym   = sym.strsym

    logger.info('Detected {:3s} symmetry.'.format(self.strsym))


  def computeData(self, tbdata, charge, mu=None):
    self.tbdata  = tbdata
    self.charge  = charge

    self.energyBandMax  = int(np.max(tbdata[:,3]))
    self.tbparams       = tbdata.shape[0]
    self.opticalBandMax = self.energyBandMax

    if self.charge <= 0 or self.charge >= self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands (out of range (0,2*bands))')

    self._setupArrays()
    self._computeDispersion()
    self._calcFermiLevel(mu)

  def _setupKmesh(self):
    '''
    Setup the kmesh in the interval [0,1) 2pi/a
    '''

    logger.info('Setting up kmesh with {} reducible kpoints'.format(self.nkp))
    logger.info('nkx = {}'.format(self.nkx))
    logger.info('nky = {}'.format(self.nky))
    logger.info('nkz = {}'.format(self.nkz))

    self._kmeshx = np.linspace(0,1,self.nkx,endpoint=False)
    self._kmeshy = np.linspace(0,1,self.nky,endpoint=False)
    self._kmeshz = np.linspace(0,1,self.nkz,endpoint=False)

    if self.kshift:
      if self.nkx > 1: self._kmeshx += 1./self.nkx/2.
      if self.nky > 1: self._kmeshy += 1./self.nky/2.
      if self.nkz > 1: self._kmeshz += 1./self.nkz/2.

    mesh = np.meshgrid(self._kmeshx,self._kmeshy,self._kmeshz)
    kpoints = np.array(list(zip(*(dim.flat for dim in mesh))))

    if self.irreducible:
      logger.info('Generating irreducible kpoints:')
      kgen = set([]) # we start with the gamma point
      kirr = collections.OrderedDict()

      for ik in range(self.nkp):
        es.ElectronicStructure.progressBar(ik+1,self.nkp)

        kfloat = kpoints[ik,:]
        kstring = tuple(['{0:.7f}'.format(s%1) for s in kfloat])

        # if its already there .. skip completely
        if kstring in kgen:
          continue

        size_before = len(kgen) # size before we start adding

        knew = np.matmul(self.symop,kfloat) # generate all new k-vectors
        kmod = knew%1                  # shift them into the BZ
        for isym in range(self.nsym):       # iterate through them and add them to the set
          kshift = ['{0:.7f}'.format(s) for s in kmod[isym]]
          kgen.add(tuple(kshift))

        size_after = len(kgen) # size after we added

        kirr[kstring] =  size_after-size_before # difference = multiplicity

        if size_after == self.nkp: # if we arrived at all possible k-points break out of loop
          es.ElectronicStructure.progressBar(self.nkp,self.nkp)
          break

      self.kpoints = []
      self.multiplicity = []
      for kk, m in kirr.items():
        self.kpoints.append(np.array(kk, dtype=np.float64))
        self.multiplicity.append(m)

      self.kpoints = np.array(self.kpoints, dtype=np.float64)
      self.nkp = self.kpoints.shape[0]
      self.multiplicity = np.array(self.multiplicity, dtype=np.float64)
      self.weights      = self.weightsum * self.multiplicity / np.sum(self.multiplicity)
      logger.info('Generated irreducible kmesh with {} irreducible kpoints'.format(self.nkp))
    else:
      self.kpoints = kpoints
      self.nkp = self.kpoints.shape[0]
      self.weights = np.full((self.nkp,), fill_value=2./self.nkp, dtype=np.float64)
      self.multiplicity = np.ones((self.nkp,), dtype=np.float64)

  def _computeDispersion(self):
    '''
    create the energy dispersion from the following arrays:
    self.e0[nbands]
    self.hopping[nbands]

    e0 : 0 point energies of the bands
    hopping : hopping parameters t in the Hubbard Hamiltonian
              one minus sign already included
    '''

    nparas = self.tbdata

    # [0,1) -> [0,2pi)
    # dispersion
    ek = np.zeros((self.nkp, self.energyBandMax,), dtype=np.complex128)
    vk = np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.complex128)
    ck = np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.complex128)

    for itb in range(self.tbparams):
      rvec = self.tbdata[itb,:3]
      band = int(self.tbdata[itb,3]) - 1
      hop  = self.tbdata[itb,4]

      if np.sum(np.abs(rvec)) == 0:
        hop = -hop

      ''' e(k) = - t * e^{i k.r}'''
      ek[:,band] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1))

      ''' v(k) = de(k)/dk= - t * e^{i k.r} * i * r'''
      for i in range(3):
        vk[:,band, i] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1)) \
                             * 1j * rvec[i] * self.spacing[i]

      ''' c(k) = dv(k)/dk= - t * e^{i k.r} * (-1) * rr '''
      for i in range(3):
        ck[:,band, i] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1)) \
                             * (-1) * rvec[i]**2 * self.spacing[i]**2
      ck[:,band, 3] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1)) \
                           * (-1) * rvec[0] * rvec[1] * self.spacing[0] * self.spacing[1]
      ck[:,band, 4] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1)) \
                           * (-1) * rvec[0] * rvec[2] * self.spacing[0] * self.spacing[2]
      ck[:,band, 5] -= hop * np.exp(1j * np.sum(2*np.pi*self.kpoints * rvec[None,:], axis=1)) \
                           * (-1) * rvec[1] * rvec[2] * self.spacing[1] * self.spacing[2]


    if np.any(np.abs(ek.imag) > 1e-5):
      logger.warn('Detected complex energies ... check tight-binding parameter set')
    if np.any(np.abs(vk.imag) > 1e-5):
      logger.warn('Detected complex velocities ... check tight-binding parameter set')
    if np.any(np.abs(ck.imag) > 1e-5):
      logger.warn('Detected complex curvatures ... check tight-binding parameter set')

    self.energies[0][:,:]     = ek.real
    self.velocities[0][:,:,:] = vk.real
    self.curvatures[0][:,:,:] = ck.real

    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)

    if self.irreducible:
      '''
          unfortunately here we have code duplication
          to calculate the energies / velocities / curavtures of the generated k-points
          which we symmetrize over
      '''
      logger.info('Symmetrizing optical elements')

      for ikp in range(self.nkp):
        es.ElectronicStructure.progressBar(ikp+1,self.nkp)

        vel     = self.velocities[0][ikp,:,:]
        cur     = self.curvatures[0][ikp,:,:]

        # put the curvatures into matrix form
        curmat  = np.zeros((self.energyBandMax,3,3), dtype=np.float64)
        curmat[:, [0,1,2,1,2,2], [0,1,2,0,0,1]] = cur[:,:]
        curmat[:, [0,0,1], [1,2,2]] = curmat[:, [1,2,2], [0,0,1]]

        # generate the transformed velocities and curvatures
        vk = np.einsum('nij,bj->bni',self.symop,vel) # bands, nsym, 3 -- active transormation
        ck = np.einsum('nij,bjk,nkl->bnil',self.invsymop,curmat,self.symop) # bands, nsym, 3, 3

        # take the mean over the squares
        vk2 = vk[:,:,[0,1,2,0,0,1]] * vk[:,:,[0,1,2,1,1,2]]
        vk2 = np.mean(vk2,axis=1)

        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bnx,bnj,bnyi->bnxyz',levmatrix,vk,vk,ck)
        mb = np.mean(mb,axis=1)

        self.opticalMoments[0][ikp,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
                                          = vk2[:,:3] # only use xyz since we enforce orthogonality in this class
        self.BopticalDiag[0][ikp,:,:,:,:] = mb

    else: # reducible

      # Mbopt _abc = epsilon cij * va * vj * M^-1 bi # << BoltrzTrap 1 CPC !
      # doi: 10.1016/j.cpc.2006.03.007
      vk = self.velocities[0][:,:,:]

      ck = np.zeros((self.nkp,self.energyBandMax,3,3), dtype=np.float64)
      ck[:,:, [0,1,2,1,2,2], [0,1,2,0,0,1]] = self.curvatures[0][:,:,:]
      ck[:,:, [0,0,1], [1,2,2]] = ck[:,:, [1,2,2], [0,0,1]]

      mb = np.einsum('zij,pbx,pbj,pbyi->pbxyz',levmatrix,vk,vk,ck)

      self.BopticalDiag[0][...] = mb
      self.opticalMoments[0][:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
                  = self.velocities[0][:,:,[0,1,2]] \
                  * self.velocities[0][:,:,[0,1,2]]

    # last thing: take the diagonal elements
    self.opticalDiag[0][:,:,:] = self.opticalMoments[0][:,np.arange(self.energyBandMax), np.arange(self.energyBandMax), :]
