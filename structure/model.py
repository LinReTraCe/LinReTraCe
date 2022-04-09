#!/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import logging
import collections
logger = logging.getLogger(__name__)

import numpy as np
try:
  import spglib
  spglib_exist = True
except ImportError:
  spglib_exist = False

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
from   structure.aux import progressBar

# python 2 & 3 compatible usage of abstract base classes
if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

from structure.es import ElectronicStructure

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

    self.nsym           = structure.symmetries.C1.nsym
    self.symop          = structure.symmetries.C1.symop
    self.invsymop       = structure.symmetries.C1.invsymop

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


class orthogonal_tightbinding(Model):
  '''
  Tightbinding class to create (ir-)reducible k-grids
  for orthogonal lattices in 1 / 2 / 3 dimensions
  Main class that gets interfaces with lotb script
  This class can be used for irreducible grids without spglib
  '''

  def __init__(self, ax=1, ay=1, az=1, nkx=1, nky=1, nkz=1, irreducible=True, kshift=False):
    super(orthogonal_tightbinding, self).__init__(nkx,nky,nkz)
    self.ax          = ax
    self.ay          = ay
    self.az          = az
    self.spacing     = [self.ax, self.ay, self.az]
    self.irreducible = irreducible  # generate irreducible grid instead of reducible
    self.kshift      = kshift       # shift by half a k-point to avoid Gamma point
    self.vol         = self.ax*self.ay*self.az
    self.ortho       = True

    logger.info('Setting up kmesh with {} reducible kpoints'.format(self.nkp))
    logger.info('nkx = {}'.format(self.nkx))
    logger.info('nky = {}'.format(self.nky))
    logger.info('nkz = {}'.format(self.nkz))
    if self.irreducible:
      try:
        self._setupKmeshSpglib()
      except Exception as e:
        self._defineSymmetries()
        self._setupKmesh()
    else:
      self._setupKmesh()

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
    self._checkFile()

    self.energyBandMax  = int(np.max(tbdata[:,3]))
    self.tbparams       = tbdata.shape[0]
    self.opticalBandMax = self.energyBandMax

    if self.charge <= 0 or self.charge >= self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands (out of range (0,2*bands))')

    self._setupArrays(self.ortho) # this takes an explicit argument
    self._computeDispersion()
    self._calcFermiLevel(mu)

  def _setupKmeshSpglib(self):
    '''
        employ spglib to create irreducible kpoints,
        multiplicity and symmetry operations
    '''
    import spglib # we are protected by a try except from outside

    ''' define k-grid, shift if required '''
    kgrid = np.array([self.nkx,self.nky,self.nkz], dtype=np.int)

    if self.kshift:
      is_shift = np.array([int(i) for i in self.dims], dtype=np.float64)
    else:
      is_shift = np.array([0,0,0], dtype=np.int)

    # primitive tight-binding lattice model
    lattice = np.array([[self.ax,0.0,0.0],\
                        [0.0,self.ay,0.0],\
                        [0.0,0.0,self.az]])
    positions = [[0.0,0.0,0.0]] # only one atomic position at origin
    numbers = [1,]              # to distinguish atomic species

    logger.info('Spglib: Generating irreducible kpoints.')

    cell = (lattice, positions, numbers)
    mapping, grid = spglib.get_ir_reciprocal_mesh(kgrid, cell, is_shift=is_shift)

    unique, counts = np.unique(mapping, return_counts=True)
    self.nkp  = len(unique)

    logger.info('Spglib: Generated irreducible kmesh with {} irreducible kpoints'.format(self.nkp))

    ''' from the mapping and counts thereof generate multiplicity and weights '''
    self.multiplicity = np.array(counts, dtype=int)
    self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))

    self.kpoints = grid[unique]
    self.kpoints = (self.kpoints + is_shift.astype(np.float64)/2.) / kgrid[None,:].astype(np.float64)


    ''' get symmetry and reduce unnecessary ones '''
    symmetry = spglib.get_symmetry(cell, symprec=1e-5)
    symsfull = symmetry['rotations']

    self.symop = []
    for ii in np.arange(symsfull.shape[0]):
      isym = symsfull[ii]
      to_add = True
      for i, dim in enumerate(self.dims):
        if dim: continue
        for j in range(3):
          if i==j:
            if isym[i,j] != 1: to_add = False
          else:
            if isym[j,i] != 0: to_add = False
            if isym[i,j] != 0: to_add = False # redundant I think

      if to_add: self.symop.append(isym)

    self.symop = np.array(self.symop)
    self.invsymop = np.linalg.inv(self.symop)
    self.nsym = self.symop.shape[0]

  def _setupKmesh(self):
    '''
    Setup the kmesh in the interval [0,1) 2pi/a
    '''

    self._kmeshx = np.linspace(0,1,self.nkx,endpoint=False)
    self._kmeshy = np.linspace(0,1,self.nky,endpoint=False)
    self._kmeshz = np.linspace(0,1,self.nkz,endpoint=False)

    if self.kshift:
      self._kmeshshift = []
      for ik in [self.nkx,self.nky,self.nkz]:
        if ik > 1:
          self._kmeshshift.append(1./ik/2.)
        else:
          self._kmeshshift.append(0.0)
      self._kmeshshift = np.array(self._kmeshshift, dtype=np.float64)

    # the way these points are ordered is important for the indexing below
    kpoints = []
    for ikx in self._kmeshx:
      for iky in self._kmeshy:
        for ikz in self._kmeshz:
          kpoints.append([ikx,iky,ikz])
    kpoints = np.array(kpoints)
    if self.kshift: kpoints += self._kmeshshift[None,:]

    unique  = np.ones((self.nkx*self.nky*self.nkz), dtype=np.int)
    mult    = np.zeros((self.nkx*self.nky*self.nkz), dtype=np.float64)
    irrk    = 0

    if self.irreducible:
      logger.info('Generating irreducible kpoints:')

      for ik in range(self.nkp):
        progressBar(ik+1,self.nkp,status='k-points')

        if unique[ik] == 0: continue # skip if we already went there via symmetry
        irrk += 1    # new point -> increase irreducible counter
        mult[ik] = 1 # reset multiplicity counter

        ''' generate all the symmetry related k-points in the Brillouin zone '''
        knew = np.einsum('nij,j->ni',self.symop,kpoints[ik,:])
        kmod = knew%1
        ''' in order to index properly and if kshift is applied , shift back '''
        if self.kshift:
          kmod -= self._kmeshshift
        ''' round to neareast integer '''
        kround = np.rint(kmod * np.array([self.nkx,self.nky,self.nkz])[None,:])
        ''' exact floating calculation '''
        kexact = kmod * np.array([self.nkx,self.nky,self.nkz])[None,:]
        ''' only use the values that transform properly on all three axes '''
        mask = np.all(np.isclose(kround,kexact),axis=1)
        ''' apply the mask to filter '''
        kmask = kround[mask]
        ''' get the hash index '''
        kindex = (kmask[:,2] + \
                  kmask[:,1] * self.nkz + \
                  kmask[:,0] * self.nkz * self.nky).astype(int)
        ''' remove the k-points connected via symmetry and increase the multiplicity accordingly '''
        for ikk in kindex:
          if ikk <= ik: continue
          if unique[ikk]:
            unique[ikk] = 0
            mult[ik] += 1

      self.nkp = irrk
      self.kpoints = kpoints[unique>0]
      self.multiplicity = mult[unique>0]
      self.weights      = self.weightsum * self.multiplicity / np.sum(self.multiplicity)
      logger.info('Generated irreducible kmesh with {} irreducible kpoints'.format(self.nkp))
    else:
      self.nkp = self.nkx * self.nky * self.nkz
      self.kpoints = []
      for ikx in np.linspace(0,1,self.nkx,endpoint=False):
        for iky in np.linspace(0,1,self.nky,endpoint=False):
          for ikz in np.linspace(0,1,self.nkz,endpoint=False):
            self.kpoints.append([ikx,iky,ikz])
      self.kpoints = np.array(self.kpoints, dtype=np.float64)
      if self.kshift:
        self.kpoints += is_shift.astype(np.float64)/2. / kgrid[None,:].astype(np.float64)
      self.weights = np.full((self.nkp,), fill_value=2./self.nkp, dtype=np.float64)
      self.multiplicity = np.ones((self.nkp,), dtype=int)

  def _checkFile(self):
    '''
        Check if bands start at 1
        Check if every band has some given value
        Raise IOError if something is wrong
        Check symmetries ?
    '''

    bandmin = int(np.min(self.tbdata[:,3:5]))
    bandmax = int(np.max(self.tbdata[:,3:5]))

    if bandmin != 1:
      raise IOError('Error: tight binding parameter set must start at band 1')

    bandcheck = np.full((bandmax,), fill_value=False)
    for itb in range(self.tbdata.shape[0]):
      band1 = int(self.tbdata[itb,3]) - 1 # band identifier
      band2 = int(self.tbdata[itb,4]) - 1 # band identifier
      bandcheck[band1] = True
      bandcheck[band2] = True

    if not np.all(bandcheck):
      raise IOError('Error: tight binding parameter set does not contain all bands')

    for itb1 in range(self.tbdata.shape[0]):
      band1  = int(self.tbdata[itb1,3]) - 1 # band identifier
      rvec1  = self.tbdata[itb1,:3].astype(int)
      hop1   = self.tbdata[itb1,4]

      rvecsym = np.einsum('nij,j->ni',self.symop,rvec1)

      for isym in range(self.nsym):
        transformed = False
        rvec_transformed = rvecsym[isym]

        for itb2 in range(self.tbdata.shape[0]):
          band2  = int(self.tbdata[itb1,3]) - 1 # band identifier
          if band1 != band2: continue

          rvec2  = self.tbdata[itb2,:3].astype(int)
          hop2   = self.tbdata[itb2,4]
          if np.allclose(rvec_transformed,rvec2) and np.abs(hop1-hop2) < 1e-6:
            transformed = True
            break

        if not transformed:
          logger.warning('Tight binding parameter set does not fulfill symmteries given by unit cell' + \
                        '\n symmetry of r-vector {} is not fulfilled'.format(rvec1) + \
                        '\n avoid irreducible calculation if this is done on purpose\n\n')
          return



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

      rvec = self.tbdata[itb,:3]         # r - vector
      band = int(self.tbdata[itb,3]) - 1 # band identifier
      hop  = self.tbdata[itb,4]          # hopping parameter

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
    self.velocities.append(vk.real)
    self.curvatures.append(ck.real)

    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)

    if self.irreducible:
      '''
          Symmetrize the velocity squares -> Mopt
          and the magnetic optical elements ~ v*v*c -> Mbop
      '''
      logger.info('Symmetrizing optical elements')

      for ikp in range(self.nkp):
        progressBar(ikp+1,self.nkp,status='k-points')

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
        vk2 = vk[:,:,[0,1,2]] * vk[:,:,[0,1,2]]
        vk2 = np.mean(vk2,axis=1)

        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bnx,bnj,bnyi->bnxyz',levmatrix,vk,vk,ck)
        mb = np.mean(mb,axis=1)

        self.opticalMoments[0][ikp,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:] \
                                          = vk2[:,:] # only use xyz since we enforce orthogonality in this class
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



class tightbinding(Model):
  '''
  Tightbinding class to create (ir-)reducible k-grids
  for arbitrary lattices.

  Main class that gets interfaces with ltb script
  This class generates irreducible grids only via spglib.
  Reducible grids do not require spglib.
  '''

  def __init__(self, nkx=1, nky=1, nkz=1, irreducible=True, kshift=False):
    super(tightbinding, self).__init__(nkx,nky,nkz)
    self.irreducible = irreducible  # generate irreducible grid instead of reducible
    self.kshift      = kshift       # shift by half a k-point to avoid Gamma point

    logger.info('Setting up tight binding with {} x {} x {} kpoints'.format(self.nkx,self.nky,self.nkz))
    if self.irreducible and not spglib_exist:
      raise IOError('Irreducible grids are only supported via spglib')


  def computeData(self, tbdata, rvecdata, atomdata, charge, mu=None):
    self.tbdata   = tbdata
    self.rvecdata = rvecdata
    self.atomdata = atomdata
    self.charge   = charge

    if self.tbdata.shape[1] == 6:
      self.imaghopping = False
    elif self.tbdata.shape[1] == 7:
      self.imaghopping = True

    self._computeOrthogonality() # sets self.ortho
    self._computeReciprocLattice() # defines reciprocal lattice vectors and volume
    self._readHr() # reads hr in formatted array, detects unique r-vectors, detects inter-band transitions

    if self.charge < 0 or self.charge > self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands : charge in [0,2*bands]')

    self._setupKmesh()
    if self.irreducible:
      self._checkSymmetries()

    ''' after getting the correct number of k-points we can setup the arrays '''
    self._setupArrays(self.ortho)

    self._computeHk()
    self._calcOptical()

    ''' setting some important variables '''
    self.opticalBandMin = 0
    self.opticalBandMax = self.energyBandMax

    self._calcFermiLevel(mu)

  def _computeOrthogonality(self):
    sum_vecs = np.sum(self.rvecdata, axis=1) # a_1 + a_2 + a_3
    max_vecs = np.array([np.max(np.abs(self.rvecdata[:,i])) for i in range(3)]) # maximal entries
    ratio = sum_vecs / max_vecs
    self.ortho = np.all(np.isclose(ratio, ratio[0]))
    logger.info('   orthogonal lattice: {}'.format(self.ortho))

  def _computeReciprocLattice(self):
    self.vol = np.dot(np.cross(self.rvecdata[0,:],self.rvecdata[1,:]),self.rvecdata[2,:])
    self.kvec = np.zeros_like(self.rvecdata, dtype=np.float64)
    self.kvec[:,0] = np.cross(self.rvecdata[1,:],self.rvecdata[2,:]) / self.vol
    self.kvec[:,1] = np.cross(self.rvecdata[2,:],self.rvecdata[0,:]) / self.vol
    self.kvec[:,2] = np.cross(self.rvecdata[0,:],self.rvecdata[1,:]) / self.vol
    self.kvec *= 2*np.pi
    logger.info('   volume [Ang^3]: {}'.format(self.vol))

  def _readHr(self):
    '''
      Check that all bands exist
      Read Hr data into array  -> self.hr
      Read all rvectors into array -> self.rpoints
    '''

    bandmin = int(np.min(self.tbdata[:,3:5]))
    bandmax = int(np.max(self.tbdata[:,3:5]))

    if bandmin != 1:
      raise IOError('Error: tight binding parameter set must start at band 1')

    inter = False
    bandcheck = np.full((bandmax,), fill_value=False)
    for itb in range(self.tbdata.shape[0]):
      band1 = int(self.tbdata[itb,3]) - 1 # band identifier
      bandcheck[band1] = True
      band2 = int(self.tbdata[itb,4]) - 1 # band identifier
      bandcheck[band2] = True
      if band1 != band2:
        inter = True

    if inter:
      self.opticfull  = True
      self.bopticfull = True
      logger.info('   Detected inter-band hopping parameters.')
    else:
      logger.info('   Detected only intra-band hopping parameters.')

    if not np.all(bandcheck):
      raise IOError('Error: tight binding parameter set does not contain all bands')

    self.energyBandMax  = bandmax
    logger.info('   Number of bands: {}'.format(self.energyBandMax))

    hrset = set()
    for i in range(self.tbdata.shape[0]):
      rx, ry, rz = self.tbdata[i,:3]
      hrset.add((rx,ry,rz))

    self.nrp = len(hrset)
    logger.info('   Number of unique r-points: {}'.format(self.nrp))
    self.rpoints = np.array(list(hrset), dtype=np.int)

    self.hr = np.zeros((self.nrp,self.energyBandMax,self.energyBandMax), dtype=np.complex128)
    for i in range(self.tbdata.shape[0]):
      rvec = self.tbdata[i,:3].astype(int)
      orb1, orb2 = self.tbdata[i,3:5].astype(int)
      hop = self.tbdata[i,5]
      if self.imaghopping:
        hop += 1j * self.tbdata[i,6]

      ir = np.argwhere(np.sum(np.abs(self.rpoints-rvec[None,:]),axis=1)==0)[0][0]
      if np.all(rvec==np.zeros((3,), dtype=np.int)) and orb1==orb2:
        hop *= (-1)
      self.hr[ir,orb1-1,orb2-1] += -hop

  def _setupKmesh(self):
    '''
      employ spglib to create irreducible kpoints,
      multiplicity and symmetry operations
      otherwise create the usual reducible cell
    '''

    kgrid = np.array([self.nkx,self.nky,self.nkz], dtype=np.int)
    if self.kshift:
      is_shift = np.array([int(i) for i in self.dims], dtype=np.float64)
    else:
      is_shift = np.array([0,0,0], dtype=np.int)

    ''' define k-grid, shift if required '''
    if self.irreducible:

      lattice   = self.rvecdata
      positions = []
      numbers   = []
      for i in range(self.atomdata.shape[0]):
        numbers.append(self.atomdata[i,0])
        positions.append(self.atomdata[i,1:])

      logger.info('Spglib: Generating irreducible kpoints.')

      cell = (lattice, positions, numbers)
      mapping, grid = spglib.get_ir_reciprocal_mesh(kgrid, cell, is_shift=is_shift)

      unique, counts = np.unique(mapping, return_counts=True)
      self.nkp  = len(unique)

      logger.info('Spglib: Generated irreducible kmesh with {} irreducible kpoints'.format(self.nkp))

      ''' from the mapping and counts thereof generate multiplicity and weights '''
      self.multiplicity = np.array(counts, dtype=int)
      self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))

      self.kpoints = grid[unique]
      self.kpoints = (self.kpoints + is_shift.astype(np.float64)/2.) / kgrid[None,:].astype(np.float64)

      ''' get symmetry and reduce unnecessary ones '''
      symmetry = spglib.get_symmetry(cell, symprec=1e-5)
      symsfull = symmetry['rotations']

      self.symop = []
      for ii in np.arange(symsfull.shape[0]):
        isym = symsfull[ii]
        to_add = True
        for i, dim in enumerate(self.dims):
          if dim: continue
          for j in range(3):
            if i==j:
              if isym[i,j] != 1: to_add = False
            else:
              if isym[j,i] != 0: to_add = False
              if isym[i,j] != 0: to_add = False # redundant I think

        if to_add: self.symop.append(isym)

      self.symop = np.array(self.symop)
      self.invsymop = np.linalg.inv(self.symop)
      self.nsym = self.symop.shape[0]
      logger.info('Spglib: Found {} symmetry operations'.format(self.nsym))
    else:
      self.nkp = self.nkx * self.nky * self.nkz
      self.kpoints = []
      for ikx in np.linspace(0,1,self.nkx,endpoint=False):
        for iky in np.linspace(0,1,self.nky,endpoint=False):
          for ikz in np.linspace(0,1,self.nkz,endpoint=False):
            self.kpoints.append([ikx,iky,ikz])
      self.kpoints = np.array(self.kpoints, dtype=np.float64)
      if self.kshift:
        self.kpoints += is_shift.astype(np.float64)/2. / kgrid[None,:].astype(np.float64)

      self.multiplicity = np.ones((self.nkp), dtype=int)
      self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))

      self.symop = np.array(self.symop)
      self.invsymop = np.linalg.inv(self.symop)
      self.nsym = self.symop.shape[0]
      logger.info('Reducible grid: {} symmetry operation'.format(self.nsym))

  def _computeHk(self):
    '''
      calculate H(k) d/dk H(k) and d2/dk2 H(k)
      via H(r) and primitive lattice vectors
    '''

    # hamiltonian H(k)
    hk = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax), dtype=np.complex128)
    # hamiltonian derivative d_dk H(k)
    # 3: x y z
    hvk = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3), dtype=np.complex128)
    # hamiltonian double derivative d2_dk2 H(k)
    # 6: xx yy zz xy xz yz
    hck = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,6), dtype=np.complex128)

    ''' FOURIERTRANSFORM hk = sum_r e^{i r.k} * weight(r) * h(r) '''

    rdotk = 2*np.pi*np.einsum('ki,ri->kr',self.kpoints,self.rpoints)
    ee = np.exp(1j * rdotk)
    hk[:,:,:] = np.einsum('kr,rij->kij', ee, self.hr)

    ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
    prefactor_r = np.einsum('di,ri->dr', self.rvecdata, self.rpoints)
    hvk[:,:,:,:] = np.einsum('dr,kr,rij->kijd',1j*prefactor_r,ee,self.hr)

    ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
    prefactor_r2 = np.zeros((6,self.nrp), dtype=np.float64)
    for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
      prefactor_r2[idir,:] = prefactor_r[i,:] * prefactor_r[j,:]
    hck[:,:,:,:] = np.einsum('dr,kr,rij->kijd',-prefactor_r2,ee,self.hr)

    ''' this transforms all k points at once '''
    ek, U = np.linalg.eig(hk)

    ''' Sort eigenvalues from smallest to largest
        Required for detection of possible gaps '''
    for ik in range(self.nkp):
      ekk, Uk = ek[ik,:], U[ik,:,:]
      idx = ekk.argsort()
      ek[ik,:] = ekk[idx]
      U[ik,:,:] = Uk[:,idx]

    ''' the velocities and curvatures are ordered according to e(k)
        due to the reordering of U '''
    Uinv = np.linalg.inv(U)
    vk = np.einsum('kab,kbci,kcd->kadi',Uinv,hvk,U)
    ck = np.einsum('kab,kbci,kcd->kadi',Uinv,hck,U)

    if np.any(np.abs(ek.imag) > 1e-5):
      print('Detected complex energies ... truncating')
    ''' apparently complex velocities and curvatures are allowed now '''

    self.energies[0][...] = ek.real
    self.velocities.append(vk)
    self.curvatures.append(ck)

  def _checkSymmetries(self):

    for itb1 in range(self.tbdata.shape[0]):
      rvec1  = self.tbdata[itb1,:3].astype(int)
      band1_1  = int(self.tbdata[itb1,3]) - 1 # band identifier
      band1_2  = int(self.tbdata[itb1,4]) - 1 # band identifier
      hop1   = self.tbdata[itb1,5]
      if self.imaghopping:
        hop1 += 1j * self.tbdata[itb1,6]

      rvecsym = np.einsum('nij,j->ni',self.symop,rvec1)

      for isym in range(self.nsym):
        transformed = False
        rvec_transformed = rvecsym[isym]

        for itb2 in range(self.tbdata.shape[0]):
          band2_1  = int(self.tbdata[itb1,3]) - 1 # band identifier
          band2_2  = int(self.tbdata[itb1,4]) - 1 # band identifier
          if band1_1 != band2_1 or band2_2 != band2_2: continue

          rvec2  = self.tbdata[itb2,:3].astype(int)
          hop2   = self.tbdata[itb2,5]
          if self.imaghopping:
            hop2 += 1j * self.tbdata[itb2,6]

          if np.allclose(rvec_transformed,rvec2) and np.abs(hop1-hop2) < 1e-6:
            transformed = True
            break

        if not transformed:
          logger.warning('\n\nTight binding parameter set does not fulfill symmteries given by unit cell' + \
                        '\n symmetry of r-vector {} does not respect full point group symmetries'.format(rvec1) + \
                        '\n break unit cell symmetry or avoid irreducible calculation if this is done on purpose\n\n')
          return
    else:
      logger.info('Tight binding parameter set fulfills full point group symmetry set determined by unit cell.')


  def _calcOptical(self):
    '''
      from the velocities and curvatures : calculate
      the optical elements and b-field optical elements
    '''

    logger.info('Calculating optical elements from derivatives')
    # devine levicity matrix so we can use it in einsteinsummations
    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)

    # symmetrizing is essential
    if self.irreducible:
      if self.ortho:
        loc_opticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3), dtype=np.float64)
      else:
        loc_opticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,9), dtype=np.float64)
      loc_BopticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3,3), dtype=np.complex128)
      for ikp in range(self.nkp):
        progressBar(ikp+1,self.nkp,status='k-points')

        vel     = self.velocities[0][ikp,:,:,:] # bands, bands, 3
        cur     = self.curvatures[0][ikp,:,:,:] # bands, bands, 6

        # put the curvatures into matrix form
        curmat  = np.zeros((self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
        curmat[:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:]
        curmat[:,:, [1,2,2], [0,0,0]] = curmat[:,:, [0,0,1], [1,2,2]]

        # generate the transformed velocities and curvatures
        vk = np.einsum('nij,bpj->bpni',self.symop,vel) # bands, bands, nsym, 3
        ck = np.einsum('nij,bpjk,nkl->bpnil',self.invsymop,curmat,self.symop) # bands, bands, nsym, 3, 3

        # take the mean over the squares
        vk2 = np.conjugate(vk[:,:,:,[0,1,2,0,0,1]]) * vk[:,:,:,[0,1,2,1,2,2]]
        vk2 = np.mean(vk2,axis=2)

        if self.ortho:
          loc_opticalMoments[ikp,...] = vk2[...,:3].real
        else:
          loc_opticalMoments[ikp,:,:,:6] = vk2.real
          loc_opticalMoments[ikp,:,:,6:] = vk2[...,3:].imag

        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bpnx,bpnj,bpnyi->bpnxyz',levmatrix,vk,vk,ck)
        mb = np.mean(mb,axis=2)
        loc_BopticalMoments[ikp,...] = mb


      self.opticalMoments[0][...]  = loc_opticalMoments
      self.opticalDiag[0][...]     = loc_opticalMoments[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:]
      self.BopticalMoments[0][...] = loc_BopticalMoments
      self.BopticalDiag[0][...]    = loc_BopticalMoments[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),...]

    else:
      vel = self.velocities[0]
      cur = self.curvatures[0]

      # transform into matrix form
      curmat  = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
      curmat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:,:]
      curmat[:,:,:, [1,2,2], [0,0,1]] = np.conjugate(curmat[:,:,:, [0,0,1], [1,2,2]])
      vel2 = np.conjugate(vel[:,:,:,[0,1,2,0,0,1]]) * vel[:,:,:,[0,1,2,1,1,2]]
      if self.ortho:
        vel2 = vel2[:,:,:,:3].real
      else:
        if np.any(np.abs(vel2.imag) > 1e-6):
          temp = vel2.copy()
          vel2 = np.empty((self.nkp,self.energyBandMax,self.energyBandMax,9), dtype=np.float64)
          vel2[:,:,:,:6] = temp
          vel2[:,:,:,6:] = temp[:,:,:,:3].imag
        else:
          vel2 = vel2.real
      self.opticalMoments = [vel2]
      vel2diag            = vel2[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:]
      self.opticalDiag    = [vel2diag]

        #           epsilon_cij v_a v_j c_bi -> abc
      mb = np.einsum('cij,knma,knmj,knmbi->knmabc',levmatrix,vel,vel,curmat)
      self.BopticalMoments[0][...] = mb
      mbdiag                       = mb[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:,:,:]
      self.BopticalDiag[0][...]    = mbdiag

    if not self.ortho:
      truncate = True
      if np.any(np.abs(self.opticalMoments[0][...,6:]) > 1e-6):
        truncate = False
      if truncate:
        self.opticalMoments[0]  = self.opticalMoments[0][...,:6]
        self.opticalDiag[0]     = self.opticalDiag[0][...,:6]
