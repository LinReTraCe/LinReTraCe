#!/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np
try:
  import spglib
  spglib_exist = True
except ImportError:
  spglib_exist = False

from   structure.model import Model

from   structure.aux import levicivita
from   structure.aux import progressBar

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

class OrthogonalTightBinding(Model):
  '''
  Tightbinding class to create (ir-)reducible k-grids
  for orthogonal lattices in 1 / 2 / 3 dimensions
  Main class that gets interfaces with lotb script
  This class can be used for irreducible grids without spglib
  '''

  def __init__(self, ax=1, ay=1, az=1, nkx=1, nky=1, nkz=1, irreducible=True, kshift=False):
    super(OrthogonalTightBinding, self).__init__(nkx,nky,nkz)
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
      self._checkSymmetriesKmesh()
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

    logger.debug('Spglib symmetries:\n{}'.format(self.symop))

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
        knew = np.einsum('nji,j->ni',self.symop,kpoints[ik,:])
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

    bandmin = int(np.min(self.tbdata[:,3]))
    bandmax = int(np.max(self.tbdata[:,3]))

    if bandmin != 1:
      raise IOError('Error: tight binding parameter set must start at band 1')

    bandcheck = np.full((bandmax,), fill_value=False)
    for itb in range(self.tbdata.shape[0]):
      band = int(self.tbdata[itb,3]) - 1 # band identifier
      bandcheck[band] = True

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

  def _checkSymmetriesKmesh(self):
    '''
      compare the momentum mesh to the unit cell symmetries
      raise an IOError if there is disagreement
    '''
    nkvec = np.array([self.nkx,self.nky,self.nkz], dtype=np.float64)
    transformed = np.abs(np.einsum('nij,j->ni',self.invsymop,nkvec))
    equal_mesh = np.isclose(nkvec,transformed)
    equal_zero = np.isclose(np.zeros((3,)),transformed)
    equal  = np.all(np.logical_or(equal_zero,equal_mesh))
    logger.info('Momentum grid symmetry check: {}'.format(str(equal)))

    if not equal:
      raise IOError('\nMomentum mesh does not respect point group symmetries of unit cell.' + \
                    '\nIrreducible calculations must have momentum grid that does.')



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
        vk = np.einsum('nij,bj->bni',self.invsymop,vel) # bands, nsym, 3 -- active transormation
        ck = np.einsum('nij,bjk,nlk->bnil',self.invsymop,curmat,self.invsymop) # bands, nsym, 3, 3

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
