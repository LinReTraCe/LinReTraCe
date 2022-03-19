#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import logging
logger = logging.getLogger(__name__)

import numpy as np
import BoltzTraP2.dft as BTP
import BoltzTraP2.bandlib as BL
import BoltzTraP2.io as IO
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import serialization
from BoltzTraP2.misc import ffloat

try:
  import ase.spacegroup
  ase_exists = True
except ImportError:
  ase_exists = False

from structure import es
from structure import dft
from structure import units
from structure.aux import levicivita

class MetaW2kLoader(BTP.GenericWien2kLoader):
  '''
  BoltzTrap Custom Wien2k Loader.
  After setting the class variables one can register the Loader
  und use the provided custom files.
  The usual Wien2kLoader can only access the energy and energyso files.
  We also want to access up and dn files in spin-polarized calculations.
  '''

  # define class variables
  weightsum = None
  fscf      = None
  fstruct   = None
  fenergy   = None

  # access them here
  def __init__(self, directory):
      super(MetaW2kLoader, self).__init__(MetaW2kLoader.case, \
                                       MetaW2kLoader.weightsum, \
                                       MetaW2kLoader.fscf, \
                                       MetaW2kLoader.fstruct, \
                                       MetaW2kLoader.fenergy)

  @classmethod
  def setfiles(cls, case, weightsum, fscf, fstruct, fenergy):
    cls.case      = case
    cls.weightsum = weightsum
    cls.fscf      = fscf
    cls.fstruct   = fstruct
    cls.fenergy   = fenergy



class BTPInterpolation(object):
  '''
  BoltzTrap Interpolation class which we initialize with a DFTcalculation subclass object.
  The interpolation is implemented for w2k and vasp
  and automatically "registers" the boltztrap loaders with the files
  saved in the DFTcalculation object.
  The data is interpolated and the velocities and curvatures
  are evaluted on the original kmesh.
  '''

  def __init__(self, dftcalc):
    self.dftcalc = dftcalc

    self.spins         = self.dftcalc.spins
    self.energies      = []
    self.velocities    = []
    self.curvatures    = []
    self.opticalDiag   = []
    self.BopticalDiag  = []

  def interpolate(self, niter):
    logger.info('BoltzTrap2: Interpolating band-structure.')
    logger.info('BoltzTrap2: Requesting interpolation parameter: {}'.format(niter))
    self.niter = niter

    if isinstance(self.dftcalc, dft.w2kcalculation):
      for ispin in range(self.dftcalc.spins):

        # set class energy files
        MetaW2kLoader.setfiles(self.dftcalc.case, self.dftcalc.weightsum, self.dftcalc.fscf, \
                               self.dftcalc.fstruct, self.dftcalc.fenergyaccess[ispin])
        # and register the new loader
        BTP.register_loader("Custom", MetaW2kLoader)
        # interpolate and save
        self._interp()
        # self._rotate()
        if self.dftcalc.spinorbit:
          self._save1_so()
        else:
          self._save1()

        # otherwise we would loop a second time over the same file
        if self.dftcalc.spinorbit:
          break

    elif isinstance(self.dftcalc, dft.vaspcalculation):

      if self.dftcalc.irreducible and not ase_exists:
        raise IOError('vasp interpolation only works for reducible kmeshes or with ase installation')

      if self.dftcalc.irreducible:
        logger.info('To continue, the spacegroup of your vasp calculation is required.')
        inputmethod = input if sys.version_info >= (3, 0) else raw_input
        spacegroup = inputmethod('Spacegroup [1-230]: ')
        try:
          sg = int(spacegroup)
          asespacegroup = ase.spacegroup.Spacegroup(sg)
          self.dftcalc.symop = asespacegroup.get_rotations()
          self.dftcalc.invsymop = np.linalg.inv(self.dftcalc.symop)
          self.dftcalc.nsym = self.dftcalc.symop.shape[0]
        except:
          raise IOError('Something went wrong')

      BTP.register_loader("VASP", BTP.VASPLoader)
      self._interp()

      if (self.dftcalc.spins == 1): # this also works for spins==1 and non-collinear!
        self._save1()
      else:
        self._save2() # we have to split the interpolated arrays we get from BTP
    else:
      raise IOError('Received object which is not a supported calculation type')

    self._symmetrize()
    logger.info('BoltzTrap2: Interpolation successful.')


  def _interp(self):
    '''
    Standard BoltzTrap2 Library interface to interpolate
    the band energies, band velocities and band curvatures
    '''

    # we disable BTP logging if we are not in DEBUG mode!
    # disable = not (logger.getEffectiveLevel() == logging.DEBUG)
    disable = True

    if disable:
      logging.disable( sys.maxsize if sys.version_info >= (3,0) else sys.maxint)

    self.data = BTP.DFTData(self.dftcalc.directory, derivatives=False) # ignore mommat2 files
    self.equivalences = sphere.get_equivalences(self.data.atoms, self.data.magmom, \
                                                self.niter * len(self.data.kpoints))

    self.coeffs = fite.fitde3D(self.data, self.equivalences)
    self.metadata = serialization.gen_bt2_metadata(self.data, self.data.magmom is not None)

    self.lattvec = self.data.get_lattvec()

    self.interp_energies, self.interp_velocities, self.interp_curvatures = \
        fite.getBands(self.data.kpoints, self.equivalences, self.lattvec, self.coeffs, curvature=True)

    if disable:
      logging.disable(logging.NOTSET)

  # def _rotate(self):

  #   def isDiag(M):
  #     M[np.abs(M)<1e-10] = 0 # truncate numerical inaccuracy
  #     i, j = np.nonzero(M)
  #     return np.all(i==j)

  #   if isDiag(self.lattvec):
  #     logger.info('BoltzTrap2: Orthogonal lattice vectors.')
  #     return # do nothing

  #   logger.info('BoltzTrap2: Non-Orthogonal lattice vectors: Rotating elements')

  #   _, vecs = np.linalg.eig(self.lattvec)
  #   print(vecs)
  #   invvecs = np.linalg.inv(vecs)
  #   print(invvecs)
  #   # vecs[:,i] contains the ith eigenvalue

  #   bands, nkp = self.interp_velocities.shape[1:]
  #   for iband in range(bands):
  #     for ikp in range(nkp):
  #       vel = self.interp_velocities[:,iband,ikp]
  #       self.interp_velocities[:,iband,ikp] = vel @ vecs
  #       cur = self.interp_curvatures[:,:,iband,ikp]
  #       self.interp_curvatures[:,:,iband,ikp] = invvecs @ cur @ vecs

  def _save1(self):
    '''
    Saving the energies, velocities and curvatures by simply
    appending them to the lists we initialized at init time
    '''

    # we get the energies on the Hartree scale
    # rescaling to eV!
    self.energies.append(self.interp_energies.transpose(1,0) * units.hartree2eV)
    self.velocities.append(self.interp_velocities.transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)

    # here we dont save unnecessary information
    # because d2/dxdy = d2/dydx
    #
    # 1 4 5
    # - 2 6
    # - - 3

    # my index array to get exactly the values above on one axis
    d2ksave = tuple((np.array([0,1,2,1,2,2]), np.array([0,1,2,0,0,1])))

    # this works but I dont know whether there are better ways to do this
    # looks rather hacky
    tmp = self.interp_curvatures.transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
    tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)
    for i in range(tmp.shape[0]):
      for j in range(tmp.shape[1]):
        view = tmp[i,j,:,:]
        tmp2[i,j,:] = view[d2ksave]

    self.curvatures.append(tmp2)

  def _save1_so(self):
    '''
    Identical to _save1 only applied to cases where we have spin orbit coupling
    Wien2K saves the data in one file (energyso or energysoup) where the 'spins' alternate
    '''

    # we get the energies on the Hartree scale
    # rescaling to eV!
    self.interp_energies = self.interp_energies.transpose(1,0) * units.hartree2eV
    self.energies.append(self.interp_energies[:,::2])
    self.energies.append(self.interp_energies[:,1::2])

    self.interp_velocities = self.interp_velocities.transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom
    self.velocities.append(self.interp_velocities[:,::2,:])
    self.velocities.append(self.interp_velocities[:,1::2,:])

    # here we dont save unnecessary information
    # because d2/dxdy = d2/dydx
    #
    # 1 4 5
    # - 2 6
    # - - 3

    # my index array to get exactly the values above on one axis
    d2ksave = tuple((np.array([0,1,2,1,2,2]), np.array([0,1,2,0,0,1])))

    # this works but I dont know whether there are better ways to do this
    # looks rather hacky
    tmp = self.interp_curvatures.transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
    tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)
    for i in range(tmp.shape[0]):
      for j in range(tmp.shape[1]):
        view = tmp[i,j,:,:]
        tmp2[i,j,:] = view[d2ksave]

    self.curvatures.append(tmp2[:,::2])
    self.curvatures.append(tmp2[:,1::2])

  def _save2(self):
    '''
    In the case of a spin-dependent VASP calculation BTP2
    creates arrays which lists the data in order energyup energydn.
    We want them to be split, which is what this routines does
    '''

    nbands = self.interp_energies.shape[0] # this is guaranteed to be even here
    self.energies.append(self.interp_energies[:nbands//2,:].transpose(1,0) * units.hartree2eV)
    self.energies.append(self.interp_energies[nbands//2:,:].transpose(1,0) * units.hartree2eV)
    self.velocities.append(self.interp_velocities[:,:nbands//2,:].transpose(2,1,0) * units.hartree2eV * untis.bohr2angstrom)
    self.velocities.append(self.interp_velocities[:,nbands//2:,:].transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)
    # the last two elements dont matter here, since its symmetric anyways

    # my index array
    d2ksave = tuple((np.array([0,1,2,1,2,2]), np.array([0,1,2,0,0,1])))

    # some numpy magic to turn the 3x3 array into the only 6 necessary entries
    tmp = self.interp_curvatures[:,:,:nbands//2,:].transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
    tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)

    for i in range(tmp.shape[0]):
      for j in range(tmp.shape[1]):
        view = tmp[i,j,:,:]
        tmp2[i,j,:] = view[d2ksave]
    self.curvatures.append(tmp2)

    tmp = self.interp_curvatures[:,:,nbands//2:,:].transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
    for i in range(tmp.shape[0]):
      for j in range(tmp.shape[1]):
        view = tmp[i,j,:,:]
        tmp2[i,j,:] = view[d2ksave]
    self.curvatures.append(tmp2)

  def _symmetrize(self):
    '''
    If we want to be able to use these elements on the irreducible k-grid
    we need to symmetrized them [akin to what wien2k does]
    '''

    nsym = self.dftcalc.nsym
    syms = self.dftcalc.symop
    symsinv = self.dftcalc.invsymop

    # symmetrize each vector and matrix

    logger.info('BoltzTrap2: Symmetrizing band derivatives.')

    d2ksave = tuple((np.array([0,1,2,1,2,2]), np.array([0,1,2,0,0,1])))

    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)

    for ispin in range(self.dftcalc.spins):

      if self.dftcalc.spins == 1:
        prefix = ''
      else:
        if ispin == 0:
          prefix = 'up:'
        else:
          prefix = 'dn:'

      nkp, nbands = self.velocities[ispin].shape[:2]

      BopticalDiag = np.zeros((nkp,nbands,3,3,3), dtype=np.float64)
      opticalDiag  = np.zeros((nkp,nbands,6), dtype=np.float64)

      for ikp in range(nkp):
        es.ElectronicStructure.progressBar(ikp+1,nkp, status='k-points', prefix=prefix)


        vel     = self.velocities[ispin][ikp,:,:] # nbands, 3
        cur     = self.curvatures[ispin][ikp,:,:] # nbands, 6

        # put the curvatures in symmetric matrix form
        curmat  = np.zeros((nbands,3,3), dtype=np.float64)
        curmat[:, [0,1,2,1,2,2], [0,1,2,0,0,1]] = cur[:,:]
        curmat[:, [0,0,1], [1,2,2]] = curmat[:, [1,2,2], [0,0,1]]

        # the k-points transform like k' = R k
        # therefore the velocities transform identical: v' = R v
        # the equivalent matrix transformation is then R^{-1} c R

        vk = np.einsum('nij,bj->bni',syms,vel) # bands, nsym, 3
        ck = np.einsum('nij,bjk,nkl->bnil',symsinv,curmat,syms) # bands, nsym, 3, 3

        vk2 = vk[:,:,[0,1,2,0,0,1]] * vk[:,:,[0,1,2,1,1,2]]
        vk2 = np.mean(vk2,axis=1) # symmetrize over the squares

        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bnx,bnj,bnyi->bnxyz',levmatrix,vk,vk,ck)
        mb = np.mean(mb,axis=1)

        opticalDiag[ikp,:,:] = vk2
        BopticalDiag[ikp,:,:,:,:] = mb

      self.opticalDiag.append(opticalDiag)
      self.BopticalDiag.append(BopticalDiag)

    # if we need the peierls approximation
    self.opticalBandMin = 0
    self.opticalBandMax = self.velocities[0].shape[1]
