#! /usr/bin/env python

#    Parts of this python file have been adopted directly from the source code
#    of BoltzTraP2. More specifically, the DFTData class, used
#    to load in DFT data, has been modified to be initializable
#    via our internal data.
#    For Licensing purposes, the following Header is required:

#    BoltzTraP2, a program for interpolating band structures and calculating
#                semi-classical transport coefficients.
#    Copyright (C) 2017-2024 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2024 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2024 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
#
#    This file is part of BoltzTraP2.
#
#    BoltzTraP2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BoltzTraP2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BoltzTraP2.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    Copyright (C) 2024 Matthias Pickem <matthias.pickem@gmail.com>
#    Adapted to be used directly for interpolation purposes
#    in the LinReTraCe code available at github.com/linretrace
#

from __future__ import print_function, division, absolute_import
import sys
import os
import logging
import math
logger = logging.getLogger(__name__)

import numpy as np

import BoltzTraP2.dft as BTP
import BoltzTraP2.bandlib as BL
import BoltzTraP2.io as IO
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import serialization
from BoltzTraP2.misc import ffloat
from BoltzTraP2.units import Angstrom
import ase.spacegroup

from structure.auxiliary import progressBar
from structure.wien2k    import Wien2kCalculation
from structure.vasp      import VaspCalculation
from structure           import units
from structure.auxiliary import levicivita



class BoltztrapInterpolation(object):
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
    self.bopticdiag    = True
    self.bopticfull    = False # we only have the intra band data

  def interpolate(self, niter = 3, mesh = None):
    logger.info('BoltzTrap2 - Licensed under GPLv3.')
    logger.info('BoltzTrap2: Interpolating band-structure.')
    logger.info('BoltzTrap2: Requesting interpolation parameter: {}'.format(niter))
    self.niter = niter
    if mesh is not None:
      self.mesh = [int(i) for i in mesh]
    else:
      self.mesh = None

    '''
      we require the spacegroup operations
      if the DFT Calculation does not provide them, we have to ask for space group
      from which they can be generated via the ase module

      the interactive interface is there as a fail safe.
    '''

    if isinstance(self.dftcalc, VaspCalculation) and self.dftcalc.irreducible:
      logger.info('\n\nASE detected spacegroup number: {}'.format(self.dftcalc.spacegroup))
      logger.info('If this is correct, skip by pressing enter.')
      logger.info('Otherwise, enter the new space group in the range 1-230.')
      inputmethod = input if sys.version_info >= (3, 0) else raw_input
      spacegroup = inputmethod('Spacegroup [1-230]: ')
      try:
        if len(spacegroup.strip()) == 0:
          sg = self.dftcalc.spacegroup
        else:
          sg = int(spacegroup)
      except Exception as e:
        raise IOError('Input invalid.')
      asespacegroup = ase.spacegroup.Spacegroup(sg)
      self.dftcalc.symop = asespacegroup.get_rotations()
      logger.info('  Spacegroup: {}'.format(sg))
      self.dftcalc.invsymop = np.linalg.inv(self.dftcalc.symop)
      self.dftcalc.nsym = self.dftcalc.symop.shape[0]
      self.dftcalc._computeMomentumSymmetries()
      logger.info('  Number of symmetry operations: {}'.format(self.dftcalc.nsym))

    for ispin in range(self.dftcalc.spins):
      ''' provide the boltztrap module with the data for spin: ispin
          and append the resulting velocities / curvatures to the internal list '''
      self._interp(ispin)

    self._symmetrize()
    logger.info('BoltzTrap2: Interpolation successful.')


  def _interp(self, spin):
    '''
    Standard BoltzTrap2 Library interface to interpolate
    the band energies, band velocities and band curvatures
    '''

    # we disable BTP logging if we are not in DEBUG mode!
    # disable = not (logger.getEffectiveLevel() == logging.DEBUG)
    disable = True

    if disable:
      logging.disable( sys.maxsize if sys.version_info >= (3,0) else sys.maxint)

    ''' Under GPLv3 licensed and modified BoltzTraP2 DFTData class
        We adopted the __init__ method which we feed with our electronic structure data
    '''
    self.data = DFTData(self.dftcalc.aseobject, self.dftcalc.weightsum, self.dftcalc.kpoints, \
                        self.dftcalc.mu, self.dftcalc.energies[spin], self.dftcalc.charge)

    # this was the old direct access
    # self.data = BTP.DFTData(self.dftcalc.directory, derivatives=False) # ignore mommat2 files

    self.equivalences = sphere.get_equivalences(self.data.atoms, self.data.magmom, \
                                                self.niter * len(self.data.kpoints))

    self.coeffs = fite.fitde3D(self.data, self.equivalences)
    self.metadata = serialization.gen_bt2_metadata(self.data, self.data.magmom is not None)

    self.lattvec = self.data.get_lattvec()

    if disable:
      logging.disable(logging.NOTSET)

    ''' If we provide a new mesh, use it instead of the original one
        for the moment: this will be reducible
        for spin polarized, avoid the second computation
    '''

    if self.mesh is not None and spin==0:
      self._generate_mesh()
    elif self.mesh is None:
      self.kpoints = self.data.kpoints

    self.interp_energies, self.interp_velocities, self.interp_curvatures = \
        fite.getBands(self.kpoints, self.equivalences, self.lattvec, self.coeffs, curvature=True)



    # we get the energies on the Hartree scale -> rescale to eV
    self.energies.append(self.interp_energies.transpose(1,0) * units.hartree2eV)
    self.velocities.append(self.interp_velocities.transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)

    # here we dont save unnecessary information
    # because d2/dxdy = d2/dydx
    #
    # 1 4 5
    # - 2 6
    # - - 3

    # my index array to get exactly the values above on one axis
    d2ksave = tuple((np.array([0,1,2,0,0,1]), np.array([0,1,2,1,2,2])))

    # this works but I dont know whether there are better ways to do this
    # looks rather hacky
    tmp = self.interp_curvatures.transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
    tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)
    for i in range(tmp.shape[0]):
      for j in range(tmp.shape[1]):
        view = tmp[i,j,:,:]
        tmp2[i,j,:] = view[d2ksave]

    self.curvatures.append(tmp2)

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

  #  old save routines, before the boltztrap2 interface was simplified

  #def _save(self):
  #  '''
  #  Saving the energies, velocities and curvatures by simply
  #  appending them to the lists we initialized at init time
  #  '''

  #  # we get the energies on the Hartree scale
  #  # rescaling to eV!
  #  self.energies.append(self.interp_energies.transpose(1,0) * units.hartree2eV)
  #  self.velocities.append(self.interp_velocities.transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)

  #  # here we dont save unnecessary information
  #  # because d2/dxdy = d2/dydx
  #  #
  #  # 1 4 5
  #  # - 2 6
  #  # - - 3

  #  # my index array to get exactly the values above on one axis
  #  d2ksave = tuple((np.array([0,1,2,0,0,1]), np.array([0,1,2,1,2,2])))

  #  # this works but I dont know whether there are better ways to do this
  #  # looks rather hacky
  #  tmp = self.interp_curvatures.transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
  #  tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)
  #  for i in range(tmp.shape[0]):
  #    for j in range(tmp.shape[1]):
  #      view = tmp[i,j,:,:]
  #      tmp2[i,j,:] = view[d2ksave]

  #  self.curvatures.append(tmp2)

  #def _save1_separate(self):
  #  '''
  #  Identical to _save1 only applied to cases where we have spin orbit coupling
  #  Wien2K saves the data in one file (energyso or energysoup) where the 'spins' alternate
  #  We perform this only for spin-polarized calculations with spin -orbit coupling where it makes
  #  sense to separate out the energies

  #  For unpolarized SOC calculations leave them as is.
  #  '''

  #  # we get the energies on the Hartree scale
  #  # rescaling to eV!
  #  self.interp_energies = self.interp_energies.transpose(1,0) * units.hartree2eV
  #  self.energies.append(self.interp_energies[:,::2])
  #  self.energies.append(self.interp_energies[:,1::2])

  #  self.interp_velocities = self.interp_velocities.transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom
  #  self.velocities.append(self.interp_velocities[:,::2,:])
  #  self.velocities.append(self.interp_velocities[:,1::2,:])

  #  # here we dont save unnecessary information
  #  # because d2/dxdy = d2/dydx
  #  #
  #  # 1 4 5
  #  # - 2 6
  #  # - - 3

  #  # my index array to get exactly the values above on one axis
  #  d2ksave = tuple((np.array([0,1,2,0,0,1]), np.array([0,1,2,1,2,2])))

  #  # this works but I dont know whether there are better ways to do this
  #  # looks rather hacky
  #  tmp = self.interp_curvatures.transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
  #  tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)
  #  for i in range(tmp.shape[0]):
  #    for j in range(tmp.shape[1]):
  #      view = tmp[i,j,:,:]
  #      tmp2[i,j,:] = view[d2ksave]

  #  self.curvatures.append(tmp2[:,::2])
  #  self.curvatures.append(tmp2[:,1::2])

  # def _save2(self):
  #   '''
  #   In the case of a spin-dependent VASP calculation BTP2
  #   creates arrays which lists the data in order energyup energydn.
  #   We want them to be split, which is what this routines does
  #   '''

  #   nbands = self.interp_energies.shape[0] # this is guaranteed to be even here
  #   self.energies.append(self.interp_energies[:nbands//2,:].transpose(1,0) * units.hartree2eV)
  #   self.energies.append(self.interp_energies[nbands//2:,:].transpose(1,0) * units.hartree2eV)
  #   self.velocities.append(self.interp_velocities[:,:nbands//2,:].transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)
  #   self.velocities.append(self.interp_velocities[:,nbands//2:,:].transpose(2,1,0) * units.hartree2eV * units.bohr2angstrom)
  #   # the last two elements dont matter here, since its symmetric anyways

  #   # my index array
  #   d2ksave = tuple((np.array([0,1,2,0,0,1]), np.array([0,1,2,1,2,2])))

  #   # some numpy magic to turn the 3x3 array into the only 6 necessary entries
  #   tmp = self.interp_curvatures[:,:,:nbands//2,:].transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
  #   tmp2 = np.zeros((tmp.shape[0], tmp.shape[1], 6), dtype=np.float64)

  #   for i in range(tmp.shape[0]):
  #     for j in range(tmp.shape[1]):
  #       view = tmp[i,j,:,:]
  #       tmp2[i,j,:] = view[d2ksave]
  #   self.curvatures.append(tmp2)

  #   tmp = self.interp_curvatures[:,:,nbands//2:,:].transpose(3,2,0,1) * units.hartree2eV * units.bohr2angstrom**2
  #   for i in range(tmp.shape[0]):
  #     for j in range(tmp.shape[1]):
  #       view = tmp[i,j,:,:]
  #       tmp2[i,j,:] = view[d2ksave]
  #   self.curvatures.append(tmp2)

  def _generate_mesh(self, shift=False):
    '''
    Generate new moentum mesh for which we generate
    energies / velocities / curvatures
    '''

    _kmeshx = np.linspace(0,1,self.mesh[0],endpoint=False)
    _kmeshy = np.linspace(0,1,self.mesh[1],endpoint=False)
    _kmeshz = np.linspace(0,1,self.mesh[2],endpoint=False)

    if shift:
      self._kmeshshift = []
      for ik in [self.mesh[0],self.mesh[1],self.mesh[2]]:
        if ik > 1:
          self._kmeshshift.append(1./ik/2.)
        else:
          self._kmeshshift.append(0.0)
      self._kmeshshift = np.array(self._kmeshshift, dtype=np.float64)

    # the way these points are ordered is important for the indexing below
    kpoints = []
    for ikx in _kmeshx:
      for iky in _kmeshy:
        for ikz in _kmeshz:
          kpoints.append([ikx,iky,ikz])
    kpoints = np.array(kpoints, dtype=np.float64)
    if shift: kpoints += self._kmeshshift[None,:]

    unique  = np.ones((self.mesh[0]*self.mesh[1]*self.mesh[2]), dtype=int)
    mult    = np.zeros((self.mesh[0]*self.mesh[1]*self.mesh[2]), dtype=int)
    irrk    = 0

    mesh_warning = False
    if self.dftcalc.irreducible and self.dftcalc.nsym > 1:
      # logger.info('Generating irreducible kpoints:')

      for ik in range(np.product(self.mesh)):
        # progressBar(ik+1,self.nkp,status='k-points')

        if unique[ik] == 0: continue # skip if we already went there via symmetry
        irrk += 1    # new point -> increase irreducible counter
        mult[ik] = 1 # reset multiplicity counter

        ''' generate all the symmetry related k-points in the Brillouin zone
            Python modulo via % is implemented as floored division -> -0.2 % 1 = 0.8
        '''
        knew = np.einsum('nji,j->ni',self.dftcalc.momsymop,kpoints[ik,:])
        kmod = knew%1
        # ''' in order to index properly and if kshift is applied , shift back '''
        if shift:
          kmod -= self._kmeshshift
        ''' round to neareast integer '''
        kround = np.rint(kmod * np.array([self.mesh[0],self.mesh[1],self.mesh[2]])[None,:])
        ''' exact floating calculation '''
        kexact = kmod * np.array([self.mesh[0],self.mesh[1],self.mesh[2]])[None,:]
        ''' only use the values that transform properly on all three axes '''
        mask = np.all(np.isclose(kround,kexact),axis=1)
        if not np.all(mask):
          mesh_warning = True
        ''' apply the mask to filter '''
        kmask = kround[mask]
        ''' get the hash index '''
        kindex = (kmask[:,2] + \
                  kmask[:,1] * self.mesh[2] + \
                  kmask[:,0] * self.mesh[2] * self.mesh[1]).astype(int)
        ''' remove the k-points connected via symmetry and increase the multiplicity accordingly '''
        for ikk in kindex:
          if ikk <= ik: continue
          if unique[ikk]:
            unique[ikk] = 0
            mult[ik] += 1

      if mesh_warning:
        logger.critical('\n\n############\nProvided momentum mesh does not conform with symmetry.\n' + \
                        'Accuracy of results cannot be guaranteed.\n############\n')

      self.nkp                     = irrk
      self.nkx, self.nky, self.nkz = self.mesh
      self.kpoints                 = kpoints[unique>0]
      self.multiplicity            = mult[unique>0]
      self.weights                 = self.dftcalc.weightsum * self.multiplicity / np.sum(self.multiplicity)
      self.nsym                    = self.dftcalc.nsym
      self.momsymop                = self.dftcalc.momsymop
      self.symop                   = self.dftcalc.symop
      self.invsymop                = self.dftcalc.invsymop
      self.irreducible             = True
      logger.info('Generated new irreducible kmesh with {} irreducible kpoints'.format(self.nkp))

    else:
      self.kpoints                 = kpoints
      self.nkp                     = np.product(self.mesh)
      self.nkx, self.nky, self.nkz = self.mesh
      self.multiplicity            = np.ones((self.nkp,), dtype=int)
      self.weightsum               = self.dftcalc.weightsum
      self.weights                 = self.dftcalc.weightsum * self.multiplicity / np.sum(self.multiplicity)
      self.irreducible             = False
      self.nsym                    = 1
      self.momsymop                = np.array([[[1,0,0],[0,1,0],[0,0,1]]], dtype=np.float64)
      self.symop                   = np.array([[[1,0,0],[0,1,0],[0,0,1]]], dtype=np.float64)
      self.invsymop                = np.array([[[1,0,0],[0,1,0],[0,0,1]]], dtype=np.float64)
      logger.info('Generated new reducible kmesh with {} irreducible kpoints'.format(self.nkp))

  def _symmetrize(self):
    '''
    If we want to be able to use these elements on the irreducible k-grid
    we need to symmetrized them [akin to what wien2k does]
    '''

    logger.info('BoltzTrap2: Symmetrizing band derivatives.')

    d2ksave = tuple((np.array([0,1,2,0,0,1]), np.array([0,1,2,1,2,2])))

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

      BopticalDiag = np.zeros((nkp,nbands,3,3,3), dtype=np.complex128)

      if self.dftcalc.opticdiag:
        # use the number of elements from the optics
        ioptical = self.dftcalc.opticalDiag[0].shape[-1]
      else:
        # use 3 or 6 according to our own rules.
        if self.dftcalc.ortho:
          ioptical = 3
        else:
          ioptical = 6
      # this is 3 6 or 9
      opticalDiag  = np.zeros((nkp,nbands,ioptical), dtype=np.float64)

      nsym = self.dftcalc.nsym
      rotsymop  = np.einsum('ij,njk,kl->nil',np.linalg.inv(self.dftcalc.kvec),self.dftcalc.invsymop,self.dftcalc.kvec)
      rotsymopT = np.einsum('ij,njk,kl->nli',np.linalg.inv(self.dftcalc.kvec),self.dftcalc.invsymop,self.dftcalc.kvec)

      for ikp in range(nkp):
        progressBar(ikp+1,nkp, status='k-points', prefix=prefix)


        vel     = self.velocities[ispin][ikp,:,:] # nbands, 3
        cur     = self.curvatures[ispin][ikp,:,:] # nbands, 6

        # put the curvatures in symmetric matrix form
        curmat  = np.zeros((nbands,3,3), dtype=np.float64)
        curmat[:, [0,1,2,1,2,2], [0,1,2,0,0,1]] = cur[:,:]
        curmat[:, [0,0,1], [1,2,2]] = curmat[:, [1,2,2], [0,0,1]]

        vk = np.einsum('nij,bj->bni',rotsymop,vel)
        vk_conj = np.conjugate(vk)
        ck = np.einsum('nij,bjk,nkl->bnil',rotsymop,curmat,rotsymopT) # bands, bands, nsym, 3, 3

        ''' these are band interpolation, nothing complex can appear here '''
        vk2 = vk_conj[:,:,[0,1,2,0,0,1]] * vk[:,:,[0,1,2,1,1,2]]
        vk2 = np.mean(vk2,axis=1).real # symmetrize over the squares

        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bnx,bnj,bnyi->bnxyz',levmatrix,vk_conj,vk,ck)
        mb = np.mean(mb,axis=1)

        if ioptical==3:
          opticalDiag[ikp,:,:] = vk2[...,:3]
        else:
          opticalDiag[ikp,:,:6] = vk2[...]

        BopticalDiag[ikp,:,:,:,:] = mb

      self.opticalDiag.append(opticalDiag)
      self.BopticalDiag.append(BopticalDiag)

    # if we need the peierls approximation
    self.opticalBandMin = 0
    self.opticalBandMax = self.velocities[0].shape[1]


class DFTData:
    """
      Objects of this class hold structural and dynamical information from DFT
      results in any supported format.
    """

    def __init__(self, aseobject, weightsum, kpoints, mu, energies, charge):
      """
        Create BoltzTraP2 DFTData object given our electronicstructure objects

         We provide the data as explicit arguments to avoid any spin related problems
         INFO: We transform our internal units (eV) to the BoltzTraP units (Ha)
         nota bene: 1 Hartree = 27.211407953 eV
      """

      self.sysname   = "DFT_to_BTP2"
      self.atoms     = aseobject
      self.dosweight = weightsum
      self.kpoints   = kpoints.copy()
      self.fermi     = mu * 0.0367492929 # eV to Hartree
      self.ebands    = energies.T.copy() * 0.0367492929 # eV to Hartree
      self.mommat    = None
      self.magmom    = None
      self.nelect    = charge
      self.source    = "LinReTraCe"

    # def __init__(self, directory, derivatives=False, *args, **kwargs):
    #     """Create a DFTData object."""
    #     for label, loader in loaders[::-1]:
    #         BoltzTraP2.misc.info("looking for a {} calculation".format(label))
    #         try:
    #             loaded = loader(directory, *args, **kwargs)
    #         except LoaderError as e:
    #             BoltzTraP2.misc.info("error in {} loader: {}".format(label, e))
    #             continue
    #         self.source = label
    #         break
    #     else:
    #         raise ValueError(
    #             "no calculation found in directory {}".format(directory)
    #         )
    #     BoltzTraP2.misc.info(
    #         "successfully loaded a {} calculation".format(self.source)
    #     )
    #     # Try to copy all relevant attributes from the loader
    #     if derivatives:
    #         try:
    #             self.mommat = loaded.mommat
    #         except AttributeError:
    #             raise ValueError(
    #                 "no derivative information found in directory {}".format(
    #                     directory
    #                 )
    #             )
    #     else:
    #         try:
    #             loaded.mommat
    #         except AttributeError:
    #             pass
    #         else:
    #             BoltzTraP2.misc.info(
    #                 "derivative information will be discarded"
    #             )
    #         self.mommat = None
    #     try:
    #         self.sysname = loaded.sysname
    #         self.atoms = loaded.atoms
    #         self.dosweight = loaded.dosweight
    #         self.kpoints = loaded.kpoints
    #         self.fermi = loaded.fermi
    #         self.ebands = loaded.ebands
    #     except AttributeError:
    #         raise ValueError(
    #             "some essential piece of information was not loaded"
    #         )

    #     # Warn the user if the spin up and spin down Fermi energies
    #     # are different in CASTEP.
    #     try:
    #         self.castep_fermi_mismatch = loaded.castep_fermi_mismatch
    #         if self.castep_fermi_mismatch:
    #             BoltzTraP2.misc.info(
    #                 "CASTEP WARNING: "
    #                 "Different spin up and spin down Fermi energy."
    #                 "\nProceeding with spin up Fermi energy. Transport results might"
    #                 " be inaccurate."
    #             )
    #     except AttributeError:
    #         pass

    #     BoltzTraP2.misc.info("Fermi energy:", self.fermi)
    #     # If no initial magnetic moments are provided by the loader, assume a
    #     # non-spin-polarized calculation.
    #     try:
    #         self.magmom = loaded.magmom
    #     except AttributeError:
    #         self.magmom = None
    #         BoltzTraP2.misc.info("Assuming a non-spin-polarized calculation")
    #     # If the number of valence electrons has not been set yet, compute it
    #     # from the bands.
    #     try:
    #         self.nelect = loaded.nelect
    #     except AttributeError:
    #         degeneracies = BoltzTraP2.sphere.calc_reciprocal_degeneracies(
    #             self.atoms, self.magmom, self.kpoints
    #         )
    #         weights = degeneracies.astype(np.float64) / degeneracies.sum()
    #         occupancy = (loaded.ebands < loaded.fermi).astype(np.intc)
    #         self.nelect = round(self.dosweight * (occupancy * weights).sum())

    def bandana(self, emin=-np.inf, emax=np.inf):
      bandmin = np.min(self.ebands, axis=1)
      bandmax = np.max(self.ebands, axis=1)
      ntoolow = np.count_nonzero(bandmax <= emin)
      accepted = np.logical_and(bandmin < emax, bandmax > emin)
      BoltzTraP2.misc.info("BANDANA output")
      for iband in range(len(self.ebands)):
          BoltzTraP2.misc.info(
              iband, bandmin[iband], bandmax[iband], accepted[iband]
          )
      self.ebands = self.ebands[accepted]
      if self.mommat is not None:
          self.mommat = self.mommat[:, accepted, :]
      # Removing bands may change the number of valence electrons
      self.nelect -= self.dosweight * ntoolow
      return accepted

    def get_lattvec(self):
      try:
          self.lattvec
      except AttributeError:
          self.lattvec = self.atoms.get_cell().T * Angstrom
      return self.lattvec

    def get_volume(self):
      try:
          self.UCVol
      except AttributeError:
          lattvec = self.get_lattvec()
          self.UCvol = np.abs(np.linalg.det(lattvec))
      return self.UCvol

    def get_formula_count(self):
     """Return the number of irreducible formulas in the unit cell.

     Useful for computing molar quantities.
     """
     counts = collections.Counter(self.atoms.get_chemical_symbols())
     return functools.reduce(math.gcd, counts.values())




# old class that defines the W2k loader for BoltzTraP2 from our internal file list
# class MetaW2kLoader(BTP.GenericWien2kLoader):
#   '''
#   BoltzTrap Custom Wien2k Loader.
#   After setting the class variables one can register the Loader
#   und use the provided custom files.
#   The usual Wien2kLoader can only access the energy and energyso files.
#   We also want to access up and dn files in spin-polarized calculations.
#   '''

#   # define class variables
#   weightsum = None
#   fscf      = None
#   fstruct   = None
#   fenergy   = None

#   # access them here
#   def __init__(self, directory):
#       super(MetaW2kLoader, self).__init__(MetaW2kLoader.case, \
#                                        MetaW2kLoader.weightsum, \
#                                        MetaW2kLoader.fscf, \
#                                        MetaW2kLoader.fstruct, \
#                                        MetaW2kLoader.fenergy)

#   @classmethod
#   def setfiles(cls, case, weightsum, fscf, fstruct, fenergy):
#     cls.case      = case
#     cls.weightsum = weightsum
#     cls.fscf      = fscf
#     cls.fstruct   = fstruct
#     cls.fenergy   = fenergy
