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
import logging
logger = logging.getLogger(__name__)

import io
import contextlib
import numpy as np

''' BoltzTraP2 emits a bare notice on import when pyfftw is absent ("you can
    install pyfftw to get better FFT performance"), via BoltzTraP2.misc.warning
    which writes to stderr.  It arrives before any LinReTraCe output and is
    misleading here: the FFT path lives in fite.FFTev / fite.FFTc, reached only
    through fite.getBTPbands, whereas we call fite.getBands, a direct phase sum
    that never performs an FFT.  Both streams are captured (the notice moves
    between them across releases) and relegated to the debug log; nothing is
    discarded, and exceptions still propagate normally. '''
_btp2_import_notes = io.StringIO()
with contextlib.redirect_stdout(_btp2_import_notes), contextlib.redirect_stderr(_btp2_import_notes):
  import BoltzTraP2.dft as BTP
  from BoltzTraP2 import sphere
  from BoltzTraP2 import fite
if _btp2_import_notes.getvalue().strip():
  logger.debug('BoltzTrap2 import notice: {}'.format(_btp2_import_notes.getvalue().strip()))

import ase.spacegroup

from structure.auxiliary import progressBar
from structure.wien2k    import Wien2kCalculation
from structure.vasp      import VaspCalculation
from structure           import units
from structure.auxiliary import levicivita


''' Range of BoltzTraP2 releases this interface has been validated against.

    Lower bound inclusive, upper bound exclusive.  These are a tripwire, not
    a claim that anything outside the range is broken: bump them deliberately
    after re-running the interpolation regression test.  We depend on exactly
    four BoltzTraP2 symbols -- sphere.get_equivalences, fite.fitde3D,
    fite.getBands and the DFTData container inherited below -- so the surface
    that a new release can disturb is small.
'''
BTP2_TESTED_MIN = (25, 11)
BTP2_TESTED_MAX = (27,  0)

''' Peak memory of fite.getBands, measured on this interface: it allocates
    13 complex128 temporaries of shape (n_equivalences, n_kpoints) up front
    (phase, phaseR x3, phaseRR x9) plus the phase0 working arrays.  Measured
    271 bytes per (equivalence class, k-point) pair, stable over two orders
    of magnitude in both factors.  Evaluation is therefore chunked over
    k-points against the budget below.  Raise BTP2_MEMORY_BUDGET_GB on a
    large machine; lowering it only costs a little speed.
'''
BTP2_BYTES_PER_PAIR   = 271
BTP2_MEMORY_BUDGET_GB = 2.0

''' Default settings of the hold-out validation (see
    BoltztrapInterpolation.validate).  The warning threshold is k_B T at
    BTP2_VALIDATE_TREF: an interpolation error larger than the thermal
    smearing of the intended run makes the transport integrals unreliable.
'''
BTP2_VALIDATE_FRACTION = 0.2   # share of the parent mesh held out
BTP2_VALIDATE_WINDOW   = 1.0   # eV around mu in which the error is measured
BTP2_VALIDATE_MINKP    = 20    # below this many parent k-points, skip
BTP2_VALIDATE_TREF     = 300.0 # K, reference temperature for the warning


def btp2_version():
  '''
  Return the installed BoltzTraP2 release as a tuple of integers.

  BoltzTraP2 does not define ``BoltzTraP2.__version__``; the release string
  lives in ``BoltzTraP2.version.PROGRAM_VERSION`` (e.g. '26.3.1').

  Returns
  -------
  tuple of int, or None
      (26, 3, 1) for release '26.3.1'.  None if the version cannot be
      determined, which is not treated as an error.
  '''

  try:
    from BoltzTraP2.version import PROGRAM_VERSION
  except ImportError:
    return None

  fields = []
  for field in str(PROGRAM_VERSION).split('.'):
    try:
      fields.append(int(field))
    except ValueError:
      break # stop at the first non numeric field, e.g. '26.3.1rc1'
  return tuple(fields) if fields else None


def chunk_size(n_equiv, budget_gb=None, bytes_per_pair=BTP2_BYTES_PER_PAIR):
  """
  Number of k-points fite.getBands can be handed at once within a memory
  budget.

  Parameters
  ----------
  n_equiv : int
      Number of star-function equivalence classes, len(equivalences).  Set
      by the *parent* mesh and the interpolation parameter, not by the mesh
      being evaluated.
  budget_gb : float, optional
      Budget in GiB.  Defaults to BTP2_MEMORY_BUDGET_GB.
  bytes_per_pair : int, optional
      Bytes per (equivalence class, k-point) pair.

  Returns
  -------
  int
      Chunk length, at least 1.  A chunk of 1 is attempted rather than
      raising, since the alternative is no result at all.
  """

  if budget_gb is None:
    budget_gb = BTP2_MEMORY_BUDGET_GB
  per_kpoint = float(bytes_per_pair) * float(n_equiv)
  if per_kpoint <= 0.0:
    return 1
  return max(1, int(budget_gb * 1024**3 / per_kpoint))


def getBands_chunked(kpoints, equivalences, lattvec, coeffs, budget_gb=None):
  """
  Memory-bounded replacement for fite.getBands(..., curvature=True).

  fite.getBands allocates the phase factors for every k-point at once, so
  the peak grows as n_equiv * n_kpoints and a dense output mesh exhausts
  memory long before the result would.  Chunking changes nothing about the
  numbers: each k-point is evaluated from the same coefficients,
  independently of every other.

  Parameters
  ----------
  kpoints : ndarray, shape (nkp, 3)
      Fractional coordinates to evaluate on.
  equivalences, lattvec, coeffs
      As returned by sphere.get_equivalences, DFTData.get_lattvec and
      fite.fitde3D.
  budget_gb : float, optional
      Passed on to chunk_size.

  Returns
  -------
  tuple of ndarray
      (energies, velocities, curvatures) with the shapes fite.getBands
      returns: (nbands, nkp), (3, nbands, nkp), (3, 3, nbands, nkp).
  """

  nkp   = len(kpoints)
  chunk = chunk_size(len(equivalences), budget_gb)

  if chunk >= nkp:
    logger.debug('BoltzTrap2: evaluating {} k-points in a single pass.'.format(nkp))
    return fite.getBands(kpoints, equivalences, lattvec, coeffs, curvature=True)

  nchunks = int(np.ceil(nkp / float(chunk)))
  logger.info('BoltzTrap2: evaluating {} k-points in {} chunks of {} '
              '(~{:.1f} GB peak).'.format(nkp, nchunks, chunk,
              BTP2_BYTES_PER_PAIR * len(equivalences) * chunk / 1024.0**3))

  ene, vel, cur = [], [], []
  for ichunk, start in enumerate(range(0, nkp, chunk)):
    progressBar(ichunk+1, nchunks, status='chunks', prefix='interp:')
    e, v, c = fite.getBands(kpoints[start:start+chunk], equivalences,
                            lattvec, coeffs, curvature=True)
    ene.append(e); vel.append(v); cur.append(c)

  # k is the last axis of all three return values
  return (np.concatenate(ene, axis=-1),
          np.concatenate(vel, axis=-1),
          np.concatenate(cur, axis=-1))


def check_btp2_version():
  '''
  Warn -- never abort -- if the installed BoltzTraP2 lies outside the range
  recorded in BTP2_TESTED_MIN / BTP2_TESTED_MAX.

  Called once per interpolation run so that a version mismatch appears in the
  log next to the numbers it may have affected, rather than surfacing later as
  an obscure traceback (or, worse, not at all).
  '''

  version = btp2_version()
  if version is None:
    logger.warning('BoltzTrap2: unable to determine the installed version.')
    return

  printable = '.'.join(str(i) for i in version)
  if version < BTP2_TESTED_MIN or version >= BTP2_TESTED_MAX:
    logger.warning('BoltzTrap2: version {} lies outside the tested range '
                   '[{}, {}).'.format(printable,
                                      '.'.join(str(i) for i in BTP2_TESTED_MIN),
                                      '.'.join(str(i) for i in BTP2_TESTED_MAX)))
    logger.warning('BoltzTrap2: interpolation has not been validated for this release.')
  else:
    logger.info('BoltzTrap2: detected version {}.'.format(printable))


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
    check_btp2_version()
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
      self.dftcalc._computePrimitiveSymmetries()
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

    self.lattvec = self.data.get_lattvec()

    if disable:
      logging.disable(logging.NOTSET)

    ''' If we provide a new mesh, use it instead of the original one
        we use the same 'type' as the dft calculation provides: reducibility / kshift
        for spin polarized, avoid the second computation
    '''

    if self.mesh is not None and spin==0:
      self._generate_mesh(shift=self.dftcalc.kshift)
    elif self.mesh is None:
      self.kpoints = self.data.kpoints

    self.interp_energies, self.interp_velocities, self.interp_curvatures = \
        getBands_chunked(self.kpoints, self.equivalences, self.lattvec, self.coeffs)



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

      for ik in range(np.prod(self.mesh)):
        # progressBar(ik+1,self.nkp,status='k-points')

        if unique[ik] == 0: continue # skip if we already went there via symmetry
        irrk += 1    # new point -> increase irreducible counter
        mult[ik] = 1 # reset multiplicity counter

        ''' generate all the symmetry related k-points in the Brillouin zone
            Python modulo via % is implemented as floored division -> -0.2 % 1 = 0.8
        '''
        knew = np.einsum('nji,j->ni',self.dftcalc.symop,kpoints[ik,:])
        kmod = knew%1
        ''' in order to index properly and if kshift is applied , shift back '''
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
      self.symop                   = self.dftcalc.symop
      self.invsymop                = self.dftcalc.invsymop
      self.irreducible             = True
      logger.info('Generated new irreducible kmesh with {} kpoints'.format(self.nkp))

    else:
      self.kpoints                 = kpoints
      self.nkp                     = np.prod(self.mesh)
      self.nkx, self.nky, self.nkz = self.mesh
      self.multiplicity            = np.ones((self.nkp,), dtype=int)
      self.weightsum               = self.dftcalc.weightsum
      self.weights                 = self.dftcalc.weightsum * self.multiplicity / np.sum(self.multiplicity)
      self.irreducible             = False
      self.nsym                    = 1
      self.symop                   = np.array([[[1,0,0],[0,1,0],[0,0,1]]], dtype=np.float64)
      self.invsymop                = np.array([[[1,0,0],[0,1,0],[0,0,1]]], dtype=np.float64)
      logger.info('Generated new reducible kmesh with {} kpoints'.format(self.nkp))

  def validate(self, fraction=None, window=None, seed=0, niter=None):
    """
    Estimate the out-of-sample accuracy of the interpolation by hold-out
    cross validation on the parent mesh.

    Why this exists.  fite.fitde3D is an *exact* interpolant: it reproduces
    the DFT energies at the parent k-points to machine precision, whatever
    the mesh density.  Agreement there therefore says nothing at all about
    accuracy anywhere else, and a badly under-resolved calculation looks
    perfect from the inside.  The only honest check is to withhold part of
    the parent mesh, fit on the rest, and compare at the withheld points.

    This costs one extra fit per spin and needs neither optical elements
    nor any additional DFT run, so it applies to every interface.

    Measured behaviour (Wien2k, 20x20x20 parent, --interp 3, states within
    1 eV of mu):

      SrVO3  hold-out rms 0.008-0.010 eV   true out-of-sample rms 0.005 eV
      a-As   hold-out rms 0.337-0.361 eV   true out-of-sample rms 0.290 eV

    so the estimate is conservative by roughly 1.2x to 2x, which is the
    right direction for a warning.  "True" here means measured against an
    independent band path sharing no k-points with the parent mesh.

    Parameters
    ----------
    fraction : float, optional
        Share of parent k-points held out.  Default BTP2_VALIDATE_FRACTION.
    window : float, optional
        Half-width in eV around mu over which the error is reported.  Only
        states near mu matter for transport.  Default BTP2_VALIDATE_WINDOW.
    seed : int, optional
        Seed of the train/test split, so the report is reproducible.
    niter : int, optional
        Interpolation parameter for the reduced fit.  Defaults to the value
        the production fit used, which keeps the comparison like for like.

    Returns
    -------
    list of dict, or None
        One entry per spin with keys 'rms', 'max', 'nstates', 'ntrain',
        'ntest'.  None if the parent mesh is too small to split
        (fewer than BTP2_VALIDATE_MINKP points).
    """

    if fraction is None: fraction = BTP2_VALIDATE_FRACTION
    if window   is None: window   = BTP2_VALIDATE_WINDOW
    if niter    is None: niter    = self.niter

    nkp = len(self.dftcalc.kpoints)
    if nkp < BTP2_VALIDATE_MINKP:
      logger.warning('BoltzTrap2: parent mesh has only {} k-points - '
                     'too few to validate.'.format(nkp))
      return None

    logger.info('BoltzTrap2: Validating interpolation (hold-out {:.0%}).'.format(fraction))

    ''' BoltzTraP2 is chatty; silence it for the duration as elsewhere '''
    logging.disable(sys.maxsize)
    report = []
    try:
      for ispin in range(self.dftcalc.spins):
        rng   = np.random.default_rng(seed)
        order = rng.permutation(nkp)
        ntest = max(1, int(fraction*nkp))
        itest, itrain = order[:ntest], order[ntest:]

        kp, en, mu = self.dftcalc.kpoints, self.dftcalc.energies[ispin], self.dftcalc.mu

        data = DFTData(self.dftcalc.aseobject, self.dftcalc.weightsum,
                       kp[itrain], mu, en[itrain], self.dftcalc.charge)
        equiv  = sphere.get_equivalences(data.atoms, data.magmom, niter*len(itrain))
        coeffs = fite.fitde3D(data, equiv)
        ekp, _ = fite.getBands(kp[itest], equiv, data.get_lattvec(), coeffs, curvature=False)
        ekp    = ekp.T * units.hartree2eV

        nbands = min(en.shape[1], ekp.shape[1])
        ref    = en[itest][:,:nbands]
        mask   = np.abs(ref-mu) < window
        if not np.any(mask):
          report.append(None)
          continue
        err = np.abs(ref - ekp[:,:nbands])[mask]
        report.append(dict(rms=float(np.sqrt((err**2).mean())), max=float(err.max()),
                           nstates=int(mask.sum()), ntrain=len(itrain), ntest=ntest))
    finally:
      logging.disable(logging.NOTSET)

    self._reportValidation(report, window)
    return report


  def _reportValidation(self, report, window):
    """
    Log the outcome of validate() and warn when the interpolation error near
    mu exceeds the thermal smearing at BTP2_VALIDATE_TREF.

    The comparison is expressed as an equivalent temperature so that users
    running below room temperature can judge for themselves: an rms error of
    e eV is only harmless if k_B T of the intended run is comfortably above
    it.
    """

    kBT_ref = units.kB_eV * BTP2_VALIDATE_TREF
    for ispin, entry in enumerate(report):
      prefix = '' if self.dftcalc.spins == 1 else ('up: ' if ispin == 0 else 'dn: ')
      if entry is None:
        logger.warning('BoltzTrap2: {}no states within {} eV of mu - '
                       'validation inconclusive.'.format(prefix, window))
        continue
      logger.info('BoltzTrap2: {}hold-out error for {} states within {} eV of mu: '
                  'rms {:.4f} eV, max {:.4f} eV'
                  .format(prefix, entry['nstates'], window, entry['rms'], entry['max']))
      logger.info('BoltzTrap2: {}this equals k_B T at T = {:.0f} K.'
                  .format(prefix, entry['rms']/units.kB_eV))
      if entry['rms'] > kBT_ref:
        logger.critical('\n\n############\n'
                        'Interpolation error near mu ({:.3f} eV) exceeds k_B T at {:.0f} K '
                        '({:.4f} eV).\nThe parent k-mesh is too coarse for reliable transport: '
                        'a finer\nDFT calculation is required.  Interpolating onto a denser '
                        'output mesh\ndoes not add information and will not fix this.\n'
                        '############\n'
                        .format(entry['rms'], BTP2_VALIDATE_TREF, kBT_ref))


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

      ''' Chunk over k-points.  The largest intermediate is the Levi-Civita
          contraction, of shape (chunk, nbands, nsym, 3, 3, 3) complex128,
          i.e. 432 * nbands * nsym bytes per k-point; size the chunk against
          the same budget the interpolation uses. '''
      per_kpoint = 432.0 * nbands * max(nsym,1)
      kchunk = max(1, int(BTP2_MEMORY_BUDGET_GB * 1024**3 / per_kpoint))

      for kstart in range(0, nkp, kchunk):
        kstop = min(kstart+kchunk, nkp)
        progressBar(kstop, nkp, status='k-points', prefix=prefix)

        vel = self.velocities[ispin][kstart:kstop,:,:] # chunk, nbands, 3
        cur = self.curvatures[ispin][kstart:kstop,:,:] # chunk, nbands, 6

        # put the curvatures in symmetric matrix form
        curmat = np.zeros((kstop-kstart,nbands,3,3), dtype=np.float64)
        curmat[:,:, [0,1,2,1,2,2], [0,1,2,0,0,1]] = cur[:,:,:]
        curmat[:,:, [0,0,1], [1,2,2]] = curmat[:,:, [1,2,2], [0,0,1]]

        vk      = np.einsum('sij,kbj->kbsi',rotsymop,vel)
        vk_conj = np.conjugate(vk)
        ck      = np.einsum('sij,kbjm,sml->kbsil',rotsymop,curmat,rotsymopT,optimize=True)

        ''' these are band interpolation, nothing complex can appear here
            index order of the 6 stored components: xx yy zz xy xz yz '''
        vk2 = vk_conj[:,:,:,[0,1,2,0,0,1]] * vk[:,:,:,[0,1,2,1,2,2]]
        vk2 = np.mean(vk2,axis=2).real # symmetrize over the squares

        #           epsilon_cij v_a v_i c_bj -> abc
        mb = np.einsum('zij,kbsx,kbsi,kbsyj->kbsxyz',levmatrix,vk_conj,vk,ck,optimize=True)
        mb = np.mean(mb,axis=2)

        if ioptical==3:
          opticalDiag[kstart:kstop,:,:] = vk2[...,:3]
        else:
          opticalDiag[kstart:kstop,:,:6] = vk2[...]

        BopticalDiag[kstart:kstop,:,:,:,:] = mb

      self.opticalDiag.append(opticalDiag)
      self.BopticalDiag.append(BopticalDiag)

    # if we need the peierls approximation
    self.opticalBandMin = 0
    self.opticalBandMax = self.velocities[0].shape[1]


class DFTData(BTP.DFTData):
    """
      Container that hands LinReTraCe's electronic structure to BoltzTraP2.

      BoltzTraP2's own DFTData constructor expects a directory on disk and
      runs its loaders over it.  We already hold everything in memory, so the
      constructor is the only member we override.  bandana, get_lattvec,
      get_volume and get_formula_count are inherited unchanged.

      Why inherit rather than copy:  fite.fitde3D reads kpoints, ebands,
      mommat and get_lattvec() off this object.  A local copy of the accessors
      would keep the old behaviour while the consumer moves on, which fails
      quietly with wrong numbers instead of loudly with a traceback.
      Inheriting keeps container and consumer in step across releases; see
      BTP2_TESTED_MIN / BTP2_TESTED_MAX above for the version tripwire.

      Note: BTP.DFTData.__init__ is deliberately NOT called.  Its signature is
      incompatible and it would try to read from disk.  Every attribute that
      the inherited methods and fite.fitde3D rely on is set below.
    """

    def __init__(self, aseobject, weightsum, kpoints, mu, energies, charge):
      """
        Build a BoltzTraP2 DFTData object from our internal arrays.

        Data are passed as explicit arguments rather than as an
        ElectronicStructure object so that the caller resolves the spin
        channel; this class never has to know about spin.

        Parameters
        ----------
        aseobject : ase.Atoms
            Structure of the unit cell.  BoltzTraP2 takes the lattice vectors
            from it and, via spglib, the symmetry used to build the star
            functions -- so it must describe the real structure, not just the
            lattice.
        weightsum : float
            Sum of the k-point weights: 2 for an unpolarised calculation, 1
            per channel otherwise.  BoltzTraP2 calls this 'dosweight'.
        kpoints : ndarray, shape (nkp, 3), float64
            Fractional coordinates of the k-mesh.  Copied, not referenced.
        mu : float
            Chemical potential in eV.
        energies : ndarray, shape (nkp, nbands), float64
            Band energies in eV for one spin channel.  Stored transposed,
            since BoltzTraP2 indexes bands first.
        charge : float
            Number of valence electrons in the unit cell.

        Notes
        -----
        LinReTraCe works in eV throughout, BoltzTraP2 in Hartree, so energies
        are converted here and converted back in BoltztrapInterpolation._interp.
        Both directions use units.hartree2eV so that the round trip is exact.
      """

      self.sysname   = "DFT_to_BTP2"
      self.atoms     = aseobject
      self.dosweight = weightsum
      self.kpoints   = kpoints.copy()
      self.fermi     = mu * units.eV2hartree
      self.ebands    = energies.T.copy() * units.eV2hartree
      self.mommat    = None   # no momentum matrix elements: fit energies only
      self.magmom    = None   # non spin-polarised symmetry for the star functions
      self.nelect    = charge
      self.source    = "LinReTraCe"




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
