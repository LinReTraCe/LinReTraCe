#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import os
import logging
logger = logging.getLogger(__name__)

# python 2 & 3 compatible usage of abstract base classes
if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

import scipy.optimize
import numpy as np

import structure.symmetries.C1

class Converged(Exception):
  def __init__(self, mu):
    super(Converged, self).__init__(self)
    self.mu = mu

class ElectronicStructure(ABC):
  '''
  Abstract Parent class for all electronic structures.
  Here we define all the common elements:
  Number of spins, k-points, multiplicity, weights, etc.

  Methods required to be ran before using the output function in inout.py:
    _calcFermiLevel
    _defineDimensions
  '''

  def __init__(self):
    # information about kmesh
    self.nkp            = 0    # number of k-points (if reducible nkp = nkx.nky.nkz)
    self.nkx            = 0    # number of k-points in x-direction
    self.nky            = 0    # number of k-points in x-direction
    self.nkz            = 0    # number of k-points in x-direction

    self.multiplicity   = None # multiplicity of the k-points
    self.weights        = None # weights of the k-points
    self.weightsum      = 0    # sum of the weights
    self.kpoints        = None # list of the k-points ... shape [nkp,3] float64
    self.kshift         = False# shifted from gamma origin
    self.irreducible    = False# REGULAR irreducible grid: weights are reconstructible
                               # from multiplicity / (nkx.nky.nkz).  This is what
                               # gets written to .kmesh/irreducible and is the ONLY
                               # thing linretrace uses it for -- see main.F90.
    self.symmetrize     = False# moments must be group-averaged over the star of
                               # every k-point (wedge that is NOT a regular grid,
                               # i.e. a refined/custom irreducible mesh)
    self.ortho          = False# orthogonal unit cell

    self.spins          = 1    # number of inequivalent spins we are considering
    self.charge         = 0    # charge in the given bands
    self.mu             = 0    # chemical potential

    self.vol            = 0    # volume of the unit cell in AA^3
    self.rvec           = None # real space lattice vectors (rows)
    self.kvec           = None # reciprocal space lattice vectors (rows)
                               # [i,j] represents the jth element of the ith vector
                               # i.e. first entry -> vector 1 2 3
                               #      second entry -> x y z

    # symmetries
    self.nsym           = structure.symmetries.C1.nsym
    self.symop          = structure.symmetries.C1.symop    # real space symmetry
    self.invsymop       = structure.symmetries.C1.invsymop # inverse of real space symmetry

    self.ndim           = -1   # number of dimensions
    self.dims           = np.array([False,False,False]) # valid dimension, i.e. k_i > 1

    # information about bands
    self.energies       = []   # list entries for energies elements are the respective spins
                               # entries themselves are shape [nkp,bands] float64
    self.energyBandMax  = 0    # band maximum for energies

    # information about optical elements
    self.opticdiag       = False # intra band optical elements
    self.opticalDiag     = []    # list for band-diagonal optical elements
    self.opticfull       = False # full optical elements (inter+intra)
    self.opticalMoments  = []    # list for full optical elements
    # same thing for the b-field quantities
    self.bopticdiag      = False
    self.BopticalDiag    = []
    self.bopticfull      = False
    self.BopticalMoments = []   # list for full optical elements
    self.opticalBandMin  = 0    # band interval minimum for optical elements
    self.opticalBandMax  = 0    # band interval maximum for optical elements

    # if used as model / wannier90 calculation we save those
    self.hk              = None
    self.hvk             = None
    self.hck             = None
    self.Ukohnsham       = None
    self.Uinvkohnsham    = None

  def _computePrimitiveSymmetries(self):
    '''
    Given the __conventional space__ symmetry operations
    we calculate the primitive symmetry operations

    this is necessary e.g. in wien2k (and ase <-> wien2k interface)
    to get the correct matrices for the k-mesh
    '''

    transform = self.kvec @ self.symop @ np.linalg.inv(self.kvec)
    transform_int = np.rint(transform).astype(int)

    ''' this selects non-primitive cells '''
    if np.allclose(transform,transform_int) and \
       np.all(np.abs(transform_int) <= 1) and \
       np.all(np.linalg.det(transform_int) == np.linalg.det(self.symop)):
      logger.debug(' Use transformed momentum matrices')
      if np.all(transform_int == self.symop):
        logger.debug('   n.b.: Transformation did not change matrices')
      self.symop = transform_int
    else:
      logger.debug(' Use unchanged momentum matrices')
      logger.debug('    n.b.: Transformation would have generated invalid matrices')

    self.invsymop = np.linalg.inv(self.symop)

  def _detectKshift(self):
    '''
    Given the current k-points, determine whether the Gamma-point is included
    if included: not shifted
    if not included: shifted
    '''

    gamma = np.zeros(3,)[None,:]
    equals_zero  = np.all(self.kpoints == gamma, axis=1)
    self.kshift = not np.any(equals_zero)

  def _defineDimensions(self):
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

  def _defineDimensionsFromKpoints(self, kpoints, dims=None, tol=1e-10):
    '''
    Set ndim/dims for a custom (irregular) k-mesh.

    Custom meshes carry no regular grid dimensions -- by convention
    nkx=nky=nkz=1 is written for them -- so _defineDimensions would report
    ndim=0 and dims=[F,F,F].  That value is propagated through h5output into
    .unitcell/ndim and is consumed both by linretrace and by postproc (which
    would then zero the off-diagonal Onsager elements and skip the tensor
    inversion).

    If *dims* is given (a length-3 boolean sequence, typically read from the
    .unitcell/dims dataset of the coarse HDF5 the mesh was derived from) it is
    adopted verbatim.  This is the preferred route: the dimensionality of the
    parent calculation is a property of the physical setup, not of the k-mesh,
    and inheriting it guarantees that a refined file is treated exactly like
    the coarse file it came from.

    Only when no *dims* is supplied are the dimensions inferred from the
    spread of the k-points along each axis.  Note that inference cannot see
    cases where the k-mesh is three-dimensional but the projection is not
    (e.g. a single planar d_x2-y2 orbital in a layered cuprate); such setups
    must pass *dims* explicitly.
    '''
    if dims is not None:
      dims = np.asarray(dims, dtype=bool).ravel()
      if dims.shape != (3,):
        raise ValueError('dims must be a length-3 boolean sequence')
      self.dims = dims
      self.ndim = int(np.sum(dims))
      logger.info('Adopted {} dimensions from the parent calculation.'.format(self.ndim))
      return

    kpoints = np.asarray(kpoints, dtype=np.float64)
    if kpoints.ndim != 2 or kpoints.shape[1] != 3:
      raise ValueError('kpoints must be an (N,3) array of fractional coordinates')

    self.ndim = 0
    self.dims = []
    for i in range(3):
      spread = np.max(kpoints[:, i]) - np.min(kpoints[:, i])
      if spread > tol:
        self.ndim += 1
        self.dims.append(True)
      else:
        self.dims.append(False)
    self.dims = np.array(self.dims)
    logger.info('Detected {} dimensions from the custom k-mesh.'.format(self.ndim))

  # Target working set for the automatic block size; see _resolveKblock.
  _WORKING_SET_TARGET = 128 * 1024**2

  @staticmethod
  def _availableMemoryBytes():
    '''Best-effort free memory in bytes; None if it cannot be determined.'''
    try:
      with open('/proc/meminfo') as fh:
        for line in fh:
          if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    except Exception:
      pass
    try:
      return os.sysconf('SC_AVPHYS_PAGES') * os.sysconf('SC_PAGE_SIZE')
    except Exception:
      return None

  def _resolveKblock(self, kblock, memory_budget_gb, nbands, nrpoints, ndir,
                     nsym=1, fullmoments=True, debug=False, label='k-blocking'):
    '''
    Choose the number of k-points evaluated per pass.

    Blocking bounds only the TRANSIENT arrays.  Their cost per k-point is

      24 * nrp                       (r.k and exp(i r.k))
      + ~1056 * nbands^2             (H(k), dH, d2H, U, Uinv, v, c, curmat,
                                      vel2 and the B-field einsum output)

    times a safety factor for the einsum contraction intermediates, and times
    *nsym* on the symmetrising path, where every irreducible k-point is
    expanded into its full star before any linear algebra happens.

    The OUTPUT arrays are independent of the block size and are reported so
    that a run which cannot fit regardless of blocking says so up front.
    '''
    per_k_transient = (24.0 * nrpoints + 1056.0 * nbands**2) * float(nsym)
    per_k_transient *= 1.5   # einsum intermediates

    per_k_output = 8.0 * nbands + nbands**2 * (48.0 + 96.0)
    if fullmoments:
      per_k_output += nbands**2 * (8.0 * ndir + 432.0)
    else:
      per_k_output += nbands * (8.0 * ndir + 432.0)
    if debug:
      per_k_output += 176.0 * nbands**2

    output_gb = per_k_output * self.nkp / 1024**3

    if memory_budget_gb is None:
      avail = self._availableMemoryBytes()
      budget = 0.25 * avail if avail else 2.0 * 1024**3
    else:
      budget = float(memory_budget_gb) * 1024**3

    # The automatic block size targets THROUGHPUT, not the memory ceiling.
    # Beyond a working set of a few hundred MB the block-local arrays stop
    # fitting in cache and the whole pass becomes memory-bandwidth bound:
    # measured on a 512-k-point, 12-orbital, nsym=48 model, 7.1 s at
    # kblock=8 against 19.3 s at kblock=128 and 22.4 s at kblock=512, for
    # bitwise identical output.  A user-supplied memory_budget_gb can only
    # lower this cap, never raise it; pass an explicit kblock to override.
    budget = min(budget, self._WORKING_SET_TARGET)

    if kblock is None:
      kblock = int(max(1, budget // per_k_transient))
    kblock = int(min(max(1, kblock), self.nkp))

    nblocks = -(-self.nkp // kblock)
    logger.info('   {}: {} k-points per pass ({} pass(es)); transient ~{:.2f} '
                'GB, output arrays ~{:.2f} GB.'.format(
                label, kblock, nblocks,
                kblock * per_k_transient / 1024**3, output_gb))
    if output_gb > 8.0:
      logger.warning('The output arrays alone need ~{:.1f} GB and must be '
                     'held in full before writing; k-blocking cannot reduce '
                     'this. Use intraonly if only band-diagonal elements are '
                     'required (drops the dominant nbands^2 B-field term to '
                     'nbands), or refine to fewer k-points.'.format(output_gb))
    return kblock

  @staticmethod
  def _ftHamiltonian(ee, hr, prefactor=None):
    '''
    Fourier transform h(r) -> h(k) as a single BLAS call per Cartesian
    component.

      prefactor is None : out[k,i,j]   = sum_r ee[k,r] hr[r,i,j]
      prefactor given   : out[k,i,j,d] = sum_r prefactor[d,r] ee[k,r] hr[r,i,j]

    Written explicitly rather than as np.einsum('dr,kr,rij->kijd', ...):
    einsum re-plans its contraction path per call, and for the batch sizes
    that arise on the symmetrising path (nkblock * nsym matrices) it picks
    orders whose intermediates grow superlinearly -- measured 33 s at
    kblock=1, 52 s at kblock=2, 7 s at kblock=8 and an out-of-memory kill at
    kblock=512 on the same 512-k-point, 12-orbital model.  Reshaping to a
    plain (nk, nrp) x (nrp, nbands^2) GEMM is deterministic in both time and
    memory.
    '''
    nrp, nb = hr.shape[0], hr.shape[1]
    hrflat  = hr.reshape(nrp, nb*nb)
    if prefactor is None:
      return ee.dot(hrflat).reshape(ee.shape[0], nb, nb)

    ndir = prefactor.shape[0]
    out  = np.empty((ee.shape[0], nb, nb, ndir), dtype=np.complex128)
    for d in range(ndir):
      out[..., d] = (ee * prefactor[d][None, :]).dot(hrflat).reshape(ee.shape[0], nb, nb)
    return out

  @staticmethod
  def _rotateToBandBasis(Uinv, mat, U):
    '''
    Batched basis rotation  out[n,i,l,...] = Uinv[n,i,j] mat[n,j,k,...] U[n,k,l]

    *mat* may carry any number of trailing (Cartesian) axes; they are moved to
    the front so that np.matmul broadcasts over them and every multiplication
    is a batched zgemm.  Replaces np.einsum('nij,njkd,nkl->nild', ...) and its
    two-index variant for the same reason as _ftHamiltonian.
    '''
    extra = mat.ndim - 3
    if extra == 0:
      return np.matmul(Uinv, np.matmul(mat, U))
    fwd = tuple(range(3, mat.ndim)) + (0, 1, 2)
    res = np.matmul(Uinv, np.matmul(np.transpose(mat, fwd), U))
    back = tuple(range(extra, extra + 3)) + tuple(range(extra))
    return np.transpose(res, back)

  @staticmethod
  def _bfieldMoment(levmatrix, velconj, vel, curmat):
    '''
    B-field optical moment

      mb[...,x,y,z] = sum_ij eps[z,i,j] velconj[...,x] vel[...,i] curmat[...,y,j]

    contracted in two steps (a tensordot and a batched matmul) followed by an
    outer product, instead of the four-operand
    np.einsum('zij,bpnx,bpni,bpnyj->bpnxyz', ...).  Same result, deterministic
    cost.
    '''
    # L[...,z,j] = sum_i eps[z,i,j] vel[...,i]
    L = np.tensordot(vel, levmatrix, axes=([-1], [1]))
    # A[...,y,z] = sum_j curmat[...,y,j] L[...,z,j]
    A = np.matmul(curmat, np.swapaxes(L, -1, -2))
    return velconj[..., :, None, None] * A[..., None, :, :]

  def _setCustomSymmetries(self, symop):
    '''
    Attach point-group operations to a custom (irregular) k-mesh, or declare
    the mesh reducible when none are given.

    Why this is needed
    ------------------
    An irreducible mesh covers only a wedge of the Brillouin zone.  Band
    energies are symmetry invariant, so the wedge reproduces the density of
    states directly -- but the optical matrix elements v_i v_j are NOT: the
    wedge is not invariant under the operations that relate the Cartesian
    directions, so a raw weighted sum over it breaks the symmetry of the
    Onsager tensor.  For a cubic single-band model at 8x8x8 the raw wedge sum
    gives (0.2148, 0.3203, 0.2148) where the correct answer is
    (0.25, 0.25, 0.25).

    The symmetrising branch of _computeHk / computeHamiltonian fixes this by
    replacing the moment at every k-point with its group average
    M_sym(k) = (1/nsym) sum_S M(S k).  On a REGULAR grid this is exact: the
    wedge sum then reproduces the full-BZ sum to machine precision, whatever
    the lattice.

    On a REFINED wedge the group average is still the right object -- the
    Onsager tensor comes out with the full point-group symmetry, the weights
    are correct and the total weight is conserved -- but the identity

      sum_c (w_p/n^3) M_sym(k_c)
        = w_p/(n^3 nsym) sum_S sum_c M(S k_c)
        = w_p * <M> over the union of the image cells

    requires that for each S the points {S k_c} are again the sub-cell centres
    of the image cell S.(parent cell), i.e. that the n x n x n sub-tiling of a
    cell is invariant under the point group.  That holds when the operations
    merely permute and flip the fractional axes (cubic, tetragonal,
    orthorhombic) -- the case this was originally validated on -- and FAILS for
    lattices whose operations mix the axes non-trivially, hexagonal above all:
    a 60 degree rotation carries the parallelogram patch of children onto a
    sheared one, so the union of the images is not the refined reducible grid.

    The consequence is a different QUADRATURE, not a wrong one.  Both routes
    converge to the same integral; they simply sample different points, so a
    refined wedge and the refined reducible mesh it descends from must not be
    expected to agree bitwise on a non-orthogonal lattice.  Measured for
    graphene (hexagonal), 8x8x1 subdivided by 3, sum_k w_k |v_xx|^2:

      refined reducible : 1.167906   (identical to the plain 24x24 grid)
      refined wedge     : 1.169197   (0.11% apart)
      converged 192x192 : 1.172925   (both ~0.4% away from this)

    See testsuite/tests/test_refined_mesh_flags.py, which asserts the exact
    equality on the coarse grid and only the quadrature-level agreement after
    refinement.

    Parameters
    ----------
    symop : (nsym,3,3) array or None
        Operations in the convention used by _computeHk, i.e. the reducible
        partners of k are obtained as k_red = P^T . k_irr.  Read back from
        .unitcell/symop of the parent HDF5.  None marks the mesh reducible.
    '''
    ''' A custom mesh is never a REGULAR grid, so its weights can never be
        reconstructed from multiplicity / (nkx.nky.nkz).  self.irreducible is
        exactly that reconstructability flag on the HDF5 side (linretrace uses
        .kmesh/irreducible for nothing else), so it must stay False here and
        the explicit .kmesh/weights must be used downstream.  The symmetrising
        code path is selected by the separate self.symmetrize flag. '''
    self.irreducible = False

    if symop is None:
      self.symmetrize = False
      return

    symop = np.asarray(symop)
    if symop.ndim != 3 or symop.shape[1:] != (3, 3):
      raise ValueError('symop must be an (nsym,3,3) array')
    ''' Point-group operations in the reciprocal-lattice basis are integer
        matrices.  Read back from HDF5 they arrive as float64, which would
        make a refined file carry a float .unitcell/symop where the coarse
        file it descends from carries an int one.  Restore the integer dtype
        whenever the entries are integral -- the einsum in _computeHk is
        agnostic, but the HDF5 schema should not drift between iterations. '''
    symint = np.rint(symop).astype(int)
    if np.allclose(symop, symint, rtol=0.0, atol=1e-10):
      symop = symint
    else:
      symop = np.asarray(symop, dtype=np.float64)

    self.symop      = symop
    self.nsym       = symop.shape[0]
    self.invsymop   = np.linalg.inv(symop)
    self.symmetrize = True
    logger.info('   Custom k-mesh treated as an irreducible wedge: moments will '
                'be group-averaged over {} symmetry operations. The mesh itself '
                'is written as non-uniform (.kmesh/irreducible = False) so that '
                'the explicit weights are used.'.format(self.nsym))

  def _checkBrillouinZone(self, kpoints, name='custom k-mesh'):
    '''
    Map fractional k-points into the primitive reciprocal cell [0,1)^3 and
    report how many had to be wrapped.

    Wrapping itself is physics-neutral -- every quantity is evaluated through
    exp(2 pi i k.R), which is periodic -- but a mesh that leaves [0,1) is a
    symptom worth surfacing: it means the caller built points outside the
    cell, and any downstream code that indexes a regular grid (or compares
    k-points for equality) would then silently disagree with this mesh.
    '''
    kpoints = np.asarray(kpoints, dtype=np.float64)
    outside = np.count_nonzero((kpoints < 0.0) | (kpoints >= 1.0))
    wrapped = np.mod(kpoints, 1.0)
    # guard against 1-eps rounding up to exactly 1.0 in the modulo
    wrapped[wrapped >= 1.0 - 1e-14] = 0.0
    if outside:
      logger.debug('{}: {} coordinate(s) outside [0,1) wrapped into the '
                   'primitive reciprocal cell.'.format(name, outside))
    return wrapped

  def _calcOccupation(self, mu, raiseError=False):
    '''
    Calculate the deviation of the occupation the the given charge in the system.
    Here we use the Fermi function at T=0 (theta function).
    i.e. energies smaller than the chemical potential are fully occupied
    while energies larger than the chemical potential are copletely empty
    '''

    energies = np.zeros((self.spins, self.nkp, self.energyBandMax), dtype=np.float64)
    for ispin in range(self.spins):
      energies[ispin] = self.energies[ispin]

    # T = 0 fermi function -> theta function
    mask = (energies-mu) < 0.0
    energies[mask] = 1.0
    energies[np.logical_not(mask)] = 0.0
    dev = np.sum(energies * self.weights[None,:,None]) - self.charge

    # this is a nasty work-around for the whole gap-thingy
    # to make sure we are where we want to be
    if raiseError and (abs(dev) < np.min(self.weights)/2.):
      raise Converged(mu)

    return dev

  def _calcFermiLevel(self, mu=None):
    '''
    Use the energies and the charge: calculate the T=0 chemical
    potential and determine if a gap exists.
    If one exists, determine the gap size and between which bands
    it lies
    '''

    # get the bisection start and end points
    x0 = np.min(self.energies[0][:,0])
    x1 = np.max(self.energies[0][:,0])

    for ispin in range(self.spins):
      for iband in range(self.energyBandMax):
        ymin = np.min(self.energies[0][:,iband])
        ymax = np.max(self.energies[0][:,iband])
        if ymin < x0:
          x0 = ymin
        if ymax > x1:
          x1 = ymax

    # safety offset
    x0 -= 0.1
    x1 += 0.1

    if mu is None:
      try:
        mu_sol = scipy.optimize.bisect(self._calcOccupation, x0, x1, args=(True,)) # raise Error if Converged
      except Converged as c: # work-around so we always get the band-gap correct
        mu_sol = c.mu
      except ValueError:
        raise Exception('Fermi Level Calculation: Bisection failed.')
    else:
      mu_sol = mu

    if (self.nkp == 1 and abs(mu_sol) <= 1e-11):
      mu_sol = 0
    self.mu = mu_sol

    self.gapped     = []
    self.gap        = []
    self.ecb        = []
    self.evb        = []
    self.cb         = []
    self.vb         = []

    # detect the spin-dependent gap
    for ispin in range(self.spins):
      locgapped = True
      for iband in range(self.energyBandMax):
        ene = self.energies[ispin][:,iband]
        enemin = np.min(ene)
        enemax = np.max(ene)
        if mu_sol > enemin and mu_sol < enemax: # it cuts through this band
          locgapped = False
          self.cb.append(iband)
          self.vb.append(iband)
          break
        if enemin > mu_sol: # band is for the first time above chemical potential
          self.cb.append(iband)
          self.vb.append(iband-1)
          break
        if mu_sol == enemin:
          self.cb.append(iband)
          self.vb.append(iband-1)
          break
        if mu_sol == enemax:
          self.vb.append(iband)
          self.cb.append(iband+1)
          break
      else:
        locgapped = False

      # save the gap data
      if locgapped:
        enevalence    = self.energies[ispin][:,self.vb[ispin]]
        eneconduction = self.energies[ispin][:,self.cb[ispin]]
        gap = np.min(eneconduction) - np.max(enevalence)
        if gap < 1e-13: # required for 'touching' bands like in graphene
          self.gapped.append(False)
          self.gap.append(np.nan)
          self.ecb.append(np.nan)
          self.evb.append(np.nan)
        else:
          self.gapped.append(True)
          self.gap.append(gap)
          self.ecb.append(np.min(eneconduction))
          self.evb.append(np.max(enevalence))
      else:
        self.gapped.append(False)
        self.gap.append(np.nan)
        self.ecb.append(np.nan)
        self.evb.append(np.nan)

    # adjust mu for 1 spin and gapped to be exactly in the middle
    # this is correct for T = 0
    if self.spins==1 and self.gapped[0]:
      self.mu = (self.ecb[0] + self.evb[0]) / 2.
      logger.info("Putting Chemical potential in the middle of the gap.")

    # adjust mu for 2 spins
    # this should be correct for T = 0
    if self.spins==2 and self.gapped[0] and self.gapped[1]:
      self.mu = (min(self.ecb) + max(self.evb)) /2.
      logger.info("Putting Chemical potential in the middle of the common gap.")

    # notify user
    logger.info("Chemical potential: {} [eV]".format(self.mu))
    for ispin in range(self.spins):
      if self.gapped[ispin]:
        logger.info('  Spin: {} / {}: Found energy gap: {} [eV]'.format(ispin+1,self.spins,self.gap[ispin]))
        logger.info('         vbmax: {} [eV] - cbmin: {} [eV]'.format(self.evb[ispin],self.ecb[ispin]))
      else:
        logger.info('  Spin: {} / {}: no energy gap'.format(ispin+1,self.spins))

  @staticmethod
  def distributeWorkLoad(datainterval, processes):
    '''
    Distrube a datainterval to the given number of processes.
    Return the working ranges in form of a list of tuples.
    Might be useful for more tasks.
    '''

    displ=[] # displacement
    displ.append(0)
    rct=[]   # receive count

    # calculate the count and displaceent by continuisly
    # dividing the remaining datarange by the remaining processes
    # until we reach the end
    for i in range(processes-1):
      rct.append((datainterval-displ[i])//(processes-i))
      displ.append(rct[i]+displ[i])
    rct.append(datainterval- displ[processes-1])

    # save the working ranges in form of a list of tuples
    # where the tuple contains the start and (not-included) stop k-point
    myranges = []
    for i in range(processes):
      myranges.append((displ[i],displ[i]+rct[i]))

    return myranges
