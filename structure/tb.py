#!/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np
import spglib

from structure.auxiliary import levicivita
from structure.auxiliary import progressBar
from structure.model     import Model
import structure.symmetries.C1

class TightBinding(Model):
  '''
  Tightbinding class to create (ir-)reducible k-grids
  for arbitrary lattices.

  Main class that gets interfaces with ltb script
  This class generates irreducible grids only via spglib.
  Reducible grids do not require spglib.
  '''

  def __init__(self, nkx=1, nky=1, nkz=1, irreducible=True, kshift=False):
    super(TightBinding, self).__init__(nkx,nky,nkz)
    self.irreducible = irreducible  # generate irreducible grid instead of reducible
    self.kshift      = kshift       # shift by half a k-point to avoid Gamma point

    logger.info('Setting up tight binding with {} x {} x {} kpoints'.format(self.nkx,self.nky,self.nkz))

  def computeData(self, tbfile, charge, mu=None, mushift=False, corronly=False):
    self.tbfile       = tbfile
    self.charge       = charge
    self.corronly     = corronly

    self._readTb()
    self._computeOrthogonality() # sets self.ortho
    self._computeReciprocLattice() # defines reciprocal lattice vectors and volume

    if self.charge < 0 or self.charge > self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands : charge in [0,2*bands]')

    self._setupKmesh()
    if self.irreducible:
      self._checkSymmetriesTightbinding()
      self._checkSymmetriesKmesh()

    ''' after getting the correct number of k-points we can setup the arrays '''
    self._setupArrays(self.ortho)

    ''' calculate hamiltonian and optical elements '''
    self._computeHk()

    ''' setting some important variables '''
    self.opticalBandMin = 0
    self.opticalBandMax = self.energyBandMax

    self._calcFermiLevel(mu)
    if mushift:
      logger.info('Shifting energies -> Chemical Potential === 0.')
      self.energies[0][...] -= self.mu
      # dont forget about possible valence and conduction energies ... np.nan + value = np.nan
      self.ecb[0] -= self.mu
      self.evb[0] -= self.mu
      self.mu = 0.0

  def _computeOrthogonality(self):
    '''
      weird way to determine orthogonality of unit cells
      that includes fcc and bcc vectors that are naturally not orthogonal
      but describe a conventionally orthogonal unit cell
      should work for all cases
    '''
    sum_vecs = np.sum(self.rvec, axis=0) # a_1 + a_2 + a_3
    max_vecs = np.array([np.max(np.abs(self.rvec[:,i])) for i in range(3)]) # maximal entries
    ratio = sum_vecs / max_vecs
    self.ortho = np.all(np.isclose(ratio, ratio[0]))
    logger.info('   orthogonal lattice: {}'.format(self.ortho))

  def _computeReciprocLattice(self):
    self.vol = np.abs(np.dot(np.cross(self.rvec[0,:],self.rvec[1,:]),self.rvec[2,:]))
    self.kvec = np.zeros_like(self.rvec, dtype=np.float64)
    self.kvec[0,:] = np.cross(self.rvec[1,:],self.rvec[2,:]) / self.vol
    self.kvec[1,:] = np.cross(self.rvec[2,:],self.rvec[0,:]) / self.vol
    self.kvec[2,:] = np.cross(self.rvec[0,:],self.rvec[1,:]) / self.vol
    self.kvec *= 2*np.pi
    logger.info('   volume [Ang^3]: {}'.format(self.vol))
    logger.debug('   real space lattice (rows) :\n{}'.format(self.rvec))
    logger.debug('   reciprocal lattice (rows) :\n{}'.format(self.kvec))
    logger.debug('   recip.T @ real / (2pi) =\n{}'.format(self.kvec.T @ self.rvec / 2/np.pi))

  def _readTb(self):
    '''
      read tight binding input
      format is akin to wannier90
    '''

    logger.info("Reading: {}".format(self.tbfile))

    ''' read atom block '''
    self.atoms = []
    with open(str(self.tbfile),'r') as data:
      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if line.startswith('begin atoms'): break
        else:
          raise IOError('ltb could not detect "begin atoms" block')

      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if len(line)==0 or line.startswith('#'): continue
          if line.startswith('end atoms'): break
        else:
          raise IOError('ltb is not at the "end atoms" marker after reading'.format(str(self.tbfile)))

        comment = line.find('#')
        if comment != -1: line = line[:comment]
        atom_id, rx, ry, rz = [float(i) for i in line.split()]
        if abs(np.rint(atom_id) - atom_id) > 1e-6: raise IOError('Invalid atom description: {}'.format(line))
        if rx < 0 or rx > 1: raise IOError('Invalid atom description: {}'.format(line))
        if ry < 0 or ry > 1: raise IOError('Invalid atom description: {}'.format(line))
        if rz < 0 or rz > 1: raise IOError('Invalid atom description: {}'.format(line))
        self.atoms.append([atom_id, rx, ry, rz])
    self.atoms = np.array(self.atoms, dtype=np.float64)
    logger.debug('   atoms: \n{}'.format(self.atoms))

    ''' read real_lattice block '''
    self.rvec = []
    with open(str(self.tbfile),'r') as data:
      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if line.startswith('begin real_lattice'): break
        else:
          raise IOError('ltb could not detect "begin real_lattice" block')

      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if len(line)==0 or line.startswith('#'): continue
          if line.startswith('end real_lattice'): break
          comment = line.find('#')
          if comment != -1: line = line[:comment]
          rx, ry, rz = [float(i) for i in line.split()]
          self.rvec.append([rx,ry,rz])
        else:
          raise IOError('ltb is not at the "end real_lattice" marker after reading'.format(str(self.tbfile)))

    self.rvec = np.array(self.rvec, dtype=np.float64)
    if self.rvec.shape != (3,3):
      raise IOError('Provided real_lattice is not a properly defined 3x3 matrix')
    logger.debug('   real_lattice (rows): \n{}'.format(self.rvec))

    ''' read hopping block '''
    self.tbdata = []
    bandmin = bandmax = None
    with open(str(self.tbfile),'r') as data:
      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if line.startswith('begin hopping'): break
        else:
          raise IOError('ltb could not detect "begin hopping" block')

      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if len(line)==0 or line.startswith('#'): continue
          if line.startswith('end hopping'): break

          comment = line.find('#')
          if comment != -1: line = line[:comment]
          linedata = line.split()
          lenline  = len(linedata)
          if lenline==6 or lenline==7:
            self.tbdata.append(np.array([float(i) for i in linedata]))
          else:
            raise IOError('Provided hopping paramater line is invald: {}'.format(line))
        else:
          raise IOError('ltb is not at the "end hopping" marker after reading'.format(str(self.tbdata)))

    ''' we do not transform the list of lists into an array since the different lines
        could have different lengths (imaginary part optional) '''

    ''' detect inter band transitions and minimum and maximum band entry '''
    inter     = False
    bandcheck = set()
    for i in range(len(self.tbdata)):
      locmin = min(self.tbdata[i][3:5].astype(int))
      locmax = max(self.tbdata[i][3:5].astype(int))
      if bandmin is None:
        bandmin = locmin
        bandmax = locmax
      else:
        if locmin < bandmin: bandmin = locmin
        if locmax > bandmax: bandmax = locmax
      if locmin != locmax: inter = True
      bandcheck.update({locmin,locmax})
    if bandmin != 1:
      raise IOError('Tight binding parameter set must start at band 1')
    if len(bandcheck) != bandmax:
      raise IOError('Tight binding parameter set does not contain the full band range')

    if inter:
      self.opticfull  = True
      self.bopticfull = True
      logger.info('   Detected inter-band hopping parameters.')
    else:
      logger.info('   Detected only intra-band hopping parameters.')

    self.energyBandMax  = bandmax
    logger.info('   Number of bands: {}'.format(self.energyBandMax))


    ''' read optional orbitals block '''
    self.orbitals = []
    orbitals_exist = False
    with open(str(self.tbfile),'r') as data:
      while True:
        line = data.readline()
        if line:
          line = line.strip()
          if line.startswith('begin orbitals'):
            orbitals_exist = True
            break
        else:
          orbitals_exist = False # optional
          break

      if orbitals_exist:
        while True:
          line = data.readline()
          if line:
            line = line.strip()
            if len(line)==0 or line.startswith('#'): continue
            if line.startswith('end orbitals'): break

            comment = line.find('#')
            if comment != -1: line = line[:comment]
            orbital_id, rx, ry, rz = [float(i) for i in line.split()]
            if abs(np.rint(orbital_id) - orbital_id) > 1e-6: raise IOError('Invalid atom description: {}'.format(line))
            if rx < 0 or rx > 1: raise IOError('Invalid atom description: {}'.format(line))
            if ry < 0 or ry > 1: raise IOError('Invalid atom description: {}'.format(line))
            if rz < 0 or rz > 1: raise IOError('Invalid atom description: {}'.format(line))
            orbital_id = int(orbital_id)
            self.orbitals.append([orbital_id, rx, ry, rz])
          else:
            raise IOError('ltb is not at the "end orbitals" marker after reading'.format(str(self.tbdata)))
      else:
        self.orbitals = None


    ''' validate orbital list '''
    if orbitals_exist:
      orbital_set = set()
      for orbital in self.orbitals:
        orbital_id = orbital[0]
        if orbital_id < 1 or orbital_id > self.energyBandMax:
          raise IOError('Provided orbital list contains orbital outside provided hopping range')
        orbital_set.add(orbital_id)
      if len(orbital_set) != self.energyBandMax:
        raise IOError('Provided orbital list does not contain all necessary orbitals.')

      ''' rewrite into ordered numpy array '''
      orbitals = np.zeros((self.energyBandMax,3), dtype=np.float64)
      for orbital in self.orbitals:
        orbital_id, rx, ry, rz = orbital
        orbitals[orbital_id-1,:] = rx, ry, rz
      self.orbitals = orbitals

    if self.orbitals is not None:
      logger.info('   Detected orbital positions.')
    logger.debug('orbital list: \n{}'.format(self.orbitals))

    hrset = set()
    rset  = set()
    for i in range(len(self.tbdata)):
      rx, ry, rz, band1, band2 = self.tbdata[i][:5].astype(int)
      rset.add((rx,ry,rz))
      hrset.add((rx,ry,rz,band1,band2))

    if len(hrset) != len(self.tbdata):
      raise IOError('Detected duplicate entries in the hopping parameters')

    self.nrp = len(rset)
    logger.info('   Number of unique r-points: {}'.format(self.nrp))
    self.rpoints = np.array(list(rset), dtype=int)
    logger.debug(' Unique r-points:\n{}'.format(self.rpoints))

    self.hr = np.zeros((self.nrp,self.energyBandMax,self.energyBandMax), dtype=np.complex128)
    for i in range(len(self.tbdata)):
      rvec = self.tbdata[i][:3].astype(int)
      orb1, orb2 = self.tbdata[i][3:5].astype(int)
      hop = self.tbdata[i][5]
      try:
        hop += 1j * self.tbdata[i][6]
      except:
        pass

      ''' re-interpret sign of local energies '''
      ir = np.argwhere(np.sum(np.abs(self.rpoints-rvec[None,:]),axis=1)==0)[0][0]
      if orb1==orb2:
        if hop.imag:
          raise IOError('\nDetected complex intra-band values \n{} in line \n{}\nFix input file.'.format(hop,self.tbdata[i]))
      if np.all(rvec==np.zeros((3,), dtype=int)) and orb1==orb2:
        hop *= (-1)
      self.hr[ir,orb1-1,orb2-1] += -hop

    ''' check sparsity of h(r) matrix '''
    zeros = None
    for i in np.arange(self.nrp):
      if zeros is None:
        zeros = (self.hr[i] == 0)
      else:
        zeros = np.logical_or(zeros, (self.hr[i] == 0))
    number_zeros = zeros.sum()
    self.sparse = number_zeros / self.energyBandMax**2

    logger.debug('Sparsity of H(R) matrix: {0:.2f} %'.format(self.sparse * 100))
    logger.debug('Tight binding parameter set:\n{}'.format(self.hr))

  def _setupKmesh(self):
    '''
      employ spglib to create irreducible kpoints,
      multiplicity and symmetry operations
      otherwise create the usual reducible cell
    '''

    kgrid = np.array([self.nkx,self.nky,self.nkz], dtype=int)
    if self.kshift:
      is_shift = np.array([int(i) for i in self.dims], dtype=np.float64)
    else:
      is_shift = np.array([0,0,0], dtype=int)

    ''' define k-grid, shift if required '''
    if self.irreducible:

      lattice   = self.rvec
      positions = []
      numbers   = []
      for i in range(self.atoms.shape[0]):
        numbers.append(self.atoms[i,0])
        positions.append(self.atoms[i,1:])


      cell = (lattice, positions, numbers)

      ''' get spacegroup'''
      spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)
      logger.info('Spglib: Detected spacegroup {}'.format(spacegroup))

      logger.info('   Generating irreducible kpoints.')
      mapping, grid = spglib.get_ir_reciprocal_mesh(kgrid, cell, is_shift=is_shift)
      unique, counts = np.unique(mapping, return_counts=True)
      self.nkp  = len(unique)
      logger.info('   Generated irreducible kmesh with {} irreducible kpoints.'.format(self.nkp))
      self.multiplicity = np.array(counts, dtype=int)
      self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))
      self.kpoints = grid[unique]
      self.kpoints = (self.kpoints + is_shift.astype(np.float64)/2.) / kgrid[None,:].astype(np.float64)

      logger.debug('kpoints:\n{}'.format(self.kpoints))
      ''' get symmetry and reduce unnecessary ones '''
      # try:
      #   lattice, scaled_positions, numbers = spglib.standardize_cell(cell, symprec=1e-5)
      #   cell_standardized = (lattice, scaled_positions, numbers)
      # except:
      #   cell_standardized = cell
      symmetry = spglib.get_symmetry(cell, symprec=1e-5)
      symsfull = symmetry['rotations']

      non_standard = False
      self.symop = []
      for ii in np.arange(symsfull.shape[0]):
        isym = symsfull[ii]

        ''' Truncate if we detect a non-standardized unit cell.
            I avoid redefining the unit cell, so the user knows what is happening
            If the truncation results in a matrix without proper determinant raise an Exception.'''
        if abs(np.linalg.det(isym)) != 1:
          non_standard = True
          isym[isym>1] = 0
          isym[isym<(-1)] = 0
          if abs(np.linalg.det(isym)) != 1:
            raise ValueError('Non-stanardized unit cell resulted in matrix with invalid determinant')

        ''' Filter out the symmetries corresponding to dimensions not in use (deactivated via nk_i = 1) '''
        to_add = True
        for i, dim in enumerate(self.dims):
          if dim: continue
          for j in range(3):
            if i==j:
              if isym[i,j] != 1: to_add = False
            else:
              if isym[j,i] != 0: to_add = False
              if isym[i,j] != 0: to_add = False # redundant I think

        ''' avoid duplicates '''
        for jsym in self.symop:
          if np.allclose(isym,jsym):
            to_add = False

        if to_add:
          self.symop.append(isym)

      self.symop = np.array(self.symop)
      self.invsymop = np.linalg.inv(self.symop)

      self.nsym = self.symop.shape[0]
      if non_standard:
        logger.critical('\n\n   Detected non-standardized unit cell. Continuing by truncating rotations.' + \
                        '\n   Validity of generated data can not be guaranteed.\n')
      logger.info('   Found {} applicable symmetry operations.'.format(self.nsym))
      logger.debug('Symmetry operations: \n{}'.format(self.symop))
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

      logger.info('Reducible grid: {} symmetry operation'.format(self.nsym))

  def _computeHk(self):
    '''
      calculate H(k) d/dk H(k) and d2/dk2 H(k)
      via H(r) and primitive lattice vectors
    '''

    logger.info('Calculating Hamiltonian and Optical elements')

    ''' r vector for hv(k) in Cartesian coordinates'''
    prefactor_r = np.einsum('id,ri->dr', self.rvec, self.rpoints)
    prefactor_r2 = np.zeros((6,self.nrp), dtype=np.float64)

    ''' r.r matrix for hc(k) in Cartesian coordinates -- we only need the upper triangle '''
    for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
      prefactor_r2[idir,:] = prefactor_r[i,:] * prefactor_r[j,:]

    # hamiltonian H(k)
    hk = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax), dtype=np.complex128)
    # hamiltonian derivative d_dk H(k)
    # 3: x y z
    hvk = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3), dtype=np.complex128)
    # hamiltonian double derivative d2_dk2 H(k)
    # 6: xx yy zz xy xz yz
    hck = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,6), dtype=np.complex128)

    # devine levicity matrix so we can use it in einsteinsummations
    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)


    if self.irreducible:
      '''
      as we deal with a general multi-orbital case here
      we cannot apply velocity / curvature rotations
      as they are only properly defined in the band basis

      in order to perform the symmetrization, we generate the reducible points
      via symmetry operations and calculate the velocities / curvatures
      on the generated points explicitly

      the symmetrized optical elements are then the average over all these points
      '''

      ''' setup arrays '''
      if self.ortho:
        loc_opticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3), dtype=np.float64)
      else:
        loc_opticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,9), dtype=np.float64)
      loc_BopticalMoments = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3,3), dtype=np.complex128)

      for ikp in range(self.nkp):
        progressBar(ikp+1,self.nkp,status='k-points')

        ''' generate reducible k-points and bring back to first BZ'''
        redk = np.einsum('nji,j->ni',self.symop,self.kpoints[ikp]) # k_red = P^T . k_irr
        redk = redk%1

        ''' generate hamiltonian '''
        red_rdotk = 2*np.pi*np.einsum('ni,ri->nr',redk,self.rpoints)
        red_ee = np.exp(1j * red_rdotk)
        red_hk = np.einsum('nr,rij->nij', red_ee, self.hr)
        red_hk[np.abs(red_hk) < 1e-14] = 0.0

        ''' generate velocity hamiltonian '''
        red_hvk = np.einsum('dr,nr,rij->nijd',1j*prefactor_r,red_ee,self.hr)

        if self.orbitals is not None:
          # Jan's code snippet
          # generalized Peierls for multi-atomic unit-cells (and, obviously, supercells)
          # correction term: -1j (rho_l^alpha - rho_l'^alpha) H_ll' (k)
          distances = self.orbitals[:,None,:] - self.orbitals[None,:,:] # nproj nproj 3 -- w.r.t basis vectors
          ri_minus_rj = np.einsum('id,abi->abd', self.rvec, distances) # cartesian directions
          hvk_correction = - 1j * red_hk[:,:,:,None] * ri_minus_rj[None,:,:,:]
          red_hvk += hvk_correction


        ''' generate curvature hamiltonian '''
        red_hck = np.einsum('dr,nr,rij->nijd',-prefactor_r2,red_ee,self.hr)
        red_hck_mat  = np.zeros((self.nsym,self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
        red_hck_mat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = red_hck[:,:,:,:]
        red_hck_mat[:,:,:, [1,2,2], [0,0,1]] = red_hck_mat[:,:,:,[0,0,1], [1,2,2]]


        ''' imporant to note: the generated hamiltonians on the reducible grid
            must not be identical to the hamiltonian of the irreducible point
            what must be identical is only the eigenvalue
        '''

        ''' to get a full one-to-one comparison
            we need to sort these now
            through the application of rotations pairs of colums and rows can be interchanged '''
        red_ek, red_U = np.linalg.eigh(red_hk)
        red_ek = red_ek.real

        for isym in range(self.nsym):
          ekk, Uk = red_ek[isym,:], red_U[isym,:,:]
          idx = ekk.argsort()
          red_ek[isym,:] = ekk[idx]
          red_U[isym,:,:] = Uk[:,idx]

        self.energies[0][ikp,:] = red_ek[0,:]
        red_Uinv = np.linalg.inv(red_U)

        ''' transform velocity hamiltonian into Kohn Sham basis '''
        vel = np.einsum('nij,njkd,nkl->nild',red_Uinv,red_hvk,red_U)
        velconj = np.conjugate(vel)
        ''' transform curvature hamiltonian into Kohn Sham basis '''
        curmat = np.einsum('nij,njkab,nkl->nilab',red_Uinv,red_hck_mat,red_U)

        ''' take the mean over the opt matrix (velocity squares) '''
        vk2 = velconj[:,:,:,[0,1,2,0,0,1]] * vel[:,:,:,[0,1,2,1,2,2]]
        vk2 = np.mean(vk2,axis=0)

        if self.ortho:
          loc_opticalMoments[ikp,...] = vk2[...,:3].real
        else:
          loc_opticalMoments[ikp,:,:,:6] = vk2.real
          loc_opticalMoments[ikp,:,:,6:] = vk2[...,3:].imag

        ''' take the mean over the optb matrix '''
        #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('zij,bpnx,bpnj,bpnyi->bpnxyz',levmatrix,velconj,vel,curmat)
        mb = np.mean(mb,axis=0)
        loc_BopticalMoments[ikp,...] = mb

      self.opticalMoments[0][...]  = loc_opticalMoments
      self.opticalDiag[0][...]     = loc_opticalMoments[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:]
      self.BopticalMoments[0][...] = loc_BopticalMoments
      self.BopticalDiag[0][...]    = loc_BopticalMoments[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),...]

      if not self.ortho:
        ''' truncate if there are only negligible imaginary terms '''
        if np.all(np.abs(self.opticalMoments[0][...,6:]) < 1e-6):
          self.opticalMoments[0]  = self.opticalMoments[0][...,:6]
          self.opticalDiag[0]  = self.opticalDiag[0][...,:6]

    else:

      ''' FOURIERTRANSFORM h(k) '''
      rdotk = 2*np.pi*np.einsum('ki,ri->kr',self.kpoints,self.rpoints)
      ee = np.exp(1j * rdotk)
      hk[:,:,:] = np.einsum('kr,rij->kij', ee, self.hr)

      ''' this solves the problem of almost singular matrices
          the eigenvalue routine then ignores a full 0 matrix -> thats what we want '''
      hk[np.abs(hk) < 1e-14] = 0.0

      ''' FOURIERTRANSFORM hvk(alpha) = sum_r complex_i . r_alpha . e^{i r.k} * weight(r) * h(r) '''
      ''' self.rvec[i,j] is the jth entry of the ith vector (aka rows)
          -> to get Cartesian we need to dotproduct the unit cell displacement (self.rpoints)
          with the x/y/z entries (columns) of the rvec matrix
      '''
      hvk[:,:,:,:] = np.einsum('dr,kr,rij->kijd',1j*prefactor_r,ee,self.hr)

      ''' FOURIERTRANSFORM hvk(alpha,beta) = - sum_r r_alpha . r_beta . e^{i r.k} * weight(r) * h(r) '''
      hck[:,:,:,:] = np.einsum('dr,kr,rij->kijd',-prefactor_r2,ee,self.hr)

      if self.corronly:
        hvk[...] = 0.0

      if self.orbitals is not None:
        # Jan's code snippet
        # generalized Peierls for multi-atomic unit-cells (and, obviously, supercells)
        # correction term: -1j (rho_l^alpha - rho_l'^alpha) H_ll' (k)
        distances = self.orbitals[:,None,:] - self.orbitals[None,:,:] # nproj nproj 3 -- w.r.t basis vectors
        ri_minus_rj = np.einsum('id,abi->abd', self.rvec, distances) # cartesian directions
        hvk_correction = - 1j * hk[:,:,:,None] * ri_minus_rj[None,:,:,:]
        hvk += hvk_correction


      ''' this transforms all k points at once '''
      ek, U = np.linalg.eigh(hk)

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

      vel = np.einsum('kab,kbci,kcd->kadi',Uinv,hvk,U)
      vel_conj = np.conjugate(vel)
      vel2 = vel_conj[:,:,:,[0,1,2,0,0,1]] * vel[:,:,:,[0,1,2,1,2,2]]

      cur = np.einsum('kab,kbci,kcd->kadi',Uinv,hck,U)
      # transform into matrix form
      curmat  = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
      curmat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:,:]
      curmat[:,:,:, [1,2,2], [0,0,1]] = curmat[:,:,:, [0,0,1], [1,2,2]]

      if np.any(np.abs(ek.imag) > 1e-5):
        logger.warning('\n\nDetected complex energies ... truncating\n')
      self.energies[0][...] = ek.real

      if self.ortho:
        vel2 = vel2[:,:,:,:3].real
      else:
        if np.any(np.abs(vel2.imag) > 1e-6):
          temp = vel2.copy()
          vel2 = np.empty((self.nkp,self.energyBandMax,self.energyBandMax,9), dtype=np.float64)
          vel2[:,:,:,:6] = temp.real
          vel2[:,:,:,6:] = temp[:,:,:,3:].imag
        else:
          vel2 = vel2.real
      self.opticalMoments = [vel2]
      vel2diag            = vel2[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:]
      self.opticalDiag    = [vel2diag]

      #           epsilon_cij v_a v_j c_bi -> abc
      mb = np.einsum('cij,knma,knmj,knmbi->knmabc',levmatrix,vel_conj,vel,curmat)
      self.BopticalMoments[0][...] = mb
      mbdiag                       = mb[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax),:,:,:]
      self.BopticalDiag[0][...]    = mbdiag

    if not self.irreducible and logging.getLogger().isEnabledFor(logging.DEBUG):
      self.hk           = hk
      self.hvk          = hvk
      self.hck          = hck
      self.Ukohnsham    = U
      self.Uinvkohnsham = Uinv

  def _checkSymmetriesTightbinding(self):

    for itb1 in range(len(self.tbdata)):
      rvec1  = self.tbdata[itb1][:3].astype(int)
      band1_1  = int(self.tbdata[itb1][3]) - 1 # band identifier
      band1_2  = int(self.tbdata[itb1][4]) - 1 # band identifier
      hop1   = self.tbdata[itb1][5]
      try:
        hop1 += 1j * self.tbdata[itb1][6]
      except:
        pass

      ''' Inter-atomic hopping cannot be safely check: skip it '''
      if band1_1 != band1_2: continue

      rvecsym = np.einsum('nij,j->ni',self.symop,rvec1)

      for isym in range(self.nsym):
        rvec_transformed = rvecsym[isym]

        for itb2 in range(len(self.tbdata)):
          band2_1  = int(self.tbdata[itb1][3]) - 1 # band identifier
          band2_2  = int(self.tbdata[itb1][4]) - 1 # band identifier
          if band1_1 != band2_1 or band2_2 != band2_2: continue

          rvec2  = self.tbdata[itb2][:3].astype(int)
          hop2   = self.tbdata[itb2][5]
          try:
            hop2 += 1j * self.tbdata[itb2][6]
          except:
            pass

          if np.allclose(rvec_transformed,rvec2) and np.abs(hop1-hop2) < 1e-6:
            break
        else:
          logger.warning('\n\nIntra-orbital tight binding symmetry check: False' + \
                        '\n    Symmetry of r-vector {} does not respect point group symmetries'.format(rvec1) + \
                        '\n    Break unit cell symmetry or avoid irreducible calculation if this is done on purpose.\n\n')
          return
    else:
      logger.info('Intra-orbital tight binding symmetry check: True')

  def _checkSymmetriesKmesh(self):
    '''
      compare the momentum mesh to the unit cell symmetries
    '''
    nkvec = 1./np.array([self.nkx,self.nky,self.nkz], dtype=np.float64) # smallest possible k-spacing
    transformed = np.abs(np.einsum('nij,j->ni',self.invsymop,nkvec))
    spacing_exact = transformed / nkvec[None,:]
    spacing_round = np.rint(spacing_exact)
    conform = np.all(np.isclose(spacing_round,spacing_exact))
    logger.info('Momentum grid symmetry check: {}'.format(str(conform)))
    if not conform:
      logger.critical('\n    Momentum mesh does not conform to point group symmetries.' + \
                      '\n    Change momentum grid or calculate on reducible grid (--red)\n\n')
