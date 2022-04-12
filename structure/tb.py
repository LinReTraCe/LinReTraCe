#!/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np
import spglib

from structure.aux   import levicivita
from structure.aux   import progressBar
from structure.model import Model
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

  def computeData(self, tbdata, rvec_spglib, atoms_spglib, charge, mu=None):
    self.tbdata       = tbdata
    self.rvec_spglib  = rvec_spglib # the vectors are the ROWS (required for spglib python interface)
    self.atoms_spglib = atoms_spglib
    self.charge       = charge

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
      self._checkSymmetriesTightbinding()
      self._checkSymmetriesKmesh()

    ''' after getting the correct number of k-points we can setup the arrays '''
    self._setupArrays(self.ortho)

    self._computeHk()
    self._calcOptical()

    ''' setting some important variables '''
    self.opticalBandMin = 0
    self.opticalBandMax = self.energyBandMax

    self._calcFermiLevel(mu)

  def _computeOrthogonality(self):
    '''
      weird way to determine orthogonality of unit cells
      that includes fcc and bcc vectors that are naturally not orthogonal
      but describe a conventionally orthogonal unit cell
      should work for all cases
    '''
    sum_vecs = np.sum(self.rvec_spglib, axis=0) # a_1 + a_2 + a_3
    max_vecs = np.array([np.max(np.abs(self.rvec_spglib[:,i])) for i in range(3)]) # maximal entries
    ratio = sum_vecs / max_vecs
    self.ortho = np.all(np.isclose(ratio, ratio[0]))
    logger.info('   orthogonal lattice: {}'.format(self.ortho))

  def _computeReciprocLattice(self):
    self.rvec = self.rvec_spglib
    self.vol = np.abs(np.dot(np.cross(self.rvec[:,0],self.rvec[:,1]),self.rvec[:,2]))
    self.kvec = np.zeros_like(self.rvec, dtype=np.float64)
    self.kvec[0,:] = np.cross(self.rvec[1,:],self.rvec[2,:]) / self.vol
    self.kvec[1,:] = np.cross(self.rvec[2,:],self.rvec[0,:]) / self.vol
    self.kvec[2,:] = np.cross(self.rvec[0,:],self.rvec[1,:]) / self.vol
    self.kvec *= 2*np.pi
    logger.info('   volume [Ang^3]: {}'.format(self.vol))
    logger.debug('   real space lattice (rows) :\n{}'.format(self.rvec))
    logger.debug('   reciprocal lattice (rows) :\n{}'.format(self.kvec))
    logger.debug('   recip.T @ real / (2pi) =\n{}'.format(self.kvec.T @ self.rvec / 2/np.pi))

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
    logger.debug(' Unique r-points:\n{}'.format(self.rpoints))

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

      lattice   = self.rvec_spglib
      positions = []
      numbers   = []
      for i in range(self.atoms_spglib.shape[0]):
        numbers.append(self.atoms_spglib[i,0])
        positions.append(self.atoms_spglib[i,1:])


      cell = (lattice, positions, numbers)

      ''' get spacegroup'''
      spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)
      logger.info('Spglib: Detected spacegroup {}'.format(spacegroup))

      logger.info('Spglib: Generating irreducible kpoints.')
      mapping, grid = spglib.get_ir_reciprocal_mesh(kgrid, cell, is_shift=is_shift)
      unique, counts = np.unique(mapping, return_counts=True)
      self.nkp  = len(unique)
      logger.info('Spglib: Generated irreducible kmesh with {} irreducible kpoints'.format(self.nkp))
      self.multiplicity = np.array(counts, dtype=int)
      self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))
      self.kpoints = grid[unique]
      self.kpoints = (self.kpoints + is_shift.astype(np.float64)/2.) / kgrid[None,:].astype(np.float64)

      logger.debug('kpoints:\n{}'.format(self.kpoints))
      ''' get symmetry and reduce unnecessary ones '''
      symmetry = spglib.get_symmetry(cell, symprec=1e-5)
      symsfull = symmetry['rotations']

      orthosym = True
      self.symop = []
      self.symopT = []
      for ii in np.arange(symsfull.shape[0]):
        isym = symsfull[ii]

        ''' Spglib gives a weird form of these rotations. this is an attempt to hotfix this mess.
            Each element in the point gropup _must_ have a determinant of either -1 or +1.
            Here we catch the problematic 'rotations' and overwrite the off-diagonal terms to
            produce the correct result.
            This is still very much experimental and I'm looking for a more robust solution.
        '''
        det = np.linalg.det(isym)
        if abs(det) != 0:
          isym[isym>1] = 0
          isym[isym<(-1)] = 0

        if abs(np.linalg.det(isym)) != 1.:
          logger.critical('Spglib: Throwing out symmetry operation without properly defined determinant:\n {}'.format(isym))

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

        if to_add:
          self.symop.append(isym)

      self.symop = np.array(self.symop)
      self.invsymop = np.linalg.inv(self.symop)

      self.nsym = self.symop.shape[0]
      logger.info('Spglib: Found {} symmetry operations'.format(self.nsym))
      logger.debug('Symmetry operation: \n{}'.format(self.symop))
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

    ''' FOURIERTRANSFORM hvk(j) = sum_r i . r_j e^{i r.k} * weight(r) * h(r) '''
    ''' self.rvec[i,j] is the jth entry of the ith vector (aka rows)
        -> to get Cartesian we need to dotproduct the unit cell displacement (self.rpoints)
        with the x/y/z entries (columns) of the rvec matrix
    '''
    prefactor_r = np.einsum('id,ri->dr', self.rvec, self.rpoints)
    hvk[:,:,:,:] = np.einsum('dr,kr,rij->kijd',1j*prefactor_r,ee,self.hr)

    ''' FOURIERTRANSFORM hvk(j) = sum_r - r_j e^{i r.k} * weight(r) * h(r) '''
    prefactor_r2 = np.zeros((6,self.nrp), dtype=np.float64)
    for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
      prefactor_r2[idir,:] = prefactor_r[i,:] * prefactor_r[j,:]
    hck[:,:,:,:] = np.einsum('dr,kr,rij->kijd',-prefactor_r2,ee,self.hr)


    if logger.isEnabledFor(logging.DEBUG):
      ''' irreducible point '''
      import random
      ik = random.randint(1,self.nkp-1) # avoid gamma point
      print('\n\n Randomly chosen k-points: {}'.format(ik))
      print('irreducible k: {}'.format(self.kpoints[ik,:]))
      print('irreducible hk: {}'.format(hk[ik,:,:]))
      print('irreducible hvk: {}'.format(hvk[ik,:,:,:]))
      print('irreducible hck: {}'.format(hck[ik,:,:,:]))
      print('multiplicity k: {}\n\n'.format(self.multiplicity[ik]))

      ''' generate all connected points via tranposed symmetry operations '''
      ''' on these explitily generated k-points .. apply the energy, velocity and curvature equations '''

      ''' yes this is the transposed symop here '''
      redk = np.einsum('nji,j->ni',self.symop,self.kpoints[ik])
      ''' bring back to BZ '''
      redk[redk<0] += 1
      redk[redk>1] -= 1
      red_rdotk = 2*np.pi*np.einsum('ni,ri->nr',redk,self.rpoints)
      red_ee = np.exp(1j * red_rdotk)
      #energy
      red_hk = np.einsum('nr,rij->nij', red_ee, self.hr)
      #velocitiy
      red_hvk = np.einsum('dr,nr,rij->nijd',1j*prefactor_r,red_ee,self.hr)
      #curvature
      red_hck = np.einsum('dr,nr,rij->nijd',-prefactor_r2,red_ee,self.hr)
      #mopt
      red_hvkvk = np.zeros((redk.shape[0],self.energyBandMax,self.energyBandMax,3,3), dtype=np.float64)
      red_hvkvk[:,:,:,[0,1,2,0,0,1], [0,1,2,1,2,2]] = (np.conjugate(red_hvk[:,:,:,[0,1,2,0,0,1]]) * red_hvk[:,:,:,[0,1,2,1,2,2]]).real
      red_hvkvk[:,:,:, [1,2,2], [0,0,1]] = red_hvkvk[:,:,:, [0,0,1], [1,2,2]]
      #curvature in matrix form
      red_hck_mat  = np.zeros((redk.shape[0],self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
      red_hck_mat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = red_hck[:,:,:,:]
      red_hck_mat[:,:,:, [1,2,2], [0,0,1]] = red_hck_mat[:,:,:,[0,0,1], [1,2,2]]

      ''' on the contrary apply the symmetry operations on the velocities of the irreducible point '''
      ''' apply the symmetry in matrix form onto curvature matrix of irreducible point '''
      testsymop = np.einsum('ij,njk,kl->nil',np.linalg.inv(self.kvec),self.invsymop,self.kvec)
      testsymopT = np.einsum('ij,njk,kl->nli',np.linalg.inv(self.kvec),self.invsymop,self.kvec)

      ''' transform into matrix form so we can apply the symmetries effectively '''
      hck_mat = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
      hck_mat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = hck[:,:,:,:]
      hck_mat[:,:,:, [1,2,2], [0,0,1]] = hck_mat[:,:,:, [0,0,1], [1,2,2]]
      hvkvk_mat = np.zeros((self.nkp,self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
      hvkvk_mat[:,:,:,[0,1,2,0,0,1], [0,1,2,1,2,2]] = np.conjugate(hvk[:,:,:,[0,1,2,0,0,1]]) * hvk[:,:,:,[0,1,2,1,2,2]]
      hvkvk_mat[:,:,:, [1,2,2], [0,0,1]] = hvkvk_mat[:,:,:, [0,0,1], [1,2,2]]

      symvk = np.einsum('nij,abj->nabi',testsymop,hvk[ik,:,:,:])
      symck = np.einsum('nij,abjk,nkl->nabil',testsymop,hck_mat[ik,:,:,:,:],testsymopT)
      symvkvk = np.einsum('nij,abjk,nkl->nabil',testsymop,hvkvk_mat[ik,:,:,:,:],testsymopT)

      print('irrk, ene(irrk) --- P^T irrk, ene(P^T irrk) # all band combinations')
      for i in range(redk.shape[0]):
        print(self.kpoints[ik,:3], '\n', hk[ik,:,:], '\n --- \n', redk[i,:3], '\n', red_hk[i,:,:])
        print()

      ''' debug output this only when it makes sense '''
      if self.energyBandMax == 1:
        print('\nTransformed vx vy vz(irrk) --- vx vy vz(P^T irrk)   # only real part of band 0')
        for i in range(redk.shape[0]):
          print(symvk[i,0,0,:].real, ' --- ', red_hvk[i,0,0,:].real)

        print('\nTransformed cxx cyy czz cxy cxz cyz (irrk) --- cxx cyy czz cxy cxz cyz (P^T irrk)   # only real part')
        for i in range(redk.shape[0]):
          print('[{} {} {} {} {} {}] --- [{} {} {} {} {} {}]'.format(symck[i,0,0,0,0].real, symck[i,0,0,1,1].real, symck[i,0,0,2,2].real, symck[i,0,0,0,1].real,symck[i,0,0,0,2].real,symck[i,0,0,1,2].real, \
                      red_hck_mat[i,0,0,0,0].real, red_hck_mat[i,0,0,1,1].real, red_hck_mat[i,0,0,2,2].real, red_hck_mat[i,0,0,0,1].real,red_hck_mat[i,0,0,0,2].real,red_hck_mat[i,0,0,1,2].real))

        print('\nTransformed vxx vyy vzz vxy vxz vyz(ikrr) --- vxx vyy vzz vxy vxz vyz (P^T irrk)   # only real part')
        for i in range(redk.shape[0]):
          print('[{} {} {} {} {} {}] --- [{} {} {} {} {} {}]'.format(symvkvk[i,0,0,0,0].real, symvkvk[i,0,0,1,1].real, symvkvk[i,0,0,2,2].real, symvkvk[i,0,0,0,1].real, symvkvk[i,0,0,0,2].real, symvkvk[i,0,0,1,2].real, \
                      red_hvkvk[i,0,0,0,0].real, red_hvkvk[i,0,0,1,1].real, red_hvkvk[i,0,0,2,2].real, red_hvkvk[i,0,0,0,1].real, red_hvkvk[i,0,0,0,2].real, red_hvkvk[i,0,0,1,2].real))



      ''' to get a full one-to-one comparison
          we need to sort these now
          through the application of rotations pairs of colums and rows can be interchanged
      '''
      ek,U = np.linalg.eig(hk[ik])      # 1 transformation
      ek = ek.real
      sortindex = ek.argsort()
      ek = ek[sortindex]
      U = U[:,sortindex]
      Uinv = np.linalg.inv(U)
      hk_sorted = hk[ik,sortindex,sortindex]

      symvk_sorted   = symvk.copy()
      symck_sorted   = symck.copy()
      symvkvk_sorted = symvkvk.copy()
      for i in range(redk.shape[0]):
        for idir in range(3):
          symvkdir = symvk[i,:,:,idir]
          symvk_sorted[i,:,:,idir] = symvkdir[sortindex,sortindex]
          for jdir in range(3):
            symckdir = symck[i,:,:,idir,jdir]
            symck_sorted[i,:,:,idir,jdir] = symckdir[sortindex,sortindex]
            symvkvkdir = symvkvk[i,:,:,idir,jdir]
            symvkvk_sorted[i,:,:,idir,jdir] = symvkvkdir[sortindex,sortindex]

      red_ek, red_U = np.linalg.eig(red_hk) # n transformations
      red_ek = red_ek.real
      red_hk_sorted    = red_hk.copy()
      red_hvk_sorted   = red_hvk.copy()
      red_hck_sorted   = red_hck_mat.copy()
      red_hvkvk_sorted = red_hvkvk.copy()
      for i in range(redk.shape[0]):
        idx = red_ek[i].argsort()
        red_ek[i,:] = red_ek[i,idx]
        red_hk_sorted[i,:,:] = red_hk[i,idx,idx]
        for idir in range(3):
          red_vkdir = red_hvk[i,:,:,idir]
          red_hvk_sorted[i,:,:,idir] = red_vkdir[idx,idx]
          for jdir in range(3):
            red_ckdir = red_hck_mat[i,:,:,idir,jdir]
            red_hck_sorted[i,:,:,idir,jdir] = red_ckdir[idx,idx]
            red_vkvkdir = red_hvkvk[i,:,:,idir,jdir]
            red_hvkvk_sorted[i,:,:,idir,jdir] = red_vkvkdir[idx,idx]

      ene_ev_check = np.allclose(ek, red_ek)
      ene_hk_check = np.allclose(hk_sorted, red_hk_sorted)
      velcheck = np.allclose(symvk_sorted, red_hvk_sorted)
      curcheck = np.allclose(symck_sorted, red_hck_sorted)
      optcheck = np.allclose(symvkvk_sorted, red_hvkvk_sorted)

      print('Eigenvalue   symmetry check (irreducible, reducible): {}'.format(ene_ev_check))
      print('Hamiltonian  symmetry check (irreducible, reducible): {}'.format(ene_hk_check))
      print('Velocitiy    symmetry check (rotated,     reducible): {}'.format(velcheck))
      print('Curvature    symmetry check (rotated,     reducible): {}'.format(curcheck))
      print('Optical      symmetry check (rotated,     reducible): {}'.format(optcheck))
      print('\n\n')



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

  def _checkSymmetriesTightbinding(self):

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
          logger.warning('\n\nTight binding symmetry check: {}'.format(transformed) + \
                        '\n    Symmetry of r-vector {} does not respect point group symmetries'.format(rvec1) + \
                        '\n    Break unit cell symmetry or avoid irreducible calculation if this is done on purpose.\n\n')
          return
    else:
      logger.info('Tight binding symmetry check: {}'.format(transformed))

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
      logger.warning('\n    Momentum mesh does not conform to point group symmetries.' + \
                     '\n    Change momentum grid or calculate on reducible grid (--red)\n\n')

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

      # we have to adjust the symmetry operations for the transformation to describe cartesian coordinates
      # i tried to mimic the way Wien2K does this. it works, but i have no idea why, good luck brave adventurer
      rotsymop  = np.einsum('ij,njk,kl->nil',np.linalg.inv(self.kvec),self.invsymop,self.kvec)
      rotsymopT = np.einsum('ij,njk,kl->nli',np.linalg.inv(self.kvec),self.invsymop,self.kvec)

      for ikp in range(self.nkp):
        progressBar(ikp+1,self.nkp,status='k-points')

        vel     = self.velocities[0][ikp,:,:,:] # bands, bands, 3
        cur     = self.curvatures[0][ikp,:,:,:] # bands, bands, 6

        # put the curvatures into matrix form
        curmat  = np.zeros((self.energyBandMax,self.energyBandMax,3,3), dtype=np.complex128)
        curmat[:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:]
        curmat[:,:, [1,2,2], [0,0,1]] = curmat[:,:, [0,0,1], [1,2,2]]

        # generate the transformed velocities and curvatures
        vk = np.einsum('nij,bpj->bpni',rotsymop,vel)
        ck = np.einsum('nij,bpjk,nkl->bpnil',rotsymop,curmat,rotsymopT) # bands, bands, nsym, 3, 3

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
      curmat[:,:,:, [1,2,2], [0,0,1]] = curmat[:,:,:, [0,0,1], [1,2,2]]
      vel2 = np.conjugate(vel[:,:,:,[0,1,2,0,0,1]]) * vel[:,:,:,[0,1,2,1,1,2]]
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

