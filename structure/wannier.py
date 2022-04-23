#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import glob
import sys
import logging
logger = logging.getLogger(__name__)

from structure.aux   import progressBar
from structure.dft   import DftCalculation
from structure.es    import ElectronicStructure
from structure.inout import h5output
from structure.aux   import levicivita
from structure.units import bohr2angstrom

import numpy as np

class Wannier90Calculation(DftCalculation):
  '''
  Wannier90 calculation class which reads all the relevant information
  hopping parameteres
  lattice vectors
  kmesh
  to _calculate_ the energy dispersion, band velocities and band curvatures
  this allows (after rotation into the Kohn Sham basis) to create
  intra + inter band optical elements ~ velocities**2
  intra + inter band magnetic optical elements ~ velocities**2 * curvature
  '''

  def __init__(self, directory, charge, **kwargs):
    logger.info("\nInitializing Wannier90 calculation.")
    super(Wannier90Calculation, self).__init__()
    self.directory = directory
    self.charge    = charge
    if isinstance(self.charge, float):
      if self.charge < 0: raise ValueError('Provided charge must be >= 0')

    if not os.path.exists(self.directory):
      raise IOError("Supplied path does not exist: " + self.directory)
    if not os.path.isdir(self.directory):
      raise IOError("Supplied path is not a directory: " + self.directory)

    self._defineCase()        # filename prefix
    self._defineFiles()       # define all possible hr file variants
    self._detectCalculation() # check which hr files we have -> calculation type
    self._checkFiles()        # check if all files are okay
    logger.info("Files successfully loaded.")

    ''' and some flags '''
    self.opticdiag  = True
    self.bopticdiag = True

    ''' define the usual spin resolved arrays '''
    self.energies        = []
    self.velocities      = []
    self.curvatures      = []
    self.opticalDiag     = []
    self.BopticalDiag    = []
    self.opticalMoments  = []
    self.BopticalMoments = []

  def readData(self):
    ''' get the lattice vectors, number of k-points
        number of projections, projection centers,
        hopping parameters '''
    try:
      self._readWout()
    except:
      logger.error('An error occured while reading case.wout*, defaulting to Angstrom')
      self.lengthscale = 1.0
    self._readNnkp()
    self._readHr()
    logger.info("Files successfully read.")

    if self.nproj > 1:
      self.opticfull  = True
      self.bopticfull = True

    # internal variables stemming from DFTcalculations
    self.energyBandMax  = self.nproj
    self.opticalBandMin = 0
    self.opticalBandMax = self.nproj

    if isinstance(self.charge, float):
      if self.charge > 2*self.nproj:
        raise IOError('Provided charge is larger than possible in projected bands.')

    minarr = np.array([np.min(self.kpoints[:,i]) for i in range(3)])
    if np.any(minarr > 0):
      self.kshift = True
    else:
      self.kshift = False
    maxarr = np.array([np.max(self.kpoints[:,i]) for i in range(3)])
    spacing = (1-maxarr+minarr) # minarr is > 0 for shifted grids
    self.nkx, self.nky, self.nkz = [int(np.around(1./i)) for i in spacing]
    logger.info("   Momentum Grid:  {} x {} x {}".format(self.nkx,self.nky,self.nkz))
    logger.info("   Momentum Shift: {}".format(self.kshift))

    if not (self.nkx*self.nky*self.nkz == self.nkp):
      raise IOError('Irreducible momentum grids not implemented')

    self._defineDimensions()
    logger.info("   Number of dimensions: {}".format(self.ndim))

    self.multiplicity = np.ones((self.nkp,), dtype=int)
    self.weights = self.multiplicity * self.weightsum / (self.nkx*self.nky*self.nkz)
    # wannier90 should always give us a reducible grid
    self.irreducible = False

  def expandKmesh(self, kmesh, kirr=False):
    '''
      Instead of the kmesh found in the Wannier90 calculation
      define a new one. also redefine multiplicities and weights
    '''

    logger.info('Overwriting kmesh information with provided one: {}'.format(kmesh))
    raise_error = False
    if self.nkx == self.nky and kmesh[0] != kmesh[1]:
        raise_error = True
    if self.nkx == self.nkz and kmesh[0] != kmesh[2]:
        raise_error = True
    if self.nky == self.nkz and kmesh[1] != kmesh[2]:
        raise_error = True

    for iknew, ikold in zip(kmesh, [self.nkx, self.nky, self.nkz]):
      if iknew < 1: raise_error = True
      if ikold == 1 and not (iknew == 1): raise_error = True

    if raise_error:
      raise IOError('Provided kmesh does not conform to original kmesh')

    self.nkx, self.nky, self.nkz = kmesh
    self.nkp = self.nkx*self.nky*self.nkz

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

    self.kpoints = []
    for ikx in self._kmeshx:
      for iky in self._kmeshy:
        for ikz in self._kmeshz:
          self.kpoints.append([ikx,iky,ikz])
    self.kpoints = np.array(self.kpoints)
    if self.kshift: self.kpoints += self._kmeshshift[None,:]

    self.multiplicity = np.ones((self.nkp,), dtype=int)
    self.weights = self.multiplicity * self.weightsum / (self.nkx*self.nky*self.nkz)

  def readWien2k(self):
    '''
      Read Wien2K struct and klist file as new kmesh
    '''

    self.fstruct = self.case+'.struct'
    self.fklist = self.case+'.klist'

    self._readWien2kStruct()
    self._readWien2kKlist()

  def _readWien2kKlist(self):
    '''
    Read the kpoint list.
    Read the multiplicity.
    Define the weights.
    '''

    logger.info("Wien2K reading: {}".format(self.fklist))
    with open(str(self.fklist), 'r') as klist:
      # first we get the divisor
      firstline = klist.readline()
      pos1 = firstline.find('(')+1
      pos2 = firstline.find(')')
      divisor = []
      divisor.append(int(firstline[pos1:pos1+3]))
      divisor.append(int(firstline[pos1+3:pos1+6]))
      divisor.append(int(firstline[pos1+6:pos1+9]))
      divisor = np.array(divisor, dtype=np.int)

      raise_error = False
      if self.nkx == self.nky and divisor[0] != divisor[1]:
          raise_error = True
      if self.nkx == self.nkz and divisor[0] != divisor[2]:
          raise_error = True
      if self.nky == self.nkz and divisor[1] != divisor[1]:
          raise_error = True

      for iknew, ikold in zip(divisor, [self.nkx, self.nky, self.nkz]):
        if iknew < 1: raise_error = True
        if ikold == 1 and not (iknew == 1): raise_error = True

      if raise_error:
        raise IOError('Provided kmesh does not conform to original kmesh')

      self.nkx, self.nky, self.nkz = divisor
      # determine dimension whether we find '1's in the divisor
      self.dims = np.logical_not(divisor == np.ones(3, dtype=np.int))
      self.ndim = 3 - np.sum(divisor == np.ones(3, dtype=np.int))

      # now we reset the file
      klist.seek(0)

      myklist = []
      mymult  = []
      for line in klist:
        if ("END" in line):
          break
        else:
          kx   = float(line[10:20])
          ky   = float(line[20:30])
          kz   = float(line[30:40])
          mult = int(float(line[50:55]))
          myklist.append([kx,ky,kz])
          mymult.append(mult)
      else: # if we did not find the break condition
        raise IOError('Wien2k {}: Did not find END statement.'.format(str(self.fklist)))

    self.kpoints = np.array(myklist, dtype=np.float64)
    self.multiplicity = np.array(mymult, dtype=int)
    self.weights = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))
    self.nkp = self.kpoints.shape[0]
    self.irreducible = not (self.nkx*self.nky*self.nkz == self.nkp)
    self.kpoints = self.kpoints / divisor[None,:]

    logger.info("   Wien2K number of dimensions: {}".format(self.ndim))
    logger.info("   Wien2K number of k-points: {}".format(self.nkp))

  def _readWien2kStruct(self):
    '''
    Read the number of inequivalent atoms.
    Necessary for the header of the energy files.
        Get the orthogonality of our system (check for 90 degree angles).
    '''

    logger.info("Wien2K Reading: {}".format(self.fstruct))
    with open(str(self.fstruct), 'r') as struct:
      # get the number of inequivalent atoms
      self.iatms = 0 # inequivalent atoms
      for line in struct:
        if line.startswith('ATOM'):
          self.iatms += 1

      if self.iatms == 0:
        raise IOError('Wien2k {}: Reading number of inequivalent atoms failed.'.format(str(self.fstruct)))
      else:
        logger.info("   Number of inequivalent atoms: {}".format(self.iatms))

      struct.seek(0)
      # get the number of symmetry operations
      for line in struct:
        if 'NUMBER OF SYMMETRY OPERATIONS' in line:
          self.nsym = int(line[:4])
          logger.info("   Number of symmetry operations: {}".format(self.nsym))
          break
      else:
        raise IOError('Wien2k {}: Reading number of symmetry operations failed.'.format(str(self.fstruct)))

      # without resetting we can continue reading the operations
      self.symop    = np.zeros((self.nsym,3,3), dtype=np.float64)
      self.invsymop = np.zeros((self.nsym,3,3), dtype=np.float64)
      for isym in range(self.nsym):
        for i in range(3):
          temp = struct.readline()
          self.symop[isym,i,:] = np.array([temp[:2],temp[2:4],temp[4:6]], dtype=np.float64)
        struct.readline() # skip
        self.invsymop[isym] = np.linalg.inv(self.symop[isym])

      if struct.readline() != "": # exactly at the EOF
        raise IOError('Wien2K {} is not at the EOF after reading'.format(str(self.fstruct)))

      logger.debug('Symmetry operations')
      logger.debug(self.symop)

  def diagData(self, peierlscorrection=True):
    ''' calculate e(r), v(r), c(r)
        we need to generate:
          self.energies[[nkp,nband]]
          self.velocties[[nkp,nband,3]]
          self.curvatures[[nkp,nband,6]]
          self.opticalDiag[[nkp,nband,3/6]]  3 if ortho
          self.BopticalDiag[[nkp,nband,3,3,3]]
          self.opticalMoments[[nkp,nband,nband,3/6]]  3 if ortho
          self.BopticalMoments[[nkp,nband,nband,3,3,3]] + write output function
    '''

    logger.info('Calculating and diagonalizing H(k), v(k), c(k)')
    for ispin in range(self.spins):

      ''' decorator for progress bar '''
      if self.spins == 1:
        prefix = ''
      else:
        if ispin == 0:
          prefix = 'up:'
        else:
          prefix = 'dn:'

      # hamiltonian H(r)
      hr = self.hrlist[ispin] # self.nrp, self.nproj, self.nproj
      # hamiltonian H(k)
      hk = np.zeros((self.nkp,self.nproj,self.nproj), dtype=np.complex128)
      # hamiltonian derivative d_dk H(k)
      # 3: x y z
      hvk = np.zeros((self.nkp,self.nproj,self.nproj,3), dtype=np.complex128)
      # hamiltonian double derivative d2_dk2 H(k)
      # 6: xx yy zz xy xz yz
      hck = np.zeros((self.nkp,self.nproj,self.nproj,6), dtype=np.complex128)

      ''' FORUIERTRANSFORM hk = sum_r e^{i r.k} * weight(r) * h(r) '''
      rdotk = 2*np.pi*np.einsum('ki,ri->kr',self.kpoints, self.rlist) # no scaling to recip vector / lattice vectors here ?
      ee = np.exp(1j * rdotk) / self.rmultiplicity[None,:] # rmultiplicity takes care of double counting issues
      hk[:,:,:] = np.einsum('kr,rij->kij', ee, hr)

      ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
      prefactor_r = np.einsum('id,ri->dr', self.rvec, self.rlist)
      hvk[:,:,:,:] = np.einsum('dr,kr,rij->kijd',1j*prefactor_r,ee,hr)

      ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
      prefactor_r2 = np.zeros((6,self.nrp), dtype=np.float64)
      for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
        prefactor_r2[idir,:] = prefactor_r[i,:] * prefactor_r[j,:]
      hck[:,:,:,:] = np.einsum('dr,kr,rij->kijd',-prefactor_r2,ee,hr)

      if peierlscorrection:
        # Jan's code snippet
        # generalized Peierls for multi-atomic unit-cells (and, obviously, supercells)
        distances = self.plist[:,None,:] - self.plist[None,:,:] # nproj nproj 3
        ri_minus_rj = np.einsum('id,abi->abd', self.rvec, distances)
        hvk_correction = 1j * hk[:,:,:,None] * ri_minus_rj[None,:,:,:]
        hvk += hvk_correction

      # eigk  : self.nkp, self.nproj
      # U     : self.nkp, self.nproj, self.nproj
      #       U[0, :, i] is the eigenvector corresponding to ek[0, i]
      # inv(U) @ hk @ U = ek

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
        logger.warn('Detected complex energies ... truncating')
      ''' apparently complex velocities and curvatures are allowed now '''
      # if np.any(np.abs(hvk.imag) > 1e-5):
      #   logger.warn('Detected complex velocities ... truncating')
      # if np.any(np.abs(hck.imag) > 1e-5):
      #   logger.warn('Detected complex curvatures ... truncating')

      ek = ek.real
      self.energies.append(ek)
      self.velocities.append(vk)
      self.curvatures.append(ck)

  def calcOptical(self):
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

    if self.irreducible: # requires symmetry operations
      for ispin in range(self.spins):
        if self.ortho:
          loc_opticalMoments = np.zeros((self.nkp,self.nproj,self.nproj,3), dtype=np.float64)
        else:
          loc_opticalMoments = np.zeros((self.nkp,self.nproj,self.nproj,9), dtype=np.float64)
        loc_BopticalMoments = np.zeros((self.nkp,self.nproj,self.nproj,3,3,3), dtype=np.complex128)

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
          vk_conj = np.conjugate(vk)
          ck = np.einsum('nij,bpjk,nkl->bpnil',rotsymop,curmat,rotsymopT) # bands, bands, nsym, 3, 3

          # take the mean over the squares
          vk2 = vk_conj[:,:,:,[0,1,2,0,0,1]] * vk[:,:,:,[0,1,2,1,2,2]]
          vk2 = np.mean(vk2,axis=2)

          if self.ortho:
            loc_opticalMoments[ikp,...] = vk2[...,:3].real
          else:
            loc_opticalMoments[ikp,:,:,:6] = vk2.real
            loc_opticalMoments[ikp,:,:,6:] = vk2[...,3:].imag

          #           epsilon_cij v_a v_j c_bi -> abc
          mb = np.einsum('zij,bpnx,bpnj,bpnyi->bpnxyz',levmatrix,vk_conj,vk,ck)
          mb = np.mean(mb,axis=2)
          loc_BopticalMoments[ikp,...] = mb


        self.opticalMoments.append(loc_opticalMoments)
        self.BopticalMoments.append(loc_BopticalMoments)
        self.opticalDiag.append(loc_opticalMoments[:,np.arange(self.nproj),np.arange(self.nproj),:])
        self.BopticalDiag.append(loc_BopticalMoments[:,np.arange(self.nproj),np.arange(self.nproj),...])

    else:
      for ispin in range(self.spins):
        vel = self.velocities[ispin] # nkp, nproj, nproj, 3
        vel_conj = np.conjugate(vel)
        cur = self.curvatures[ispin] # nkp, nproj, nproj, 6

        # transform into matrix form
        curmat  = np.zeros((self.nkp,self.nproj,self.nproj,3,3), dtype=np.complex128)
        curmat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:,:]
        curmat[:,:,:, [1,2,2], [0,0,1]] = curmat[:,:,:, [0,0,1], [1,2,2]]
        vel2 = vel_conj[:,:,:,[0,1,2,0,0,1]] * vel[:,:,:,[0,1,2,1,2,2]]
        if self.ortho:
          vel2 = vel2[:,:,:,:3].real
        else:
          temp = vel2.copy()
          vel2 = np.empty((self.nkp,self.nproj,self.nproj,9), dtype=np.float64)
          vel2[:,:,:,:6] = temp.real
          vel2[:,:,:,6:] = temp[:,:,:,:3].imag
        self.opticalMoments.append(vel2)
        vel2diag = vel2[:,np.arange(self.nproj),np.arange(self.nproj),:]
        self.opticalDiag.append(vel2diag)

          #           epsilon_cij v_a v_j c_bi -> abc
        mb = np.einsum('cij,knma,knmj,knmbi->knmabc',levmatrix,vel_conj,vel,curmat)
        self.BopticalMoments.append(mb)
        mbdiag = mb[:,np.arange(self.nproj),np.arange(self.nproj),:,:,:]
        self.BopticalDiag.append(mbdiag)
    if not self.ortho:
      truncate = True
      for ispin in range(self.spins):
        if np.any(np.abs(self.opticalMoments[...,6:]) > 1e-6):
          truncate = False
      if truncate:
        for ispin in range(self.spins):
          self.opticalMoments[ispin]  = self.opticalMoments[ispin][...,:6]
          self.opticalDiag[ispin]     = self.opticalDiag[ispin][...,:6]


  def outputData(self, fname, mu=None, intraonly=False):
    ''' call the output function '''
    if intraonly:
      self.opticfull = False
      self.bopticfull = False

    ''' we either get mu=None (with self.charge = number)
        or we get mu=0.0 with self.charge = None '''
    if mu is None:
      self._calcFermiLevel(mu)
    else:
      self.charge = 0
      self.charge_old = self._calcOccupation(mu)
      logger.warning('\n\nDetermined charge in the projection: {}'.format(self.charge_old))
      self.charge = np.around(self.charge_old)
      if self.charge != self.charge_old:
        logger.warning('Rounding to nearest integer: {}\nIf this behavior is not intended provide the desired charge with --charge\n\n'.format(self.charge))
        self._calcFermiLevel()
    h5output(fname, self, self)

  def _defineCase(self):
    '''
    Get the Path prefix for a Wannier90 calculation
    case = /path/to/folder/calc/calc
    the files then can be described with case.fileending
    custom case -> path/custom.filending
    '''

    nnkp = os.path.join(self.directory,'*.nnkp')
    files = glob.glob(nnkp)
    if len(files) >= 1:
      if len(files) > 1:
        logger.warn('Detected more than 1 nnkp file in provided folder: Choosing {}'.format(files[0]))
      purennkp = os.path.basename(files[0])
      temp = purennkp.split('.')  # abc.def.nnkp -> [abc,def,nnkp]
      self.case = '.'.join(temp[:-1]) # abc.def
    else:
      raise IOError('Could not detect Wannier90 calculation')
    self.case = os.path.join(self.directory,self.case)

  def _defineFiles(self):
    '''
      static file endings with varialbe file beginnings
    '''
    self.fhr     = self.case + '_hr.dat'
    self.fhrup   = self.case + '_hr.datup'
    self.fhrdn   = self.case + '_hr.datdn'

    self.fnnkp   = self.case + '.nnkp'

    self.fwout   = self.case + '.wout'
    self.fwoutup = self.case + '.woutup'
    self.fwoutdn = self.case + '.woutdn'

  def _detectCalculation(self):
    '''
      Depending on the existance of hr file variations
      detect the type of calculation
      currently restricted to polarized and unpolarized
      self.calctype:  1 - spin polarized
                      2 - spin unpolarized
    '''

    if os.path.isfile(self.fhrup) and os.stat(str(self.fhrup)).st_size != 0 \
       and os.path.isfile(self.fhrdn) and os.stat(str(self.fhrdn)).st_size != 0:
      self.calctype = 1
      logger.info("Detected spin polarized calculation.")
    elif os.path.isfile(self.fhr) and os.stat(str(self.fhr)).st_size != 0:
      self.calctype = 2
      self.weightsum = 2
      self.spins = 1
      logger.info("Detected spin unpolarized calculation.")
    else:
      raise IOError("Error: No matching hr file combinations found")

  def _checkFiles(self):
    '''
      Check that all necessary files are existing files and are not empty
      Also set the number of spins and weightsum accordingly.
    '''
    if not os.path.isfile(self.fnnkp):
      raise IOError('Error: case_hr.dat is missing')
    if os.stat(str(self.fnnkp)).st_size == 0:
      raise IOerror('Error: case_hr.dat is empty.')

    if self.calctype == 1:
      self.weightsum = 1
      self.spins = 2
      self.fhraccess = [self.fhrup, self.fhrdn]
      self.fwoutaccess = [self.fwoutup, self.fwoutdn]
    elif self.calctype == 2:
      self.weightsum = 2
      self.spins = 1
      self.fhraccess = [self.fhr]
      self.fwoutaccess = [self.fwout]

    for i in self.fhraccess:
      if not os.path.isfile(i):
        raise IOError('Error: case_hr.dat* is missing')
      if os.stat(i).st_size == 0:
        raise IOerror('Error: case_hr.dat* is empty.')

  def _readWout(self):
    '''
      get length scale
    '''

    self.lengthscale = []

    for fwout in self.fwoutaccess:
      logger.info("Reading: {}".format(fwout))
      with open(str(fwout), 'r') as wout:
        for line in wout:
          if 'Length Unit' in line:
            if 'Ang' in line:
              logger.info("   Length scale: Ang")
              self.lengthscale.append(1.0)
            else:
              logger.info("   Length scale: Bohr")
              self.lengthscale.append(bohr2angstrom)
            break
        else:
          raise IOError('Length scale not detected in case.wout* file')

    if len(self.lengthscale) == 2:
      if self.lengthscale[0] != self.lengthscale[1]:
        raise IOError('Up and Dn files have two different length scales')

    self.lengthscale = self.lengthscale[0]


  def _readNnkp(self):
    '''
      Retrieve all possible information from case.nnkp
      real lattice vector
      reciprocal lattice vector
      number of k-points
      array of k-points
      number of projections and the atomic position of each of them
    '''
    self.rvec = []
    logger.info("Reading: {}".format(self.fnnkp))
    with open(str(self.fnnkp),'r') as nnkp:
      while ( not nnkp.readline().startswith('begin real_lattice')):
        pass
      for i in range(3):
        line = nnkp.readline()
        rx, ry, rz = float(line[:12]), float(line[12:24]), float(line[24:36])
        self.rvec.append([rx,ry,rz])
      if not nnkp.readline().startswith('end real_lattice'):
        raise IOError('Wannier90 {} is not at the end of real_lattice after reading'.format(str(self.fnnkp)))
    self.rvec = np.array(self.rvec, dtype=np.float64)
    self.rvec *= self.lengthscale
    logger.debug('   real_lattice (rows): \n{}'.format(self.rvec))

    ''' if a_1 + a_2 + a_3 is a scaled vector of the maximal entries -> ortho '''
    sum_vecs = np.sum(self.rvec, axis=0) # a_1 + a_2 + a_3
    max_vecs = np.array([np.max(np.abs(self.rvec[:,i])) for i in range(3)]) # maximal entries
    ratio = sum_vecs / max_vecs
    self.ortho = np.all(np.isclose(ratio, ratio[0]))
    logger.info('   Orthogonal lattice: {}'.format(self.ortho))

    # V = (axb . c)
    self.vol = np.abs(np.dot(np.cross(self.rvec[0,:], self.rvec[1,:]),self.rvec[2,:]))
    self.vol *= self.lengthscale**3
    logger.info('   Deduced volume: {}'.format(self.vol))


    self.kvec = []
    with open(str(self.fnnkp),'r') as nnkp:
      while ( not nnkp.readline().startswith('begin recip_lattice')):
        pass
      for i in range(3):
        line = nnkp.readline()
        kx, ky, kz = float(line[:12]), float(line[12:24]), float(line[24:36])
        self.kvec.append([kx,ky,kz])
      if not nnkp.readline().startswith('end recip_lattice'):
        raise IOError('Wannier90 {} is not at the end of recip_lattice after reading'.format(str(self.fnnkp)))
    self.kvec = np.array(self.kvec, dtype=np.float64)
    self.kvec /= self.lengthscale
    logger.debug('   recip_lattice (rows): \n{}'.format(self.kvec))

    with open(str(self.fnnkp),'r') as nnkp:
      self.kpoints = []
      while ( not nnkp.readline().startswith('begin kpoints')):
        pass
      self.nkp = int(nnkp.readline())
      logger.info('   Number of kpoints: {}'.format(self.nkp))
      for i in range(self.nkp):
        line = nnkp.readline()
        kx, ky, kz = float(line[:14]), float(line[14:28]), float(line[28:42])
        self.kpoints.append([kx,ky,kz])
      if not nnkp.readline().startswith('end kpoints'):
        raise IOError('Wannier90 {} is not at the end of kpoints after reading'.format(str(self.fnnkp)))
    self.kpoints = np.array(self.kpoints, dtype=np.float64)

    with open(str(self.fnnkp),'r') as nnkp:
      while ( not nnkp.readline().startswith('begin projections')):
        pass
      self.nproj = int(nnkp.readline())
      logger.info('   Number of projections: {}'.format(self.nproj))
      self.plist = []
      for i in range(self.nproj):
        line = nnkp.readline()
        rx, ry, rz = float(line[:10]), float(line[10:21]), float(line[21:32])
        self.plist.append([rx,ry,rz])
        nnkp.readline() # skip the one after
      if not nnkp.readline().startswith('end projections'):
        raise IOError('Wannier90 {} is not at the end of projections after reading'.format(str(self.fnnkp)))
    self.plist = np.array(self.plist, dtype=np.float64)

  def _readHr(self):
    '''
      Retrieve all possible information from case_hr.dat*
      number of r-points
      r-point multiplicity
      and the full hopping matrix [nrp, projections, projections] complex128
    '''

    for filehr in self.fhraccess:
      logger.info("Reading: {}".format(filehr))
      with open(str(filehr),'r') as hr:
        hr.readline() # 'written on ...' line
        nproj = int(hr.readline())
        self.nrp = int(hr.readline())
        logger.info('   Number of rpoints: {}'.format(self.nrp))
        if (self.nproj != nproj):
          raise IOError('Inconsistent number of projections between case.nnkp anc case_hr.dat')

        self.rmultiplicity = []
        while (len(self.rmultiplicity) != self.nrp):
          line = hr.readline()
          mult = line.split() # I think this is save
          for imult in mult:
            self.rmultiplicity.append(imult)
        self.rmultiplicity = np.array(self.rmultiplicity, dtype=np.float64)

        self.hrlist = []
        self.rlist = []
        matrix = np.zeros((self.nrp,self.nproj,self.nproj), dtype=np.complex128)
        for ir in range(self.nrp):
          for ip in range(self.nproj**2):
            line = hr.readline()
            rx, ry, rz, p1, p2 = [int(line[0+i*5:5+i*5]) for i in range(5)]
            tr, tj = float(line[25:37]), float(line[37:])
            matrix[ir,p1-1,p2-1] = tr + 1j * tj
          self.rlist.append([rx,ry,rz])
        self.hrlist.append(matrix) # one spin .. so we can loop and append
        self.rlist = np.array(self.rlist, dtype=np.float64)

        if hr.readline() != "": # exactly at the EOF
          raise IOError('Wannier90 {} is not at the EOF after reading'.format(str(filehr)))



class wannier90hamiltonian(ElectronicStructure):
  '''
  Wannier90 H(k) file
  '''

  def __init__(self, hk_file, charge):
    '''
    Constructor from wannier90 hk file
    We enforce one spin here for the time being
    Calculation of intra/intra band optical elements are left for someone else (:
    '''
    raise NotImplementedError("Interface to Wannier90 H(k) Hamiltonianos not yet fully implemented")
    super(hamiltonian_wannier90, self).__init__()
    self.charge = charge

    ''' read hamiltonian and save parameters '''
    hk, klist, nbands = hamiltonian_wannier90.read_wannier90(hk_file)
    self.hk  = hk
    self.nkp = hk.shape[0]
    self.kpoints = klist
    self.energyBandMax = nbands

    ''' diagonalize hamiltonian and save the eigenvalues in hkdiag'''
    self.hkdiag, _ = np.linalg.eig(hk)

    ''' parameter for reducible grids and unpolarized spin calculations '''
    self.spins = 1
    self.vol   = 1
    self.multiplicity = np.ones((self.nkp,), dtype=np.int)
    self.weights = 2 * np.ones((self.nkp,), dtype=np.int) / float(self.nkp)
    self.weightsum = 2

    self.optic       = True
    self.opticalfull = False
    self.opticdiag   = True

    self.opticalBandMin = 0 # convention (see structure/dft.py) -- gets corrected in output
    self.opticalBandMax = self.energyBandMax

    # self.energies       = [self.hk[:,np.arange(self.energyBandMax),np.arange(self.energyBandMax)].real] # dummy
    self.energies       = [self.hkdiag[:,:].real] # should be correct
    self.velocities     = [np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64)]
    self.curvatures     = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalMoments = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalDiag    = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.BopticalDiag   = [np.zeros((self.nkp, self.energyBandMax, 3, 3, 3), dtype=np.float64)]

    self._calcFermiLevel()

  """
  Code adapted from https://github.com/w2dynamics/w2dynamics
  GNU General Public License Version 3
  """
  @staticmethod
  def read_wannier90(hk_file):
    """
    Reads a Hamiltonian f$ H_{bb'}(k) f$ from a text file.
    """

    with open(hk_file,'r') as hkdata:
      def nextline():
        line = hkdata.readline()
        return line[:line.find('#')].split()

      # parse header
      header = nextline()
      if header[0] == 'VERSION':
        warn("Version 2 headers are obsolete (specify in input file!)")
        nkpoints, natoms = map(int, nextline())
        lines = np.array([nextline() for _ in range(natoms)], np.int)
        nbands = np.sum(lines[:,:2])
        del lines, natoms
      elif len(header) != 3:
        warn("Version 1 headers are obsolete (specify in input file!)")
        header = list(map(int, header))
        nkpoints = header[0]
        nbands = header[1] * (header[2] + header[3])
      else:
        nkpoints, nbands, _ = map(int, header)
      del header

      hk = np.fromfile(hkdata, sep=" ") # shape: nkp nbands, nspins, nbands, nspins

    nspins = 1 # enforced here
    hk = hk.reshape(-1, 3 + 2 * nbands**2 * nspins**2)
    kpoints_file = hk.shape[0]
    if kpoints_file > nkpoints:
      print("truncating Hk points")
    elif kpoints_file < nkpoints:
      raise ValueError("problem! %d < %d" % (kpoints_file, nkpoints))
    kpoints = hk[:nkpoints, :3]

    hk = hk[:nkpoints, 3:].reshape(nkpoints, nbands, nbands, 2)
    hk = hk[...,0] + 1j * hk[...,1]
    print("Hamiltonian shape: {}".format(hk.shape))
    if not np.allclose(hk, hk.transpose(0,2,1).conj()):
      print("Hermiticity violation detected in Hk file")
    print("Hermiticity check of Hamiltonian successful.")

    return hk, kpoints, nbands

class hamiltonian_matrix(ElectronicStructure):
  '''
  Class to describe arbitrary square matrix Hamiltonian
  '''

  def __init__(self, hk, vk, ck, charge):
    '''
    Constructor from arbitrary input arrays
    We enforce one spin here for the time being
    Calculation of intra/intra band optical elements are left for someone else (:

    hk shape [nkp, nband, nband]
    vk shape [nkp, nband, nband, 3] ## x y z
    ck shape [nkp, nband, nband, 6] ## xx yy zz xy xz yz
    '''
    super(hamiltonian_wannier90, self).__init__()
    self.charge = charge
    self.nkp = hk.shape[0]
    self.energyBandMax = hk.shape[1]

    self.kpoints = klist

    self.spins = 1
    self.vol   = 1
    self.multiplicity = np.ones((self.nkp,), dtype=np.int)
    self.weights = 2 * np.ones((self.nkp,), dtype=np.int) / float(self.nkp)
    self.weightsum = 2

    self.optic       = True
    self.opticalfull = False # flag for inter band
    self.opticdiag   = True # flag for intra band

    self.opticalBandMin = 0 # convention (see structure/dft.py) -- gets corrected in output
    self.opticalBandMax = self.energyBandMax

    self.energies       = [hk[:,np.arange(nbands),np.arange(nbands)].real] # dummy
    self.velocities     = [np.zeros((self.nkp, self.energyBandMax, 3), dtype=np.float64)]
    self.curvatures     = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalMoments = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.opticalDiag    = [np.zeros((self.nkp, self.energyBandMax, 6), dtype=np.float64)]
    self.BopticalDiag   = [np.zeros((self.nkp, self.energyBandMax, 3, 3, 3), dtype=np.float64)]

    self._calcFermiLevel()

  def outputData(self, fname):
    '''
    Output the loaded data
    '''
    h5output(fname, self, self, peierls=True)
