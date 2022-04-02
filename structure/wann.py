#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import glob
import sys
import logging
logger = logging.getLogger(__name__)

from structure.aux   import progressBar
from structure.dft   import DFTcalculation
from structure.es    import ElectronicStructure
from structure.inout import h5output
from structure.aux   import levicivita

import scipy.optimize
import numpy as np

try:
  import ase
  import ase.io
  import ase.spacegroup
  from ase import Atoms
  from ase.spacegroup import get_spacegroup
  mymethods = {'get_volume', 'get_global_number_of_atoms'}
  Atommethods = set([method for method in dir(Atoms) if callable(getattr(Atoms, method))])
  ase_exists = True if mymethods.issubset(Atommethods) else False
except ImportError:
  ase_exists = False

class wannier90calculation(DFTcalculation):
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

  def __init__(self, directory, **kwargs):
    logger.info("Initializing Wannier90 calculation.")
    super(wannier90calculation, self).__init__(directory)

    self._checkDirectory()    # is this actually a directory
    self._defineCase()        # filename prefix
    self._defineFiles()       # define all possible hr file variants
    self._detectCalculation() # check which hr files we have -> calculation type
    self._checkFiles()        # check if all files are okay
    logger.info("Files successfully loaded.")

    ''' define the usual spin resolved arrays '''
    self.energies        = []
    self.velocities      = []
    self.curvatures      = []
    self.opticalDiag     = []
    self.BopticalDiag    = []
    self.opticalMoments  = []
    self.BopticalMoments = []

    ''' and some flags '''
    self.opticfull  = True
    self.opticdiag  = True
    self.bopticfull = True


  def readData(self):
    ''' get the lattice vectors, number of k-points
        number of projections, projection centers,
        hopping parameters '''
    self._readNnkp()
    self._readHr()
    logger.info("Files successfully read.")

    # internal variables stemming from DFTcalculations
    self.energyBandMax  = self.nproj
    self.opticalBandMin = 0
    self.opticalBandMax = self.nproj

    minarr = np.array([np.min(self.kpoints[:,i]) for i in range(3)])
    maxarr = np.array([np.max(self.kpoints[:,i]) for i in range(3)])
    spacing = (1-maxarr+minarr) # minarr is > 0 for shifted grids
    self.nkx, self.nky, self.nkz = [int(np.around(1./i)) for i in spacing]
    logger.info("   Momentum Grid:  {} x {} x {}".format(self.nkx,self.nky,self.nkz))

    if not (self.nkx*self.nky*self.nkz == self.nkp):
      raise IOError('Irreducible momentum grids not implemented')

    greaterThanOne = ([self.nkx,self.nky,self.nkz] == np.ones(3, dtype=np.int))
    self.dims = np.logical_not(greaterThanOne)
    self.ndim = 3 - np.sum(greaterThanOne)
    logger.info("   Number of dimensions: {}".format(self.ndim))

    self.multiplicity = np.ones((self.nkp,), dtype=int)
    self.weights = self.multiplicity * self.weightsum / (self.nkx*self.nky*self.nkz)
    self.irreducible = False

  def newKmesh(self, kmesh, kshift=False):
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
    if raise_error:
      raise IOError('Provided kmesh does not conform to original kmesh')

    self.nkx, self.nky, self.nkz = kmesh
    self.nkp = self.nkx*self.nky*self.nkz

    self._kmeshx = np.linspace(0,1,self.nkx,endpoint=False)
    self._kmeshy = np.linspace(0,1,self.nky,endpoint=False)
    self._kmeshz = np.linspace(0,1,self.nkz,endpoint=False)

    if kshift:
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
    if kshift: self.kpoints += self._kmeshshift[None,:]

    self.multiplicity = np.ones((self.nkp,), dtype=int)
    self.weights = self.multiplicity * self.weightsum / (self.nkx*self.nky*self.nkz)

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

      for ik in range(self.nkp):
        # do not take the k-iteration into the einstein summation to show progress
        progressBar(ik+1,self.nkp, status='k-points', prefix=prefix)

        ''' FORUIERTRANSFORM hk = sum_r e^{i r.k} * weight(r) * h(r) '''
        rdotk = 2*np.pi*np.einsum('i,ri->r',self.kpoints[ik,:], self.rlist) # no scaling to recip vector / lattice vectors here ?
        ee = np.exp(1j * rdotk) / self.rmultiplicity # rmultiplicity takes care of double counting issues
        hk[ik,:,:] = np.einsum('r,rij->ij', ee, hr)

        ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
        prefactor_r = np.einsum('di,ri->dr', self.rvec, self.rlist)
        hvk[ik,:,:,:] = np.einsum('dr,r,rij->ijd',1j*prefactor_r,ee,hr)

        ''' FORUIERTRANSFORM hvk(j) = sum_r r_j e^{i r.k} * weight(r) * h(r) '''
        prefactor_r2 = np.zeros((6,self.nrp), dtype=np.float64)
        for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
          prefactor_r2[idir,:] = prefactor_r[i,:] * prefactor_r[j,:]
        hck[ik,:,:,:] = np.einsum('dr,r,rij->ijd',-prefactor_r2,ee,hr)

        if peierlscorrection:
          # Jan's code snippet
          # generalized Peierls for multi-atomic unit-cells (and, obviously, supercells)
          distances = self.plist[:,None,:] - self.plist[None,:,:] # nproj nproj 3
          ri_minus_rj = np.einsum('di,abi->abd', self.rvec, distances)
          hvk_correction = 1j * hk[:,:,:,None] * ri_minus_rj[None,:,:,:]


      # eigk  : self.nkp, self.nproj
      # U     : self.nkp, self.nproj, self.nproj
      #       U[0, :, i] is the eigenvector corresponding to ek[0, i]
      # inv(U) @ hk @ U = ek

      ''' this transforms all k points at once '''
      ek, U = np.linalg.eig(hk)
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

    if self.irreducible:
      logger.critical('Detected irreducible grid: optical moments are not symmetrized - will produce wrong results')

    logger.info('Calculating optical elements from derivatives')
    # devine levicity matrix so we can use it in einsteinsummations
    levmatrix = np.zeros((3,3,3), dtype=np.float64)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          levmatrix[i,j,k] = levicivita(i,j,k)

    for ispin in range(self.spins):
      vel = self.velocities[ispin] # nkp, nproj, nproj, 3
      cur = self.curvatures[ispin] # nkp, nproj, nproj, 6

      # transform into matrix form
      curmat  = np.zeros((self.nkp,self.nproj,self.nproj,3,3), dtype=np.complex128)
      curmat[:,:,:, [0,1,2,0,0,1], [0,1,2,1,2,2]] = cur[:,:,:,:]
      curmat[:,:,:, [1,2,2], [0,0,1]] = np.conjugate(curmat[:,:,:, [0,0,1], [1,2,2]])

      # no symmetrization necessary since we _currently_ restrict to reducible grids
      # we exploit the LRTC Fortran differentiation between 3 6 and 9 optical elements
      # this way the array stays float64
      vel2 = np.conjugate(vel[:,:,:,[0,1,2,0,0,1]]) * vel[:,:,:,[0,1,2,1,1,2]]
      if self.ortho:
        vel2 = vel2[:,:,:,:3].real
      else:
        if np.any(np.abs(vel2.imag) > 1e-6):
          temp = vel2.copy()
          vel2 = np.empty((self.nkp,self.nproj,self.nproj,9), dtype=np.float64)
          vel2[:,:,:,:6] = temp
          vel2[:,:,:,6:] = temp[:,:,:,:3].imag
        else:
          vel2 = vel2.real
      self.opticalMoments.append(vel2)
      vel2diag = vel2[:,np.arange(self.nproj),np.arange(self.nproj),:]
      self.opticalDiag.append(vel2diag)

        #           epsilon_cij v_a v_j c_bi -> abc
      mb = np.einsum('cij,knma,knmj,knmbi->knmabc',levmatrix,vel,vel,curmat)
      if np.any(np.abs(mb.imag) > 1e-6):
        # FIXME: do this properly
        logger.info('Detected complex magnetic optical elements. Truncating.')
      mb = mb.real
      self.BopticalMoments.append(mb)
      mbdiag = mb[:,np.arange(self.nproj),np.arange(self.nproj),:,:,:]
      self.BopticalDiag.append(mbdiag)

  def outputData(self, fname, charge, mu=None):
    ''' call the output function '''
    self.charge = charge
    self._calcFermiLevel(mu)
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
    self.fhr         = self.case + '_hr.dat'
    self.fhrup       = self.case + '_hr.datup'
    self.fhrdn       = self.case + '_hr.datdn'
    self.fnnkp       = self.case + '.nnkp'

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
    elif self.calctype == 2:
      self.weightsum = 2
      self.spins = 1
      self.fhraccess = [self.fhr]

    for i in self.fhraccess:
      if not os.path.isfile(i):
        raise IOError('Error: case_hr.dat* is missing')
      if os.stat(i).st_size == 0:
        raise IOerror('Error: case_hr.dat* is empty.')

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
    logger.info('   real_lattice: \n{}'.format(self.rvec))

    ''' FIXME: this does not detect orthogonal bcc and fcc lattices '''
    if np.abs(np.dot(self.rvec[0,:],self.rvec[1,:]) < 1e-6) and \
       np.abs(np.dot(self.rvec[0,:],self.rvec[2,:]) < 1e-6) and \
       np.abs(np.dot(self.rvec[1,:],self.rvec[2,:]) < 1e-6):
      self.ortho = True
    else:
      self.ortho = False
    logger.info('   orthogonal lattice: {}'.format(self.ortho))

    # V = (axb . c)
    self.vol = np.dot(np.cross(self.rvec[0,:], self.rvec[1,:]),self.rvec[2,:])
    logger.info('   deduced volume: {}'.format(self.vol))


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
    logger.info('   recip_lattice: \n{}'.format(self.kvec))

    with open(str(self.fnnkp),'r') as nnkp:
      self.kpoints = []
      while ( not nnkp.readline().startswith('begin kpoints')):
        pass
      self.nkp = int(nnkp.readline())
      logger.info('   number of kpoints: {}'.format(self.nkp))
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
      logger.info('   number of projections: {}'.format(self.nproj))
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
        logger.info('   number of rpoints: {}'.format(self.nrp))
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
