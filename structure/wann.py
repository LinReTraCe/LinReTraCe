#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import glob
import sys
import logging
logger = logging.getLogger(__name__)

from structure.aux import progressBar
from structure.dft  import DFTcalculation
from structure.es  import ElectronicStructure

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

    self._checkDirectory()
    self._defineCase()
    self._defineFiles()
    # self._detectCalculation() # restrict spin 1 for the time being
    self.spins = 1
    self._checkFiles()
    logger.info("Files successfully loaded.")

  def readData(self):
    ''' get the lattice vectors, number of k-points
        number of projections, projection centers,
        hopping parameters '''
    self._readNnkp()
    self._readHr()
    logger.info("Files successfully read.")

  def transformData(self):
    ''' calculate e(r), v(r), c(r) ? '''
    pass
    self._calcFermiLevel()

  def _defineFiles(self):
    self.fhr         = self.case + '_hr.dat'
    self.fnnkp       = self.case + '.nnkp'

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

  def _checkFiles(self):
    if not os.path.isfile(self.fhr):
      raise IOError('Error: case_hr.dat is missing')
    if os.stat(str(self.fhr)).st_size == 0:
      raise IOerror('Error: case_hr.dat is empty.')
    if not os.path.isfile(self.fnnkp):
      raise IOError('Error: case_hr.dat is missing')
    if os.stat(str(self.fnnkp)).st_size == 0:
      raise IOerror('Error: case_hr.dat is empty.')

  def _readNnkp(self):
    self.rvec = []
    logger.info("Reading: {}".format(self.fnnkp))
    with open(str(self.fnnkp),'r') as nnkp:
      read_rvec = False
      for line in nnkp:
        if line.startswith('begin real_lattice'):
          read_rvec = True
          continue
        if read_rvec:
          rx, ry, rz = float(line[:12]), float(line[12:24]), float(line[24:36])
          self.rvec.append([rx,ry,rz])
        if len(self.rvec) == 3:
          break
    self.rvec = np.array(self.rvec, dtype=np.float64)
    logger.info('   real_lattice: {}'.format(self.rvec))

    self.kvec = []
    with open(str(self.fnnkp),'r') as nnkp:
      read_kvec = False
      for line in nnkp:
        if line.startswith('begin recip_lattice'):
          read_kvec= True
          continue
        if read_kvec:
          kx, ky, kz = float(line[:12]), float(line[12:24]), float(line[24:36])
          self.kvec.append([kx,ky,kz])
        if len(self.kvec) == 3:
          break
    self.kvec = np.array(self.kvec, dtype=np.float64)
    logger.info('   recip_lattice: {}'.format(self.kvec))

    with open(str(self.fnnkp),'r') as nnkp:
      read_kp = False
      self.nkp = None
      self.klist = []
      for line in nnkp:
        if line.startswith('begin kpoints'):
          read_kp = True
          continue
        if read_kp:
          if self.nkp is None:
            self.nkp = int(line)
            logger.info('{}'.format(self.nkp))
          else:
            kx, ky, kz = float(line[:14]), float(line[14:28]), float(line[28:42])
            self.klist.append([kx,ky,kz])
          if len(self.klist) == self.nkp:
            break
    self.klist = np.array(self.klist, dtype=np.float64)
    logger.info('   number of kpoints: {}'.format(self.nkp))

    with open(str(self.fnnkp),'r') as nnkp:
      read_proj = False
      self.nproj = None
      self.plist = []
      skip = False
      for line in nnkp:
        if line.startswith('begin projections'):
          read_proj = True
          continue
        if skip:
          skip = False
          continue
        if read_proj:
          if self.nproj is None:
            self.nproj = int(line)
          else:
            rx, ry, rz = float(line[:10]), float(line[10:21]), float(line[21:32])
            skip = True # only every other line
            self.plist.append([rx,ry,rz])
          if len(self.plist) == self.nproj:
            break
    self.plist = np.array(self.plist, dtype=np.float64)
    logger.info('   number of projections: {}'.format(self.nproj))

  def _readHr(self):
    pass


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
