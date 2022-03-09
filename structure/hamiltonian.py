#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import logging
import collections
logger = logging.getLogger(__name__)

import numpy as np

from structure.es import ElectronicStructure
from structure.model import Model


class hamiltonian_wannier90(ElectronicStructure):
  '''
  Class to describe Wannier like Hamiltonians
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
