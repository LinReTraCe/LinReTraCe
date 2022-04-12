#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import xml.etree.ElementTree as et
import os
import glob
import sys
import abc
import logging
logger = logging.getLogger(__name__)

if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

import numpy as np
import ase.spacegroup

from structure.aux import progressBar
from structure.es  import ElectronicStructure
from structure     import units

class CustomError(Exception):
  def __init__(self, message):
    super(CustomError, self).__init__(self)
    self.message = message
  def __str__(self):
    return self.message

class DftCalculation(ElectronicStructure, ABC):
  '''
  Abstract Base Class for a generic DFTCalculation.
  Only the methods readData and truncateBands should be called in the main program.
  Everything else should be internal routines not to be accessed by the user.
  These methods are all denoted with a leading underline.
  '''

  def __init__(self, directory):
    super(DftCalculation, self).__init__()
    self.directory  = os.path.abspath(directory) # directory of our calculation
    self.version    = None      # DFT version as string
    self.aseobject  = None      # ase Atoms object
    self.spacegroup = None      # ase Spacegroup object

  @abc.abstractmethod
  def readData(self):
    '''
    Abstract method which is always called when reading in the data
    for all kinds of DFT calculations.
    '''
    pass

  def _checkDirectory(self):
    '''
    Simple method to check if provided directory is actually a directory.
    Also adds the trailing backslash ('/') to the path so we can simply
    add the filenames to the directory to get the full path.
    '''
    if os.path.isfile(self.directory):
      raise IOError("Supplied directory is a file, not a directory: " + self.directory)
    if not os.path.isdir(self.directory):
      raise IOError("Supplied Directory does not exist: " + self.directory)

  def _extractASEinformation(self):
    '''
    Use the aseobject, saved in self.aseobject and extract information
    which is relevant to us.
    '''

    self.vol = self.aseobject.cell.volume # volume of unit cell
    logger.info('Extracting ASE information:')
    logger.info('  unit cell volume [Ang^3] : {}'.format(self.vol))

    ''' This object resembles a 3x3 array whose [i, j]-th element is the jth Cartesian coordinate of the ith unit vector. '''
    self.kvec = self.aseobject.cell.reciprocal()[()]
    self.kvec *= 2*np.pi # since it is not included in the method according to documentation
    self.rvec = self.aseobject.cell[()]
    logger.debug('  real space lattice [Ang]    (rows) :\n{}'.format(self.rvec))
    logger.debug('  reciprocal lattice [Ang^-1] (rows) :\n{}'.format(self.kvec))
    logger.debug('  recip.T @ real / (2pi)=\n{}'.format(self.kvec.T @ self.rvec / 2. / np.pi))

    self.spacegroup = int(ase.spacegroup.get_spacegroup(self.aseobject).no)
    logger.info('  Space group: {}'.format(self.spacegroup))
    self.symop_ase = ase.spacegroup.Spacegroup(self.spacegroup).get_rotations()
    self.nsym_ase  = self.symop_ase.shape[0]
    logger.info('  Symmetry operations: {}'.format(self.nsym_ase))
    logger.debug('  Symmetry matrices:\n{}'.format(self.symop_ase))

    self._detectOrthogonality()
    logger.info('  Orthogonal crystal structure: {}'.format(str(self.ortho)))

  def _detectOrthogonality(self):
    '''
    Create list of symmetry groups which are orthogonal (cubic, tetragonal, orthorhombic).
    cubic:        195 - 230
    tetragonal:   75 - 142
    orthorhombic: 16 - 74
    Check if our extracted spacegroup is in this range of nubmers.
    '''

    orthogroups = list(range(16,143))
    orthogroups += list(range(195,231))
    self.ortho = (self.spacegroup in orthogroups)
