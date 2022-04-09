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

  def _getAuxiliaryInformation(self):
    '''
    Use the aseobject and extract information
    which is relevant to us.
    '''

    self.vol = self.aseobject.get_volume() # volume of unit cell
    self.atms = self.aseobject.get_global_number_of_atoms() # total number of atoms
    # laa = self.aseobject.get_cell_lengths_and_angles()
    laa = self.aseobject.cell.cellpar()
    self.latLength = np.array(laa[:3])
    self.latAngle = np.array(laa[3:])
    self.spacegroup = ase.spacegroup.get_spacegroup(self.aseobject)
    # logger.info('ASE: Detected space group {}'.format(self.spacegroup))

    self._detectOrthogonality()
    if self.ortho:
      logger.info('ASE: Detected orthogonal crystal structure')
    else:
      logger.info('ASE: Detected non-orthogonal crystal structure')

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
    if self.spacegroup.no in orthogroups:
      self.ortho = True
    else:
      self.ortho = False
