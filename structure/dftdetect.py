#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import glob
import logging
logger = logging.getLogger(__name__)

from structure.wien2k import Wien2kCalculation
from structure.vasp   import VaspCalculation

class DftDetection(object):
  def __init__(self, directory):
    self.directory      = os.path.abspath(directory) # directory of our calculation
    self._checkDirectory()
    self.methods = [self._checkW2Kcalculation, self._checkVASPcalculation]

  def __repr__(self):
    return ('DFTdetection(directory={0.directory})'.format(self))

  def __len__(self):
    return len(self.methods)

  def __call__(self):
    '''
    Functor that loops all the possible available DFT calculations.
    Return a dictionary if it finds something.
    None otherwise
    '''
    for imethod in self.methods:
      returndic = imethod()
      if returndic is not None:
        return returndic
    return None

  def _checkDirectory(self):
    ''' Simple method to check if provided directory is actually a directory '''
    if os.path.isfile(self.directory):
      raise IOError("Supplied directory is a file, not a directory: " + self.directory)
    if not os.path.isdir(self.directory):
      raise IOError("Supplied Directory does not exist: " + self.directory)

  def _checkW2Kcalculation(self):
    '''
    Detect the Wien2K calculation by the existance of the
    scf file -> also extract the case prename
    '''

    scf = os.path.join(self.directory,'*.scf')
    files = glob.glob(scf)
    if len(files) >= 1:
      if len(files) > 1:
        logger.warn('Detected more than 1 scf file in provided folder: Choosing {}'.format(files[0]))

      purescf = os.path.basename(files[0])
      temp = purescf.split('.')  # abc.def.scf -> [abc,def,scf]
      case = '.'.join(temp[:-1]) # abc.def
      return {'dft': Wien2kCalculation, 'case': case}
    else:
      return None

  def _checkVASPcalculation(self):
    '''
    Detect the vasp calculation by the existance of the vasprun.xml file
    '''

    vasprun = os.path.join(self.directory,'vasprun.xml')
    files = glob.glob(vasprun)
    if len(files) == 1: # there is only one vasprun.xml
      return {'dft': VaspCalculation}
    else:
      return None
