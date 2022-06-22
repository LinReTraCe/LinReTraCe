#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import xml.etree.ElementTree as et
import os
import glob
import logging
logger = logging.getLogger(__name__)

from structure.wien2k import Wien2kCalculation
from structure.vasp   import VaspCalculation

class DftDetection(object):
  def __init__(self, path):
    self.path      = os.path.abspath(path) # directory of our calculation
    self._checkPath()
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

  def _checkPath(self):
    ''' Simple method to check if provided path exists '''
    if not os.path.exists(self.path):
      raise IOError("Supplied path does not exist: " + self.path)

  def _checkW2Kcalculation(self):
    '''
    Detect the Wien2K calculation by the existance of the
    struct file -> also extract the case prename
    '''

    if not os.path.isdir(self.path):
      return None

    scf = os.path.join(self.path,'*.scf')
    files = glob.glob(scf)
    if len(files) >= 1:
      if len(files) > 1:
        files.sort() # this is necessary for minimization jobs where there are multiple scfs: a.scf - a_1.scf - a_2.scf - etc.
        logger.warn('Detected more than 1 scf file in provided folder: Choosing {}'.format(files[0]))

      purescf = os.path.basename(files[0])
      temp = purescf.split('.')  # abc.def.scf -> [abc,def,scf]
      case = '.'.join(temp[:-1]) # abc.def
      return {'dft': Wien2kCalculation, 'case': case}

    in2 = os.path.join(self.path,'*.in2*')
    files = glob.glob(in2)
    if len(files) >= 1:
      if len(files) > 1:
        logger.warn('Detected more than 1 scf file in provided folder: Choosing {}'.format(files[0]))

      purein2 = os.path.basename(files[0])
      temp = purein2.split('.')  # abc.def.scf -> [abc,def,scf]
      case = '.'.join(temp[:-1]) # abc.def
      return {'dft': Wien2kCalculation, 'case': case}

    return None

  def _checkVASPcalculation(self):
    '''
    Detect the vasp calculation by the existance of the vasprun.xml file
    '''

    if os.path.isdir(self.path):
      directory = os.path.abspath(self.path)
      vasprun    = os.path.join(directory,'vasprun.xml')
    elif os.path.isfile(self.path):
      vasprun    = os.path.abspath(self.path)

    try:
      xmlroot = et.parse(vasprun).getroot()
      version = xmlroot.find('generator/i[@name="version"]').text
    except:
      return None

    return {'dft': VaspCalculation}
