#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import xml.etree.ElementTree as et
import os
import logging
logger = logging.getLogger(__name__)

import numpy as np
import ase.io


from structure.aux import progressBar
from structure.dft import DftCalculation
from structure     import units

class VaspCalculation(DftCalculation):
  '''
  VASP calculation derived from our abstract base class.
  Since Vasp does not calculate transition dipole moments (optical elements)
  we can default the optical flag to False.
  All the data is read-in from the vasprun.xml file with the help of ElementTree.
  '''

  def __init__(self, directory, **kwargs):
    logger.info("\nInitializing VASP calculation.")
    super(VaspCalculation, self).__init__(directory)

    self.version = None
    self._checkDirectory() # from base class
    self._defineFiles()
    self._checkFiles()
    logger.info("Files sucessfully loaded.")
    self._extractASEinformation()

  def __repr__(self):
    return ('vaspcalculation(directory={0.directory!r}'.format(self))

  # we only care about the vasprun file
  def _defineFiles(self):
    self.fvasprun    = os.path.join(self.directory,'vasprun.xml')
    self.foutcar     = os.path.join(self.directory,'OUTCAR')
    self.aseobject   = ase.io.read(self.fvasprun)

  def _checkFiles(self):
    if not os.path.isfile(self.fvasprun):
      raise IOError("Error: vasprun.xml file missing.")
    if os.stat(str(self.fvasprun)).st_size == 0:
      raise IOError("Error: vasprun.xml is empty.")

  def readData(self):
    logger.info("Reading: {}".format(self.fvasprun))
    self.xmlroot = et.parse(self.fvasprun).getroot()
    self._read_version()
    self._read_ispin()
    self._read_lnoncollinear()
    self._read_nelect()
    self._read_vol()
    self._read_kpointlist()
    self._read_kdivisors()
    self._read_weights()
    self._read_energies()
    logger.info("Files successfully read.")
    self._calcFermiLevel()
    self._checkFermiLevel() # TODO: decide if I want to remove this =?

  def _checkFermiLevel(self):
    try:
      efer = float(self.xmlroot.find('calculation/dos/i[@name="efermi"]').text)
      logger.info("  Cross-check: Fermi level from vasp DOS: {}".format(efer))
    except Exception:
      pass # this value is not important

  def _read_version(self):
    try:
      self.version = self.xmlroot.find('generator/i[@name="version"]').text
    except Exception as s:
      pass
    logger.info("  VASP version: {}".format(self.version))

  def _read_ispin(self):
    # because VASP uses an illegal output by using non-xml-conform tags with spaces in them ...
    # please fucking fix that
    self.spins = np.nan
    for paras in self.xmlroot.findall('parameters/separator[@name="electronic"]'):
      for value in paras.iter():
        # print(value.attrib, value.text, value.tag)
        if value.attrib['name'] == 'ISPIN':
          self.spins = int(value.text)
    if self.spins == np.nan:
      raise IOError('Could not find ISPIN variable in vasprun.xml')

    if (self.spins == 1):
      self.weightsum = 2
    else:
      self.weightsum = 1
    logger.info("  Number of inequivalent spins: {}".format(self.spins))

  def _read_lnoncollinear(self):
    # same xml access problems as with the spins
    self.lnc = np.nan
    for paras in self.xmlroot.findall('parameters/separator[@name="electronic"]'):
      for value in paras.iter():
        if value.attrib['name'] == 'LNONCOLLINEAR':
          self.lnc = value.text.strip()
    if self.lnc == np.nan:
      raise IOError('Could not find ISPIN variable in vasprun.xml')

    if self.lnc == 'T' and self.spins == 1:
      self.weightsum = 1 # overwrite it
    logger.info("  Non-collinear calculation: {}".format(self.lnc))

  def _read_nelect(self):
    try:
      self.charge = float(self.xmlroot.find('parameters/separator[@name="electronic"]/i[@name="NELECT"]').text)
    except Exception:
      raise IOError('Could not find NELECT variable in vasprun.xml')
    logger.info("  Number of electrons: {}".format(self.charge))

  def _read_vol(self):
    try:
      self.vol = float(self.xmlroot.find('structure[@name="finalpos"]/crystal/i[@name="volume"]').text)
    except Exception:
      raise IOError('Could not find volume variable in vasprun.xml')
    logger.info('  Unit cell volume: {} [Angstrom^3]'.format(self.vol))

  def _read_kpointlist(self):
    try:
      arr = self.xmlroot.find('kpoints/varray[@name="kpointlist"]')
      kpointlist = []
      for ikp in arr:
        kpointlist.append(ikp.text.split())

      self.nkp = len(kpointlist)
      self.kpoints = np.array(kpointlist, dtype=np.float64)
    except Exception:
      raise IOError('Error occured during read-in of k-points in {}'.format(str(self.fvasprun)))

    logger.info("  Number of k-points: {}".format(self.nkp))

  def _read_kdivisors(self):
    try:
      arr = self.xmlroot.find('kpoints/generation/v[@name="divisions"]')
      divisor = arr.text.split()
      divisor = np.array(divisor, dtype=np.int)
      self.nkx, self.nky, self.nkz = divisor
      logger.info("  Momentum grid: {} {} {}".format(*divisor))
      self.dims = np.logical_not(divisor == np.ones(3, dtype=np.int))
      self.ndim = 3 - np.sum(divisor == np.ones(3, dtype=np.int))

      self.irreducible = not (self.nkx*self.nky*self.nkx == self.nkp)
      if not self.irreducible:
        self.nsym     = 1
        self.symop    = np.zeros((1,3,3), dtype=np.float64)
        self.invsymop = np.zeros((1,3,3), dtype=np.float64)
        self.symop[0,:,:] = np.diag((1,1,1))
        self.invsymop[0,:,:] = np.diag((1,1,1))
    except Exception:
      raise IOError('Error occured during read-in of k-point divisors {}'.format(str(self.fvasprun)))

    logger.info("  Irreducible grid: {}".format(str(self.irreducible)))
    logger.info("  Number of dimensions: {}".format(self.ndim))

  def _read_weights(self):
    try:
      arr = self.xmlroot.find('kpoints/varray[@name="weights"]')
      weightlist = []
      for ikp in arr:
        weightlist.append(ikp.text)

      # weight_i = mult_i / sum_i mult_i
      # sum_i mult_i = nkx * nky * nkz

      self.weights = np.array(weightlist, dtype=np.float64) # always normalized to 1
      self.multiplicity = np.around(self.weights * self.nkx * self.nky * self.nkz).astype(int)

      self.weights *= self.weightsum # adjust for spin unpolarized calculations


    except Exception:
      raise IOError('Error occured during read-in of weights in {}'.format(str(self.fvasprun)))

  def _read_energies(self):
    try:
      for ispin in range(self.spins):
        energylist = []
        spinarr = self.xmlroot.find('calculation/eigenvalues/array/set/set[@comment="spin {}"]'.format(ispin+1))
        for ikp in range(self.nkp):
          energylist.append([])
          kpointarr = spinarr.find('set[@comment="kpoint {}"]'.format(ikp+1))
          for ene in kpointarr:
            energylist[ikp].append(ene.text.split()[0])
        enearray = np.array(energylist, dtype=np.float64)
        self.energyBandMax = enearray.shape[-1]
        self.energies.append(enearray)
    except Exception:
      raise IOError('Error occured during read-in of energies in {}'.format(str(self.fvasprun)))
    logger.info("  Number of bands: {}".format(self.energyBandMax))

