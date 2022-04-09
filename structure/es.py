#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import abc
import logging
logger = logging.getLogger(__name__)

# python 2 & 3 compatible usage of abstract base classes
if sys.version_info >= (3, 4):
  ABC = abc.ABC
else:
  ABC = abc.ABCMeta('ABC', (), {})

import scipy.optimize
import numpy as np

class Converged(Exception):
  def __init__(self, mu):
    super(Converged, self).__init__(self)
    self.mu = mu

class ElectronicStructure(ABC):
  '''
  Abstract Parent class for all electronic structures.
  Here we define all the common elements:
  Number of spins, k-points, multiplicity, weights, etc.

  This class provides the internal methods to find the 'DFT' mu, i.e.
  the chemical potential for T=0: _calcFermiLevel(mu=None)
  Before using the output function h5out in inout.py this must be called
  in order to set the required internal variables that are output
  '''

  def __init__(self):
    # information about kmesh
    self.nkp            = 0    # number of k-points
    self.nkx            = 0    # number of k-points in x-direction
    self.nky            = 0    # number of k-points in x-direction
    self.nkz            = 0    # number of k-points in x-direction

    self.multiplicity   = None # multiplicity of the k-points
    self.weights        = None # weights of the k-points
    self.weightsum      = 0    # sum of the weights
    self.kpoints        = None # list of the k-points ... shape [nkp,3] float64
    self.irreducible    = False# irreducible or reducible k-grid
    self.ortho          = False

    self.spins          = 1    # number of spins we are considering
    self.vol            = 0    # volume of the unit cell in AA^3
    self.charge         = 0    # charge in the given bands
    self.mu             = 0    # chemical potential

    # symmetries
    self.nsym           = 0    # number of symmetry operations
    self.symop          = None # symmetry operations
    self.invsymop       = None # inverse symmetre operations

    self.ndim           = -1   # number of dimensions
    self.dims           = np.array([False,False,False]) # valid dimension, i.e. k_i > 1

    # information about bands
    self.energies       = []   # list entries for energies elements are the respective spins
                               # entries themselves are shape [nkp,bands] float64
    self.energyBandMax  = 0    # band maximum for energies

    # information about optical elements
    self.opticdiag       = False # intra band optical elements
    self.opticalDiag     = []    # list for band-diagonal optical elements
    self.opticfull       = False # full optical elements (inter+intra)
    self.opticalMoments  = []    # list for full optical elements
    # same thing for the b-field quantities
    self.bopticdiag      = False
    self.bopticalDiag    = []
    self.bopticfull      = False
    self.bopticalMoments = []   # list for full optical elements
    self.opticalBandMin  = 0    # band interval minimum for optical elements
    self.opticalBandMax  = 0    # band interval maximum for optical elements

  def _defineDimensions(self):
    '''
    Count dimension as every k-axis with more than one k-point
    1-dimension k-axis is enforced to be kx
    2-dimension k-axis are enforced to be kx and ky
    '''
    self.ndim = 0
    self.dims = []
    for i in [self.nkx, self.nky, self.nkz]:
      if i < 1:
        raise ValueError("Number of kpoints in each direction have to be positive")
      if i > 1:
        self.ndim += 1
        self.dims.append(True)
      else:
        self.dims.append(False)
    self.dims = np.array(self.dims)
    logger.info('Detected {} dimensions.'.format(self.ndim))

  def _calcOccupation(self, mu):
    '''
    Calculate the deviation of the occupation the the given charge in the system.
    Here we use the Fermi function at T=0 (theta function).
    i.e. energies smaller than the chemical potential are fully occupied
    while energies larger than the chemical potential are copletely empty
    '''

    energies = np.zeros((self.spins, self.nkp, self.energyBandMax), dtype=np.float64)
    for ispin in range(self.spins):
      energies[ispin] = self.energies[ispin]

    # T = 0 fermi function -> theta function
    mask = (energies-mu) < 0.0
    energies[mask] = 1.0
    energies[np.logical_not(mask)] = 0.0
    dev = np.sum(energies * self.weights[None,:,None]) - self.charge

    # this is a nasty work-around for the whole gap-thingy
    # to make sure we are where we want to be
    if abs(dev) < np.min(self.weights)/2.:
      raise Converged(mu)

    return dev

  def _calcFermiLevel(self, mu=None):
    '''
    Use the energies and the charge: calculate the T=0 chemical
    potential and determine if a gap exists.
    If one exists, determine the gap size and between which bands
    it lies
    '''

    # get the bisection start and end points
    x0 = np.min(self.energies[0][:,0])
    x1 = np.max(self.energies[0][:,0])

    for ispin in range(self.spins):
      for iband in range(self.energyBandMax):
        ymin = np.min(self.energies[0][:,iband])
        ymax = np.max(self.energies[0][:,iband])
        if ymin < x0:
          x0 = ymin
        if ymax > x1:
          x1 = ymax

    # safety offset
    x0 -= 0.1
    x1 += 0.1

    if mu is None:
      try:
        mu_sol = scipy.optimize.bisect(self._calcOccupation, x0, x1)
      except Converged as c: # work-around so we always get the band-gap correct
        mu_sol = c.mu
      except ValueError:
        raise Exception('Fermi Level Calculation: Bisection failed.')
    else:
      mu_sol = mu

    self.mu = mu_sol

    self.gapped     = []
    self.gap        = []
    self.ecb        = []
    self.evb        = []
    self.cb         = []
    self.vb         = []

    # detect the spin-dependent gap
    for ispin in range(self.spins):
      locgapped = True
      for iband in range(self.energyBandMax):
        ene = self.energies[ispin][:,iband]
        enemin = np.min(ene)
        enemax = np.max(ene)
        if mu_sol > enemin and mu_sol < enemax: # it cuts through this band
          locgapped = False
          self.cb.append(iband)
          self.vb.append(iband)
          break
        if enemin > mu_sol: # band is for the first time above chemical potential
          self.cb.append(iband)
          self.vb.append(iband-1)
          break
        if mu_sol == enemin:
          self.cb.append(iband)
          self.vb.append(iband-1)
          break
        if mu_sol == enemax:
          self.vb.append(iband)
          self.cb.append(iband+1)
          break
      else:
        locgapped = False

      # save the gap data
      if locgapped:
        enevalence    = self.energies[ispin][:,self.vb[ispin]]
        eneconduction = self.energies[ispin][:,self.cb[ispin]]
        gap = np.min(eneconduction) - np.max(enevalence)
        if gap < 1e-13: # required for 'touching' bands like in graphene
          self.gapped.append(False)
          self.gap.append(np.nan)
          self.ecb.append(np.nan)
          self.evb.append(np.nan)
        else:
          self.gapped.append(True)
          self.gap.append(gap)
          self.ecb.append(np.min(eneconduction))
          self.evb.append(np.max(enevalence))
      else:
        self.gapped.append(False)
        self.gap.append(np.nan)
        self.ecb.append(np.nan)
        self.evb.append(np.nan)

    # adjust mu for 1 spin and gapped to be exactly in the middle
    # this is correct for T = 0
    if self.spins==1 and self.gapped[0]:
      self.mu = (self.ecb[0] + self.evb[0]) / 2.
      logger.info("Putting Chemical potential in the middle of the gap.")

    # adjust mu for 2 spins
    # this should be correct for T = 0
    if self.spins==2 and self.gapped[0] and self.gapped[1]:
      self.mu = (min(self.ecb) + max(self.evb)) /2.
      logger.info("Putting Chemical potential in the middle of the common gap.")

    # notify user
    logger.info("Chemical potential: {} [eV]".format(self.mu))
    for ispin in range(self.spins):
      if self.gapped[ispin]:
        logger.info('  Spin: {} / {}: Found energy gap: {} [eV]'.format(ispin+1,self.spins,self.gap[ispin]))
        logger.info('         vbmax: {} [eV] - cbmin: {} [eV]'.format(self.evb[ispin],self.ecb[ispin]))
      else:
        logger.info('  Spin: {} / {}: no energy gap'.format(ispin+1,self.spins))

  @staticmethod
  def distributeWorkLoad(datainterval, processes):
    '''
    Distrube a datainterval to the given number of processes.
    Return the working ranges in form of a list of tuples.
    Might be useful for more tasks.
    '''

    displ=[] # displacement
    displ.append(0)
    rct=[]   # receive count

    # calculate the count and displaceent by continuisly
    # dividing the remaining datarange by the remaining processes
    # until we reach the end
    for i in range(processes-1):
      rct.append((datainterval-displ[i])//(processes-i))
      displ.append(rct[i]+displ[i])
    rct.append(datainterval- displ[processes-1])

    # save the working ranges in form of a list of tuples
    # where the tuple contains the start and (not-included) stop k-point
    myranges = []
    for i in range(processes):
      myranges.append((displ[i],displ[i]+rct[i]))

    return myranges
