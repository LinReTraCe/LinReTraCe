#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import logging
logger = logging.getLogger(__name__)

import numpy as np
import h5py

class LRTCscat(object):
  '''
    class to define the most general scattering dependencies
    is instantiates with a LRTC energy input that sets most variables
  '''

  def __init__(self, fname):
    self.efile = fname
    self._loadparams()
    self.inputmethod = input if sys.version_info >= (3, 0) else raw_input

    self.scattering = None
    self.qpweight   = None
    self.bandshift  = None

    self.mumode   = False
    self.tempmode = False

  def _loadparams(self):
    ''' check energy file and load in variables '''
    with h5py.File(self.efile,'r') as h5in:
      if h5in.attrs['identifier'] == 'LRTCinput':
        logger.info('Detected LRTC energy file: {}'.format(self.efile))
      self.spins  = h5in['.bands/ispin'][()]
      self.nbands = h5in['.bands/energyBandMax'][()]
      self.nkp    = h5in['.kmesh/nkp'][()]
      self.mudft  = h5in['.bands/mu'][()] # chemical potential from electronic structure
      self.kgrid  = h5in['.kmesh/points'][()]

      self.energies = np.zeros((self.spins,self.nkp,self.nbands), dtype=np.float64)
      if self.spins == 1:
        self.energies[0,...] = h5in['energies'][()]
      else:
        self.energies[0,...] = h5in['up/energies'][()]
        self.energies[1,...] = h5in['dn/energies'][()]

  def getDependencies(self):
    return self.spins, self.nkp, self.nbands

  def getMomentumGrid(self):
    return self.kgrid

  def getEnergies(self):
    return self.energies

  def defineTemperatures(self, tmin, tmax, nt, tlog=False):
    ''' Define temperature range given by the provided temperatures
        also set the interla mumode switch
    '''
    self.mumode   = False
    self.tempmode = True
    kB            = 8.61733034e-5 # eV / K
    self.nt       = nt
    self.tmin     = tmin
    self.tmax     = tmax

    if not tlog:
      self.temps = np.linspace(tmin,tmax,nt, endpoint=True) # T interval on K
    else:
      self.temps = np.logspace(np.log10(tmin),np.log10(tmax),nt, endpoint=True) # T interval on K
    self.betas = 1.0/(self.temps * kB) # beta interval on eV^-1
    self.mus = np.ones_like(self.temps, dtype=np.float64) * self.mudft

    return self.temps

  def defineChemicalPotentials(self, temp, mumin, mumax, nmu, mlog=False, muabs=False):
    '''
        Define temperature which to run on and
        define chemical potential range given by the provided parameters
        also set the internal mumode switch
    '''
    self.tempmode = False
    self.mumode   = True
    kB            = 8.61733034e-5 # eV / K

    self.nmu      = nmu
    self.mumin    = mumin
    self.mumax    = mumax

    self.nt       = nmu
    self.tmin     = temp
    self.tmax     = temp

    if not muabs:
      self.mumin += self.mudft
      self.mumax += self.mudft

    if not mlog:
      self.mus   = np.linspace(self.mumin, self.mumax, self.nmu, endpoint=True) # eV
    else:
      self.mus   = np.logspace(np.log10(mumin),np.log10(mumax),nmu, endpoint=True)

    self.temps = np.ones_like(self.mus, dtype=np.float64) * temp
    self.betas = 1.0/(self.temps * kB) # beta interval on eV^-1

    return self.mus

  def defineScatteringRates(self, scattering, qpweight=None, bandshift=None):
    '''
        readin scattering rate array
        and optionally quasiparticle weight array
        and bandshift array
        check array sizes via assertions to shape and dtype

        input arrays have the generic shape [nt, spins, nkp, nbands]
        where nt ... number of temp or mu steps
              spins ... 1 or 2
              nkp   ... number of k-points
              nbands ... number of bands

        however we allow nkp==1 and nbands==1
        this is interpreted as no k and no band dependence, respectively
    '''

    self.scattering = scattering
    self.qpweight   = qpweight
    self.bandshift  = bandshift

    try:
      assert(self.scattering.shape[0] == self.nt)
      assert(self.scattering.shape[1] == self.spins)
      assert(self.scattering.shape[2] == 1 or self.scattering.shape[2] == self.nkp)
      assert(self.scattering.shape[3] == 1 or self.scattering.shape[3] == self.nbands)
      assert(self.scattering.dtype == np.float64)
    except:
      raise IOError('Scattering array shape mismatch')

    if np.any(self.scattering < 1e-14):
      logger.warning('\n\nScattering rates below 10^{-14} eV detected.\nPlease introduce cutoff to avod numberical problems.\n\n')

    try:
      if qpweight is not None:
        assert(self.qpweight.shape == self.scattering.shape)
        assert(self.qpweight.dtype == np.float64)
      else:
        self.qpweight = np.ones_like(self.scattering, dtype=np.float64)
    except:
      raise IOError('Quasi-particle weight array shape mismatch')

    try:
      if bandshift is not None:
        assert(self.bandshift.shape == self.scattering.shape)
        assert(self.bandshift.dtype == np.float64)
    except:
      raise IOError('Bandshift array shape mismatch')

  def createOutput(self, outname):
    if (not self.mumode) and (not self.tempmode):
      raise IOError('Neither temperature or mu mode have been set')

    if os.path.isfile(outname):
      print('Output file already exists!')
      overwrite = self.inputmethod('Overwrite? (y/N) ').strip().lower()
      if overwrite not in ['y','yes']:
        logger.info('File will not be overwritten.\nExiting.')
        return

    with h5py.File(outname,'w') as h5out:
      h5out['.quantities/iSpin']    = self.spins
      h5out['.quantities/nbands']   = self.nbands
      h5out['.quantities/nkp']      = self.nkp

      h5out['.quantities/tempmode'] = self.tempmode
      h5out['.quantities/mumode']   = self.mumode

      h5out['.quantities/tempAxis'] = self.temps
      h5out['.quantities/betaAxis'] = self.betas
      h5out['.quantities/muAxis']   = self.mus
      h5out['.quantities/nT']       = self.nt
      h5out['.quantities/Tmin']     = self.tmin
      h5out['.quantities/Tmax']     = self.tmax

      # common output routines

      for ispin in range(self.spins):
        if self.spins == 1:
          prefix = '/'
        else:
          if ispin == 0:
            prefix = '/up/'
          else:
            prefix = '/dn/'

        independence = self._checkIndependence(ispin)

        if independence:
          h5out['{}step/{:06}/scatrate'.format(prefix,1)] = self.scattering[0,ispin,:,:]
          h5out['{}step/{:06}/qpweight'.format(prefix,1)] = self.qpweight[0,ispin,:,:]
          if self.bandshift is not None:
            h5out['{}step/{:06}/bandshift'.format(prefix,1)] = self.bandshift[0,ispin,:,:]

          ''' hard links '''
          for iT in range(1,self.nt):
            h5out['{}step/{:06}'.format(prefix,iT+1)] = h5out['{}step/{:06}'.format(prefix,1)]
        else:
          for iT in range(self.nt):
            h5out['{}step/{:06}/scatrate'.format(prefix,iT+1)] = self.scattering[iT,ispin,:,:]
            h5out['{}step/{:06}/qpweight'.format(prefix,iT+1)] = self.qpweight[iT,ispin,:,:]
            if self.bandshift is not None:
              h5out['{}step/{:06}/bandshift'.format(prefix,iT+1)] = self.bandshift[iT,ispin,:,:]

    logger.info('Successfully created scattering file: {}'.format(outname))

  def _checkIndependence(self, ispin):
    '''
    check whether we got a temperature / mu independent input for the provided spin
    if all temperature / mu steps are identical to the first one return True
    otherwise return False
    '''

    steps = self.nt
    if steps == 1: return True

    compare_scattering = self.scattering[0,ispin,...]
    compare_qpweight   = self.qpweight[0,ispin,...]
    if self.bandshift is not None:
      compare_bandshift = self.bandshift[0,ispin,...]

    independence = True
    for i in range(1,self.nt):
      if not np.allclose(compare_scattering, self.scattering[i,ispin,...]):
        independence = False
        break
      if not np.allclose(compare_qpweight, self.qpweight[i,ispin,...]):
        independence = False
        break
      if self.bandshift is not None:
        if not np.allclose(compare_bandshift, self.bandshift[i,ispin,...]):
          independence = False
          break

    return independence
