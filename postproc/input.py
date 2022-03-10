#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import argparse
import warnings
from collections import OrderedDict

import numpy as np
with warnings.catch_warnings():
  warnings.filterwarnings("ignore",category=FutureWarning)
  import h5py

from structure import es

class LRTCinput(object):
  '''
  Input class for the energy input files of LRTC
  '''

  def __init__(self, fname):
    self.fname    = fname.strip()
    self._parse_file()

    if sys.version_info >= (3, 0): #  for numpy
      self.textpipe = sys.stdout.buffer
    else:
      self.textpipe = sys.stdout

  def __repr__(self):
    return ('LRTCinput(fname={0.fname!r})'.format(self))

  def _parse_file(self):
    '''
    Parse provided file to check the energies group
    if it exists we implicitly detected an input file
    If nothing is found raise an IOError
    '''
    try:
      with h5py.File(self.fname,'r') as h5:
        try:
          ene =  h5['energies'] # check if its there
        except:
          eneup =  h5['up/energies'] # check if its there
    except:
      raise IOError('Provided file is not an LRTC input file.')

  def outputList(self):
    '''
    List all output options
    '''
    barlength = 55
    print('\n{:<12}  {}'.format('Key', 'Description'))
    print(barlength*u'\u2500')
    print('{:<12}  {}'.format('info', 'band structure, optical information'))
    print('{:<12}  {}'.format('dos',  'density of states / number of states'))
    print(barlength*u'\u2500')

  def outputStructure(self):
    '''
    List input information ... just hardcoded stuff
    '''
    with h5py.File(self.fname,'r') as h5:
      print('\nINPUT FILE')
      print('  name: {}'.format(self.fname))

      print('\nKMESH')
      print('  k-points: {}'.format(h5['.kmesh/nkp'][()]))
      weights = h5['.kmesh/weights'][()]
      print('  irreducible: {}'.format(h5['.kmesh/irreducible'][()]))

      print('\nBANDS')
      spins = h5['.bands/ispin'][()]
      print('  spins: {}'.format(spins))
      print('  bands: {}'.format(h5['.bands/energyBandMax'][()]))
      print('  charge: {}'.format(h5['.bands/charge'][()]))
      if spins==1:
        gapped = h5['.bands/bandgap/gapped'][()]
        print('  gapped: {}'.format('yes' if gapped else 'no'))
        if gapped:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/cband'][()]))
      else:
        gappedup = h5['.bands/bandgap/up/gapped'][()]
        print('  spin up gapped: {}'.format('yes' if gappedup else 'no'))
        if gappedup:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/up/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/up/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/up/cband'][()]))
        gappeddn = h5['.bands/bandgap/dn/gapped'][()]
        print('  spin dn gapped: {}'.format('yes' if gappeddn else 'no'))
        if gappeddn:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/dn/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/dn/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/dn/cband'][()]))
      print('  chemical potential: {0:.4f} eV'.format(h5['.bands/mu'][()]))

      print('\nUNIT-CELL')
      print('  volume: {0:.4f} angstrom**3'.format(h5['.unitcell/volume'][()]))

      print('\nOPTICAL DATA')
      if spins==1:
        print('  intra elements   : {}'.format('yes' if 'momentsDiagonal' in h5 else 'no'))
        print('  inter elements   : {}'.format('yes' if 'kPoint/{:010}/moments'.format(1) in h5 else 'no'))
        print('  band derivatives : {}'.format('yes' if 'derivatives' in h5 else 'no'))
        print('  band curvatures  : {}'.format('yes' if 'curvatures' in h5 else 'no'))
      else:
        print('  intra elements   : {}'.format('yes' if 'up/momentsDiagonal' in h5 else 'no'))
        print('  inter elements   : {}'.format('yes' if 'up/kPoint/{:010}/moments'.format(1) in h5 else 'no'))
        print('  band derivatives : {}'.format('yes' if 'up/derivatives' in h5 else 'no'))
        print('  band curvatures  : {}'.format('yes' if 'up/curvatures' in h5 else 'no'))


  def outputDOS(self, plot, broadening=0.02):
    '''
    Print DOS/NOS
    '''

    import matplotlib.pyplot as plt

    with h5py.File(self.fname,'r') as h5:
      mu  = h5['.bands/mu'][()]
      spins = h5['.bands/ispin'][()]

      if spins==1:
        ene = h5['energies'][()]
        weights = h5['.kmesh/weights'][()]
        dosaxis, dos, nos = LRTCinput.calcDOS(ene, weights, gamma=broadening)
        del ene
      else:
        eneup = h5['up/energies'][()]
        weights = h5['.kmesh/weights'][()]
        dosaxisup, dosup, nosup = LRTCinput.calcDOS(eneup, weights, gamma=broadening)
        del eneup
        enedn = h5['dn/energies'][()]
        dosaxisdn, dosdn, nosdn = LRTCinput.calcDOS(enedn, weights, gamma=broadening)
        del enedn

      if plot:
        if spins==1:
          plt.plot(dosaxis, dos, color='black', lw=2, label='DOS')
          plt.legend(loc='upper left')
          plt.ylabel('DOS [eV^-1]')
          plt.xlabel(r'$\varepsilon$ [eV]')

          plt.twinx()
          plt.ylabel('NOS')
          plt.plot(dosaxis, nos, color='gray', lw=2, label='NOS')
          plt.legend(loc='upper right')
        else:
          plt.plot(dosaxisup, dosup, color='blue', lw=2, label='DOS up')
          plt.plot(dosaxisdn, -dosdn, color='red', lw=2, label='DOS dn')
          ylim = max(np.max(dosup),np.max(dosdn)) * 1.1
          plt.ylim(-ylim,ylim)
          plt.legend(loc='upper left')
          plt.ylabel('DOS [eV^-1]')
          plt.xlabel(r'$\varepsilon$ [eV]')

          plt.twinx()
          plt.plot(dosaxisup, nosup, color='deepskyblue', lw=2, label='NOS up')
          plt.plot(dosaxisup, -nosdn, color='indigo', lw=2, label='NOS dn')
          ylim = max(np.max(nosup),np.max(nosdn)) * 1.1
          plt.ylim(-ylim,ylim)
          plt.legend(loc='upper right')
          plt.ylabel('NOS')
        plt.axvline(x=mu, ls='--', color='gray', lw=2)
        plt.xlabel('energy')
      else:
        if spins==1:
          np.savetxt(self.textpipe, np.hstack((dosaxis[:,None], dos[:,None], nos[:,None])), \
                     fmt='%25.15e %25.15e %25.15e', comments='', \
                     header='# energy [eV], DOS [eV^-1], NOS]')
        else:
          np.savetxt(self.textpipe, np.hstack((dosaxis[:,None], dosup[:,None], dosdn[:,None], nosup[:,None], nosdn[:,None])), \
                     fmt='%25.15e %25.15e %25.15e %25.15e %25.15e', comments='', \
                     header='#  energy [eV], DOSup [eV^-1], DOSdn [eV^-1], NOSup, NOSdn')

    print('') # empty line before next CLI input

  @staticmethod
  def calcDOS(energies, weights, windowsize=1.4, npoints=1000, gamma=0.02):
    # first we find the energy interval
    globmin = np.min(energies)
    globmax = np.max(energies)

    if windowsize < 1:
      windowsize = 1
    increase = (windowsize-1)/2.

    # extend it a bit outwards
    interval = globmax-globmin
    globmin -= interval*increase
    globmax += interval*increase

    dosaxis = np.linspace(globmin,globmax,npoints)

    dosresolved = es.ElectronicStructure.Lorentzian(dosaxis[None,None,:],energies[:,:,None],gamma)
    dos = np.sum(dosresolved*weights[:,None,None],axis=(0,1))

    nosresolved = es.ElectronicStructure.IntLorentzian(dosaxis[None,None,:],energies[:,:,None],gamma)
    nos = np.sum(nosresolved*weights[:,None,None],axis=(0,1))
    return dosaxis, dos, nos
