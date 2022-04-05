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
from structure.dos import calcDOS

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
        h5.attrs['identifier'].encode('utf-8') == 'LRTCinput'
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
      print('  irreducible: {}'.format(h5['.kmesh/irreducible'][()]))
      print('  kx ky kz: {} {} {}'.format(h5['.kmesh/nkx'][()],h5['.kmesh/nky'][()],h5['.kmesh/nkz'][()]))

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
        print('  intra optical elements   : {}'.format('yes' if 'momentsDiagonal' in h5 else 'no'))
        print('  inter optical elements   : {}'.format('yes' if 'kPoint/{:010}/moments'.format(1) in h5 else 'no'))
        print('  intra magnetic elements  : {}'.format('yes' if 'momentsDiagonalBfield' in h5 else 'no'))
        print('  inter magnetic elements  : {}'.format('yes' if 'kPoint/{:010}/momentsBfield'.format(1) in h5 else 'no'))
      else:
        print('  intra optical elements   : {}'.format('yes' if 'up/momentsDiagonal' in h5 else 'no'))
        print('  inter optical elements   : {}'.format('yes' if 'up/kPoint/{:010}/moments'.format(1) in h5 else 'no'))
        print('  intra magnetic elements  : {}'.format('yes' if 'up/momentsDiagonalBfield' in h5 else 'no'))
        print('  inter magnetic elements  : {}'.format('yes' if 'up/kPoint/{:010}/momentsBfield'.format(1) in h5 else 'no'))


  def outputDOS(self, plot, broadening=0.02, npoints=1001):
    '''
    Print DOS/NOS
    '''

    import matplotlib.pyplot as plt

    broadening = float(broadening)
    npoints = int(npoints)
    if broadening < 0: raise IOError('DOS broadening must be > 0')
    if npoints < 1: raise IOError('DOS points must be >= 1')

    with h5py.File(self.fname,'r') as h5:
      mu  = h5['.bands/mu'][()]
      spins = h5['.bands/ispin'][()]

      if spins==1:
        ene = h5['energies'][()]
        weights = h5['.kmesh/weights'][()]
        dosaxis, dos, nos = calcDOS(ene, weights, npoints=npoints, gamma=broadening, windowsize=1.1)
        del ene
      else:
        eneup = h5['up/energies'][()]
        weights = h5['.kmesh/weights'][()]
        dosaxisup, dosup, nosup = calcDOS(eneup, weights, npoints=npoints, gamma=broadening, windowsize=1.1)
        del eneup
        enedn = h5['dn/energies'][()]
        dosaxisdn, dosdn, nosdn = calcDOS(enedn, weights, npoints=npoints, gamma=broadening, windowsize=1.1)
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
