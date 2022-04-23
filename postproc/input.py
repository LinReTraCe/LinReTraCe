#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import argparse
import warnings
import logging
logger = logging.getLogger(__name__)

import numpy as np
with warnings.catch_warnings():
  warnings.filterwarnings("ignore",category=FutureWarning)
  import h5py
import ase

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
    print('{:<12}  {}'.format('path', 'band structure along high symmetry path'))
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
        print('  gapped: {}'.format(str(gapped)))
        if gapped:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/cband'][()]))
      else:
        gappedup = h5['.bands/bandgap/up/gapped'][()]
        print('  spin up gapped: {}'.format(str(gappedup)))
        if gappedup:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/up/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/up/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/up/cband'][()]))
        gappeddn = h5['.bands/bandgap/dn/gapped'][()]
        print('  spin dn gapped: {}'.format(str(gappeddn)))
        if gappeddn:
          print('    bandgap: {0:.4f} eV'.format(h5['.bands/bandgap/dn/gapsize'][()]))
          print('    valence    band: {}'.format(h5['.bands/bandgap/dn/vband'][()]))
          print('    conduction band: {}'.format(h5['.bands/bandgap/dn/cband'][()]))
      print('  chemical potential: {0:.6f} eV'.format(h5['.bands/mu'][()]))

      print('\nUNIT-CELL')
      print('  lattice (rows) [Ang]:\n{}'.format(h5['.unitcell/rvec'][()]))
      print('  volume: {0:.4f} [Ang^3]'.format(h5['.unitcell/volume'][()]))
      print('  orthogonal: {}'.format(h5['.unitcell/ortho'][()]))
      print('  symmetries: {}'.format(h5['.unitcell/nsym'][()]))

      print('\nOPTICAL DATA')
      if spins==1:
        print('  intra optical elements   : {}'.format(str('momentsDiagonal' in h5)))
        print('  inter optical elements   : {}'.format(str('kPoint/{:010}/moments'.format(1) in h5)))
        print('  intra magnetic elements  : {}'.format(str('momentsDiagonalBfield' in h5)))
        print('  inter magnetic elements  : {}'.format(str('kPoint/{:010}/momentsBfield'.format(1) in h5)))
      else:
        print('  intra optical elements   : {}'.format(str('up/momentsDiagonal' in h5)))
        print('  inter optical elements   : {}'.format(str('up/kPoint/{:010}/moments'.format(1) in h5)))
        print('  intra magnetic elements  : {}'.format(str('up/momentsDiagonalBfield' in h5)))
        print('  inter magnetic elements  : {}'.format(str('up/kPoint/{:010}/momentsBfield'.format(1) in h5)))


  def outputDOS(self, plot, broadening=0.02, npoints=1001):
    '''
    Print DOS/NOS
    '''

    if plot:
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

  def outputPath(self, plot, pathstring=None):
    '''
      output or plot high-symmetry path
      if pathstring is None: return all avaiable high-symmetry points
      if pathstring is provided we check for consistency and generate the result.
      since for irreducible calculations we only have limited k-points available
      we employ symmetries to map back into the avaible points
      the path length is automatically generated via the greates common denominator
      to 'hit' the largest number of possible k-point in the BZ.
    '''

    if plot:
      import matplotlib.pyplot as plt

    with h5py.File(self.fname,'r') as h5:
      spins    = h5['.bands/ispin'][()]
      rvec     = h5['.unitcell/rvec'][()]
      kvec     = h5['.unitcell/kvec'][()]
      nkx      = h5['.kmesh/nkx'][()]
      nky      = h5['.kmesh/nky'][()]
      nkz      = h5['.kmesh/nkz'][()]
      kpoints  = h5['.kmesh/points'][()] % 1
      dims     = h5['.unitcell/dims'][()]
      nsym     = h5['.unitcell/nsym'][()]
      symop    = h5['.unitcell/symop'][()]
      mudft    = h5['.bands/mu'][()]

      if spins == 1:
        energies = h5['energies'][()]
      else:
        energiesup = h5['up/energies'][()]
        energiesdn = h5['dn/energies'][()]
        energies = np.zeros((2,energiesup.shape[0],energiesup.shape[1]), dtype=np.float64)
        energies[0,...] = energiesup
        energies[1,...] = energiesdn

    ''' help arrays '''
    nk = np.array([nkx,nky,nkz], dtype=int)
    nkindex = np.array([nkx*nky*nkz,nky*nkz,nkz])
    dims = list(dims.astype(int))

    ''' create "hashing" index '''
    index = []
    for ikpt in kpoints:
      index.append(np.sum(np.around(ikpt*nkindex)))
    index = np.array(index,dtype=int)


    ''' retrieve the special points via ASE '''
    cell = ase.cell.Cell(rvec)
    special = cell.bandpath('', pbc=dims)
    fullspecial = []
    for descr, coord in special.special_points.items():
      fullspecial.append(descr)
      logger.info('special point: {} {}'.format(descr, coord))
    fullspecial = ''.join(fullspecial)
    if pathstring is None:
      logger.info('Provide k-path in form of trailing argument, e.g. -- path {}'.format(fullspecial))
      return

    ''' check for consistency '''
    for ispecial in pathstring:
      if ispecial not in fullspecial:
        logger.critical('Provided special point not available.')
        return

    substrings = [pathstring[i:i+2] for i in range(len(pathstring)-1)]
    k_distances = []
    kptsindex  = []
    kptsplotx  = []
    xticks     = []

    ''' iterate over all the sub strings that connect 2 special points '''
    for istring, stringpath in enumerate(substrings):

      special1, special2 = special.special_points[stringpath[0]], special.special_points[stringpath[1]]
      kspecial1 = np.einsum('ji,j->i',kvec,special1)
      kspecial2 = np.einsum('ji,j->i',kvec,special2)
      k_distances.append(np.sqrt(np.sum((kspecial1-kspecial2)**2)))

      ''' get optimal spacing '''
      distance = np.abs(special1-special2) * nk
      npoints = np.gcd.reduce(distance[dims].astype(int))+1

      ''' generate path with optimal spacing and dimension selection '''
      path       = cell.bandpath(stringpath, pbc=dims, npoints=npoints)
      kpts       = np.array(path.kpts) % 1 # back to BZ
      npts       = len(kpts)

      if istring == 0:
        xoffset = 0
      else:
        xoffset += k_distances[istring-1]

      xticks.append(xoffset)
      for i, ikpt in enumerate(kpts):
        if istring > 0 and i == 0: continue # to avoid double points from path overlap
        ''' only continue with point on the grid '''
        if not np.allclose((ikpt*nk),np.around(ikpt*nk), atol=0.001):
          continue

        ''' map back into available points '''
        symkpt = np.einsum('nji,j->ni',symop,ikpt) % 1
        for isym in range(symkpt.shape[0]):
          if np.allclose((symkpt[isym]*nk),np.around(symkpt[isym]*nk)):
            symindex = np.sum(np.around(symkpt[isym]*nkindex))
            if symindex in index:
              kptsindex.append(symindex)
              kptsplotx.append(xoffset + i/(float(npts)-1) * k_distances[istring])
              break

    ''' last point of the ticks '''
    xticks.append(xoffset+k_distances[-1])
    kptsindex  = np.array(kptsindex, dtype=int)
    arrayindex = np.where(index[None,:]==kptsindex[:,None])[1]

    ''' get energies and plot / output to stdout '''
    for ispin in range(spins):
      bandpath = []
      for iarray in arrayindex:
        if spins==1:
          bandpath.append(energies[iarray,:])
        else:
          bandpath.append(energies[ispin,iarray,:])
      bandpath = np.array(bandpath, dtype=np.float64)

      if plot:
        if ispin == 0:
          style = {'lw':3, 'ls':'-', 'color':'k'}
        else:
          style = {'lw':3, 'ls':'-', 'color':'red'}
        for iband in range(bandpath.shape[1]):
          plt.plot(kptsplotx,bandpath[:,iband], **style)
        plt.xticks(xticks, [str(i) for i in pathstring])
        for xval in xticks:
          plt.axvline(x=xval, color='gray', lw=0.5, ls='-', zorder=-1)
      else:
        np.savetxt(self.textpipe, np.hstack((np.array(kptsplotx)[:,None],bandpath)), \
                   comments='', \
                   header='# k [Ang^(-1)], bands [eV] - {}'.format('' if spins==1 else ['up','dn'][ispin]))
    if plot:
      plt.ylabel(r'$\varepsilon(\mathbf{k})$ [eV]', fontsize=12)
      plt.axhline(y=mudft, color='gray', lw=1, ls='--', zorder=-1)
      plt.show()
