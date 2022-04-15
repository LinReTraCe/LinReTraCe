#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import logging
import argparse

import numpy as np
import spglib

def parse_args(args=None):
  parser = argparse.ArgumentParser(
    description='Argument parser for symmetrization of tb data ')
  # mandatory
  parser.add_argument('tb_file', type=str, help='tight binding file')
  parser.add_argument('rvec_file', type=str, help='primitive lattice vectors file')
  parser.add_argument('atom_file', type=str, help='atomic positions file')
  parser.add_argument('--debug', help=argparse.SUPPRESS, default=False, action='store_true')
  return parser.parse_args(args)

if __name__ == '__main__':
  logging.basicConfig()
  logger = logging.getLogger()
  args = parse_args()
  logger.setLevel(logging.DEBUG if args.debug else logging.INFO)

  tbdata = np.genfromtxt(args.tb_file)
  if len(tbdata.shape) == 1:
    tbdata = tbdata[None,:]

  if tbdata.shape[1] == 6:
    imaghopping = False
  else:
    imaghopping = True

  rvecdata = np.genfromtxt(args.rvec_file)

  atomdata = np.genfromtxt(args.atom_file)
  if len(atomdata.shape) == 1:
    atomdata = atomdata[None,:]

  ''' Define spglib cell object '''
  lattice   = rvecdata
  positions = []
  numbers   = []
  for i in range(atomdata.shape[0]):
    numbers.append(atomdata[i,0])
    positions.append(atomdata[i,1:])

  cell = (lattice, positions, numbers)

  symmetry = spglib.get_symmetry(cell, symprec=1e-5)
  symsfull = symmetry['rotations']

  ''' Filter out weird symmetry operations due to non-standardized rvectors '''
  non_standard = False
  symop = []
  for ii in np.arange(symsfull.shape[0]):
    isym = symsfull[ii]
    if abs(np.linalg.det(isym)) != 1:
      non_standard = True
      isym[isym>1] = 0
      isym[isym<(-1)] = 0
      if abs(np.linalg.det(isym)) != 1:
        raise ValueError('Non-stanardized unit cell resulted in matrix with invalid determinant')
    symop.append(isym)
  symop = np.array(symop)
  nsym = symop.shape[0]
  logger.info('Spglib: Found {} symmetry operations'.format(nsym))
  logger.debug('Symmetry operation: \n{}'.format(symop))


  symmetrized_tbdata = []
  for itb1 in range(tbdata.shape[0]):
    rvec1  = tbdata[itb1,:3].astype(int)
    band1_1  = int(tbdata[itb1,3]) - 1 # band identifier
    band1_2  = int(tbdata[itb1,4]) - 1 # band identifier
    hop1real = tbdata[itb1,5]
    if imaghopping:
      hop1imag = tbdata[itb1,6]
      hop1 = hop1real + 1j*hop1imag
      symmetrized_tbdata.append([*rvec1,band1_1,band1_2,hop1real,hop1imag])
    else:
      hop1 = hop1real
      symmetrized_tbdata.append([*rvec1,band1_1,band1_2,hop1real])

    rvecsym = np.einsum('nij,j->ni',symop,rvec1)

    for isym in range(nsym):
      rvec_transformed = rvecsym[isym]
      for itb2 in range(tbdata.shape[0]):
        band2_1  = int(tbdata[itb1,3]) - 1 # band identifier
        band2_2  = int(tbdata[itb1,4]) - 1 # band identifier
        if band1_1 != band2_1 or band2_2 != band2_2: continue

        rvec2  = tbdata[itb2,:3].astype(int)
        hop2real   = tbdata[itb2,5]
        if imaghopping:
          hop2imag = tbdata[itb2,6]
          hop2 = hop2real + 1j*hop2imag
        else:
          hop2 = hop2real

        if np.allclose(rvec_transformed,rvec2) and np.abs(hop1-hop2) < 1e-6:
          break

      else: # if we did not break
        if imaghopping:
          data = [*rvec_transformed,band2_1,band2_2,hop2real,hop2imag]
        else:
          data = [*rvec_transformed,band2_1,band2_2,hop2real]

        if data not in symmetrized_tbdata:
          symmetrized_tbdata.append(data)

  print('\nSYMMETRIZED DATA:\n-----------------------')
  for line in symmetrized_tbdata:
    print('{:2} {:2} {:2}   {:2} {:2}    {} {}'.format(*line[:6],line[6] if imaghopping else ''))
