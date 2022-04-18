#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import logging
import argparse

import numpy as np
import spglib

def parse_args(args=None):
  parser = argparse.ArgumentParser(
    description='Argument parser for symmetrization of tb data ', \
    epilog=''' ./symmetrize_tightbinding.py <tb_file> xy ''')
  # mandatory
  parser.add_argument('tb_file',    type=str, help='tight binding file')
  parser.add_argument('dimensions', type=str, help='dimensions to apply symmetries to "xy" "xyz"')
  parser.add_argument('--debug', help=argparse.SUPPRESS, default=False, action='store_true')
  return parser.parse_args(args)

if __name__ == '__main__':
  logging.basicConfig()
  logger = logging.getLogger()
  args = parse_args()
  logger.setLevel(logging.DEBUG if args.debug else logging.INFO)

  dims = np.array(['x' in args.dimensions,'y' in args.dimensions, 'z' in args.dimensions], dtype=np.bool)

  atoms = []
  with open(str(args.tb_file),'r') as data:
    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if line.startswith('begin atoms'): break
      else:
        raise IOError('ltb could not detect "begin atoms" block')

    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if len(line)==0 or line.startswith('#'): continue
        if line.startswith('end atoms'): break
      else:
        raise IOError('ltb is not at the "end atoms" marker after reading'.format(str(args.tb_file)))

      comment = line.find('#')
      if comment != -1: line = line[:comment]
      atom_id, rx, ry, rz = [float(i) for i in line.split()]
      if abs(np.rint(atom_id) - atom_id) > 1e-6: raise IOError('Invalid atom description: {}'.format(line))
      if rx < 0 or rx > 1: raise IOError('Invalid atom description: {}'.format(line))
      if ry < 0 or ry > 1: raise IOError('Invalid atom description: {}'.format(line))
      if rz < 0 or rz > 1: raise IOError('Invalid atom description: {}'.format(line))
      atoms.append([atom_id, rx, ry, rz])
  atoms = np.array(atoms, dtype=np.float64)
  logger.debug('   atoms: \n{}'.format(atoms))

  rvec = []
  with open(str(args.tb_file),'r') as data:
    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if line.startswith('begin real_lattice'): break
      else:
        raise IOError('ltb could not detect "begin real_lattice" block')

    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if len(line)==0 or line.startswith('#'): continue
        if line.startswith('end real_lattice'): break
        comment = line.find('#')
        if comment != -1: line = line[:comment]
        rx, ry, rz = [float(i) for i in line.split()]
        rvec.append([rx,ry,rz])
      else:
        raise IOError('ltb is not at the "end real_lattice" marker after reading'.format(str(args.tb_file)))

  rvec = np.array(rvec, dtype=np.float64)
  if rvec.shape != (3,3):
    raise IOError('Provided real_lattice is not a properly defined 3x3 matrix')
  logger.debug('   real_lattice (rows): \n{}'.format(rvec))

  tbdata = []
  bandmin = bandmax = None

  with open(str(args.tb_file),'r') as data:
    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if line.startswith('begin hopping'): break
      else:
        raise IOError('ltb could not detect "begin hopping" block')

    while True:
      line = data.readline()
      if line:
        line = line.strip()
        if len(line)==0 or line.startswith('#'): continue
        if line.startswith('end hopping'): break

        comment = line.find('#')
        if comment != -1: line = line[:comment]
        linedata = line.split()
        lenline  = len(linedata)
        if lenline==6 or lenline==7:
          tbdata.append(np.array([float(i) for i in linedata]))
        else:
          raise IOError('Provided hopping paramater line is invald: {}'.format(line))
      else:
        raise IOError('ltb is not at the "end hopping" marker after reading'.format(str(tbdata)))

  ''' we do not transform the list of lists into an array since the different lines
      could have different lengths (imaginary part optional) '''

  ''' detect inter band transitions and minimum and maximum band entry '''
  inter     = False
  bandcheck = set()
  for i in range(len(tbdata)):
    locmin = min(tbdata[i][3:5].astype(int))
    locmax = max(tbdata[i][3:5].astype(int))
    if bandmin is None:
      bandmin = locmin
      bandmax = locmax
    else:
      if locmin < bandmin: bandmin = locmin
      if locmax > bandmax: bandmax = locmax
    if locmin != locmax: inter = True
    bandcheck.update({locmin,locmax})
  if bandmin != 1:
    raise IOError('Tight binding parameter set must start at band 1')
  if len(bandcheck) != bandmax:
    raise IOError('Tight binding parameter set does not contain the full band range')



  ''' Define spglib cell object '''
  lattice   = rvec
  positions = []
  numbers   = []
  for i in range(atoms.shape[0]):
    numbers.append(atoms[i,0])
    positions.append(atoms[i,1:])

  cell = (lattice, positions, numbers)

  symmetry = spglib.get_symmetry(cell, symprec=1e-5)
  symsfull = symmetry['rotations']
  logger.debug('Symmetry operation: \n{}'.format(symsfull))

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
    to_add = True
    for i, dim in enumerate(dims):
      if dim: continue
      for j in range(3):
        if i==j:
          if isym[i,j] != 1: to_add = False
        else:
          if isym[j,i] != 0: to_add = False
          if isym[i,j] != 0: to_add = False # redundant I think
    if to_add: symop.append(isym)
  symop = np.array(symop)
  nsym = symop.shape[0]
  logger.info('Spglib: Found {} symmetry operations'.format(nsym))
  logger.debug('Symmetry operation: \n{}'.format(symop))


  symmetrized_tbdata = []
  for itb1 in range(len(tbdata)):
    rvec1  = tbdata[itb1][:3].astype(int)
    incell = np.allclose(rvec1,np.zeros((3,)))
    band1_1  = int(tbdata[itb1][3]) # band identifier
    band1_2  = int(tbdata[itb1][4]) # band identifier
    hop1real = tbdata[itb1][5]
    try:
      hop1imag = tbdata[itb1][6]
      hop1 = hop1real + 1j*hop1imag
      symmetrized_tbdata.append([*rvec1,band1_1,band1_2,hop1real,hop1imag])
    except:
      hop1 = hop1real
      symmetrized_tbdata.append([*rvec1,band1_1,band1_2,hop1real])

    rvecsym = np.einsum('nij,j->ni',symop,rvec1)

    if incell:
      if band1_1 != band1_2:
        try:
          data = [*rvec1,band1_2,band1_1,hop1real,-hop1imag]
        except:
          data = [*rvec1,band1_2,band1_1,hop1real]
        if data not in symmetrized_tbdata:
          symmetrized_tbdata.append(data)
    else:
      for isym in range(nsym):
        rvec_transformed = rvecsym[isym]
        for itb2 in range(len(tbdata)):
          band2_1  = int(tbdata[itb2][3]) # band identifier
          band2_2  = int(tbdata[itb2][4]) # band identifier
          if band1_1 == band2_1 and band2_2 == band2_2: ## 12 == 12
            pass
          elif band1_1 == band2_2 and band1_2 == band2_1: ## 12 == 21
            pass
          else:
            continue

          rvec2  = tbdata[itb2][:3].astype(int)
          hop2real   = tbdata[itb2][5]
          try:
            hop2imag = tbdata[itb2][6]
            if np.allclose(rvec_transformed,rvec2) and np.abs(hop1real-hop2real) < 1e-6:
              if np.allclose(rvec_transformed,rvec2) and np.abs(hop1imag-hop2imag) < 1e-6:
                break
          except:
            if np.allclose(rvec_transformed,rvec2) and np.abs(hop1real-hop2real) < 1e-6:
              break

        else: # if we did not break
          try:
            data = [*rvec_transformed,band2_1,band2_2,hop1real,hop1imag]
          except:
            data = [*rvec_transformed,band2_1,band2_2,hop1real]

          if data not in symmetrized_tbdata:
            symmetrized_tbdata.append(data)

  print('\nSYMMETRIZED DATA:\n-----------------------')
  for line in symmetrized_tbdata:
    try:
      print('{:2} {:2} {:2}   {:2} {:2}    {} {}'.format(*line[:7]))
    except:
      print('{:2} {:2} {:2}   {:2} {:2}    {}'.format(*line[:6]))

  if bandmax > 1:
    logger.critical('\n\nInter-band hopping across unit cells cannot be safely generated this way.' +
                    '\nFilter out unwanted hopping parameters!')
