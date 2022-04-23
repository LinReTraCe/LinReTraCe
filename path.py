#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import pprint
import ase

# desired results:
# lprint <LRTCinput file> path
# -> special point + index
# lprint (-p) <LRTCinput file> path GKMG

with h5py.File('lrtc-tb-300-300-1-irr.hdf5','r') as h5:
  rvec = h5['.unitcell/rvec'][()]
  kvec = h5['.unitcell/kvec'][()]
  nkx     = h5['.kmesh/nkx'][()]
  nky     = h5['.kmesh/nky'][()]
  nkz     = h5['.kmesh/nkz'][()]
  kpoints = h5['.kmesh/points'][()] % 1
  dims    = h5['.unitcell/dims'][()]
  nsym    = h5['.unitcell/nsym'][()]
  symop   = h5['.unitcell/symop'][()]
  energies = h5['energies'][()]

nk = np.array([nkx,nky,nkz], dtype=int)
nkindex = np.array([nkx*nky*nkz,nky*nkz,nkz])
dims = list(dims.astype(int))

''' create "hashing" index '''
index = []
for ikpt in kpoints:
  index.append(np.sum(np.around(ikpt*nkindex)))
index = np.array(index,dtype=int)


cell = ase.cell.Cell(rvec)
special= cell.bandpath('', pbc=dims)
for descr, coord in special.special_points.items():
  print('special points: ', descr, coord)

pathstring = 'GKMG'
substrings = [ pathstring[i:i+2] for i in range(len(pathstring)-1) ]

kptsindex  = []
kptsplotx  = []
xticks = []
xticklabels = []
for istring, stringpath in enumerate(substrings):
  xticks.append(len(kptsindex)-1)
  if xticks[-1] < 0: xticks[-1] = 0
  path       = cell.bandpath(stringpath, pbc=dims, npoints=101)
  print(path)
  kpts       = np.array(path.kpts) % 1

  xbegin = np.sqrt(np.sum(kpts[0]**2))
  if istring == 0:
    xoffset = 0
  else:
    xoffset  += np.sqrt(np.sum(kpts[-1]**2))

  if istring != (len(substrings)-1):
    kpts = kpts[:-1,:]


  print('points on path: {}'.format(len(kpts)))

  kptsfilter = []
  ''' corresponding to an actual grid position '''
  for ikpt in kpts:
    if np.allclose((ikpt*nk),np.around(ikpt*nk), atol=0.001):
      kptsfilter.append(ikpt)

  before = len(kptsindex)
  for ikpt in kptsfilter:
    ikptindex = np.sum(np.around(ikpt*nkindex))
    if ikptindex in index:
      kptsindex.append(ikptindex)
    else:
      symkpt = np.einsum('nji,j->ni',symop,ikpt) % 1
      for isym in range(symkpt.shape[0]):
        if np.allclose((symkpt[isym]*nk),np.around(symkpt[isym]*nk), rtol = 0.01, atol = 0.01):
          symindex = np.sum(np.around(symkpt[isym]*nkindex))
          if symindex in index:
            kptsindex.append(symindex)
            kptsplotx.append(xoffset + np.sqrt(np.sum((np.array(ikpt)-kpts[0])**2)))
            break
  after = len(kptsindex)

  print('points added {}'.format(after-before))

xticks.append(len(kptsindex)-1)

kptsindex  = np.array(kptsindex, dtype=int)
arrayindex = np.where(index[None,:]==kptsindex[:,None])[1]
bandpath = []
for iarray in arrayindex:
  bandpath.append(energies[iarray,:])
bandpath = np.array(bandpath, dtype=np.float64)

style = {'lw':3, 'ls':'-', 'color':'k'}
for iband in range(bandpath.shape[1]):
  plt.plot(bandpath[:,iband], **style)

plt.xticks(xticks, [str(i) for i in pathstring])
plt.show()
