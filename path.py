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
k_distances = []

kptsindex  = []
kptsplotx  = []
xticks     = []
for istring, stringpath in enumerate(substrings):

  ''' get optimal spacing '''
  special1, special2 = special.special_points[stringpath[0]], special.special_points[stringpath[1]]
  kspecial1 = np.einsum('ji,j->i',kvec,special1)
  kspecial2 = np.einsum('ji,j->i',kvec,special2)


  print(special1, ' --> ', special2)
  k_distances.append(np.sqrt(np.sum((kspecial1-kspecial2)**2)))

  distance = np.abs(special1-special2) * nk
  npoints = np.gcd.reduce(distance[dims].astype(int))+1

  path       = cell.bandpath(stringpath, pbc=dims, npoints=npoints)
  kpts       = np.array(path.kpts) % 1
  npts       = len(kpts)

  if istring == 0:
    xoffset = 0
  else:
    xoffset += k_distances[istring-1]

  xticks.append(xoffset)

  print('points on path: {}'.format(npts))
  for i, ikpt in enumerate(kpts):
    ''' only continue with point on the grid '''
    if not np.allclose((ikpt*nk),np.around(ikpt*nk), atol=0.001):
      continue

    symkpt = np.einsum('nji,j->ni',symop,ikpt) % 1
    for isym in range(symkpt.shape[0]):
      if np.allclose((symkpt[isym]*nk),np.around(symkpt[isym]*nk)):
        symindex = np.sum(np.around(symkpt[isym]*nkindex))
        if symindex in index:
          kptsindex.append(symindex)
          kptsplotx.append(xoffset + i/float(npts) * k_distances[istring])
          break


''' last point of the ticks '''
xoffset += k_distances[istring-1]
xticks.append(xoffset)

kptsindex  = np.array(kptsindex, dtype=int)
arrayindex = np.where(index[None,:]==kptsindex[:,None])[1]
bandpath = []
for iarray in arrayindex:
  bandpath.append(energies[iarray,:])
bandpath = np.array(bandpath, dtype=np.float64)

style = {'lw':3, 'ls':'-', 'color':'k'}
for iband in range(bandpath.shape[1]):
  plt.plot(kptsplotx,bandpath[:,iband], **style)

plt.xticks(xticks, [str(i) for i in pathstring])
plt.show()
