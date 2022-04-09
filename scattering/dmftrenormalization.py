#! /usr/bin/env python

from __future__ import print_function, division, absolute_import

import numpy as np
import h5py


class DmftRenormalization(object):
  '''
  Class to load DMFT renormalization data
  We get a renormalization factor, a pole shift (given by Re Selfenergy(w=0))
  and a scattering rate (-Im Selfenergy(w=0)  in the local Wannier basis.
  After backtransforming to the Kohn-Sham basis we get band- and k-dependent
  data which we load here.
  '''

  def __init__(self, filename):
    self.filename      = filename
    self.nkp           = None
    self.efermi        = None
    self.energyBandMax = None
    self.multiplicity  = []
    self.bands         = []

  def __repr__(self):
    return ('DMFTrenormalization(filename={0.filename!r})'.format(self))

  def load(self, cutoff=1e-7):
    '''
    Load the data (Jan's format) into lists;
    truncate the data and cast into arrays.
    '''

    myZqp         = []
    myEnergyshift = []
    myGamma       = []
    myEnergy      = []

    with open(self.filename, 'r') as dmft:
      firstline = dmft.readline().split()
      self.efermi = float(firstline[0])
      self.nkp    = int(firstline[1])

      for ikp in range(self.nkp):
        kpointHeader = dmft.readline().split()
        self.bands.append(int(kpointHeader[1]))
        self.multiplicity.append(float(kpointHeader[2]))

        myZqp.append([])
        myEnergyshift.append([])
        myGamma.append([])
        myEnergy.append([])

        for _ in range(self.bands[ikp]):
          data = np.array(dmft.readline().split(), dtype=np.float64)
          myEnergy[ikp].append(data[0])
          myEnergyshift[ikp].append(data[1])
          myZqp[ikp].append(data[2])
          if (data[3] < cutoff):
            data[3] = cutoff
          myGamma[ikp].append(data[3])

      self.bands = np.array(self.bands, dtype=np.int)
      self.multiplicity = np.array(self.multiplicity, dtype=np.int)
      self.energyBandMax = self.bands.min() # only take common bands

      self.zqp         = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)
      self.energyshift = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)
      self.gamma       = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)
      self.energy      = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)

      for ikp in range(self.nkp):
        self.zqp[ikp,:]         = np.array(myZqp[ikp][:self.energyBandMax])
        self.energyshift[ikp,:] = np.array(myEnergyshift[ikp][:self.energyBandMax])
        self.gamma[ikp,:]       = np.array(myGamma[ikp][:self.energyBandMax])
        self.energy[ikp,:]      = np.array(myEnergy[ikp][:self.energyBandMax])

# # raw data and empty list which we will fill up
# dfile = 'fort.111'
# dlist = []

# # load everything into a convenient array
# # we can do this because the expected size is small
# f = open(dfile, 'r')
# farray = f.readlines()
# print('Length of file: ',len(farray))
# print()

# efermi = float(farray[0]) # fermi energy
# ikp, emin, emax = [int(x) for x in farray[1].split()]
# globmax = emax

# skip = 2 # first two lines
# while True:
#   f.seek(0) # reset
#   print(skip,emax)
#   dlist.append(np.genfromtxt(f, usecols=(0,1,2), dtype=np.float64, skip_header=skip, max_rows=emax).T)
#   if skip+emax < len(farray):
#     skip += emax+1 # all the previous plus the old data + 1 (ikp information)
#     ikp, emin, emax = [int(x) for x in farray[skip-1].split()] # -1 since we need exactly the ikp information
#     if emax > globmax: # check for global maximum
#       globmax = emax
#   else:
#     f.close()
#     break

# nkp = len(dlist)
# print('\nnkp: ', nkp)
# print('global emax: ', globmax)

# # set up full array with np.nans
# dfinal = np.full((3, nkp, globmax), np.nan)

# # fill in the array, missing values compared to the global maximum
# # are represented as np.nan
# for i in xrange(nkp):
#   dfinal[:,i,:dlist[i].shape[-1]] = dlist[i]

# # write out to h5py
# with h5py.File('result.hdf5', 'w') as f:
#   f['energies'] = dfinal[0,...]
#   f['zqp'] = dfinal[1,...]
#   f['gamma'] = dfinal[2,...]

# print('Done.')

if __name__ == '__main__':
  mydmft = DMFTrenormalization('../INTERFACE/Ce3Bi4Pt3/beta400/input_lrt.dat')
  mydmft.load()
  print(mydmft.bands)
  print(mydmft.multiplicity)
  print(mydmft.zqp.shape)
