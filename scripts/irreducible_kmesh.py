#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


nsym = 4
symop = np.zeros((nsym,3,3), dtype=np.float64)
symop[0,:,:] = np.array([[+1, 0, 0], \
                         [ 0,+1, 0], \
                         [ 0, 0,+1]])
symop[1,:,:] = np.array([[-1, 0, 0], \
                         [ 0,-1, 0], \
                         [ 0, 0,+1]])
symop[2,:,:] = np.array([[-1, 0, 0], \
                         [ 0,+1, 0], \
                         [ 0, 0,+1]])
symop[3,:,:] = np.array([[+1, 0, 0], \
                         [ 0,-1, 0], \
                         [ 0, 0,+1]])
invsymop = np.linalg.inv(symop)


class IrreducibleMesh(object):
  def __init__(self, nkx=1,nky=1,nkz=1,shift=False):
    self.nkx = nkx
    self.nky = nky
    self.nkz = nkz
    self.nkp = self.nkx*self.nky*self.nkz

    logger.info('Received {} reducible kpoints.'.format(self.nkp))
    self.shift = False

    self.weightsum = 2.0

    self._kmeshx = np.linspace(0,1,self.nkx,endpoint=False)
    self._kmeshy = np.linspace(0,1,self.nky,endpoint=False)
    self._kmeshz = np.linspace(0,1,self.nkz,endpoint=False)

    if self.shift:
      self._kmeshshift = []
      for ik in [self.nkx,self.nky,self.nkz]:
        if ik > 1:
          self._kmeshshift.append(1./ik/2.)
        else:
          self._kmeshshift.append(0.0)
      self._kmeshshift = np.array(self._kmeshshift, dtype=np.float64)

    # the way these points are ordered is important for the indexing below
    kpoints = []
    for ikx in self._kmeshx:
      for iky in self._kmeshy:
        for ikz in self._kmeshz:
          kpoints.append([ikx,iky,ikz])
    kpoints = np.array(kpoints)
    if self.shift: kpoints += self._kmeshshift[None,:]

    unique  = np.ones((self.nkx*self.nky*self.nkz), dtype=np.int)
    mult    = np.zeros((self.nkx*self.nky*self.nkz), dtype=np.float64)
    irrk    = 0

    logger.info('Generating irreducible kpoints...')

    for ik in range(self.nkp):

      if unique[ik] == 0: continue # skip if we already went there via symmetry
      irrk += 1    # new point -> increase irreducible counter
      mult[ik] = 1 # reset multiplicity counter

      ''' generate all the symmetry related k-points in the Brillouin zone '''
      knew = np.einsum('nji,j->ni',symop,kpoints[ik,:])
      kmod = knew%1
      ''' in order to index properly and if shift is applied , shift back '''
      if self.shift:
        kmod -= self._kmeshshift
      ''' round to neareast integer '''
      kround = np.rint(kmod * np.array([self.nkx,self.nky,self.nkz])[None,:])
      ''' exact floating calculation '''
      kexact = kmod * np.array([self.nkx,self.nky,self.nkz])[None,:]
      ''' only use the values that transform properly on all three axes '''
      mask = np.all(np.isclose(kround,kexact),axis=1)
      ''' apply the mask to filter '''
      kmask = kround[mask]
      ''' get the hash index '''
      kindex = (kmask[:,2] + \
                kmask[:,1] * self.nkz + \
                kmask[:,0] * self.nkz * self.nky).astype(int)
      ''' remove the k-points connected via symmetry and increase the multiplicity accordingly '''
      for ikk in kindex:
        if ikk <= ik: continue
        if unique[ikk]:
          unique[ikk] = 0
          mult[ik] += 1

    self.nkp = irrk
    self.kpoints = kpoints[unique>0]
    self.multiplicity = mult[unique>0]
    self.weights      = self.weightsum * self.multiplicity / np.sum(self.multiplicity)
    logger.info('Generated irreducible kmesh with {} irreducible kpoints.'.format(self.nkp))
    logger.info('kpoints:\n{}'.format(self.kpoints))

if __name__ == '__main__':
  logger.info('Testing irreducible mesh generation')
  mesh = IrreducibleMesh(50,20,1,shift=False)
