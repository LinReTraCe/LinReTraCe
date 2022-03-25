#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np

def calcDOS(energies, weights, windowsize=1.2, npoints=1001, gamma=0.02):
  '''
    calculate the density of states and number of states for a given energy array
    input:
      energies: array shape [nkp,nbands]
      weights:  array shape[nkp]
      windowsize: dos size w.r.t. energy interval
      npoints: number of discretized energies
      gamma: energy broadening in the lorentzian used [eV]
    output:
      dosaxis: energy axis shape [npoints]
      dos: density of states shape [npoints]
      nos: number of states shape [npoints]
  '''

  globmin = np.min(energies)
  globmax = np.max(energies)

  if windowsize < 1:
    windowsize = 1
  increase = (windowsize-1)/2.

  # extend it a bit outwards
  interval = globmax-globmin
  globmin -= interval*increase
  globmax += interval*increase

  dosaxis = np.linspace(globmin,globmax,npoints,endpoint=True)

  dosresolved = Lorentzian(dosaxis[None,None,:],energies[:,:,None],gamma)
  dos = np.sum(dosresolved*weights[:,None,None],axis=(0,1))

  nosresolved = IntLorentzian(dosaxis[None,None,:],energies[:,:,None],gamma)
  nos = np.sum(nosresolved*weights[:,None,None],axis=(0,1))

  return dosaxis, dos, nos

def Lorentzian(x,x0,gamma):
  '''
  spectral function of an electron with constant scattering rate gamma
  '''
  return 1/np.pi * gamma / ((x-x0)**2 + gamma**2)

def IntLorentzian(x,x0,gamma):
  '''
  integrated spectral function of an electron with constant scattering rate gamma
  the integration constant is set such that the lowest energy ( e to -inf)
  corresponds to 0 occupation
  '''
  return np.arctan((x-x0)/gamma) / np.pi + 0.5
