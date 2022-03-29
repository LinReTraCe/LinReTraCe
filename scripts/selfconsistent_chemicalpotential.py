#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from   scipy.special import digamma
import scipy.optimize


''' digamma occupation '''
def occ(e,gamma,mu,beta):
  return 0.5 + digamma(0.5 + beta/2.0/np.pi*(gamma - 1j*(e-mu))).imag/np.pi

''' LRTC energy file ... one spin '''
with h5py.File('tb-80-80-1-irr.hdf5','r') as h5:
  ene = h5['energies'][()] # nkp bands=2
  weights = h5['.kmesh/weights'][()] # nkp
  mudft   = h5['.bands/mu'][()]
  charge  = h5['.bands/charge'][()]
  gapped  = h5['.bands/bandgap/gapped'][()]
  if gapped:
    ene_cband = h5['/.bands/bandgap/ene_cband'][()]
    ene_vband = h5['/.bands/bandgap/ene_vband'][()]



def deviation(mu,ene,gamma,beta):
  nocc = np.sum(occ(ene,gamma,mu,beta) * weights[:,None])
  return charge - nocc

''' temperature / beta range '''
TT = np.linspace(5,300,100, endpoint=True)
kB = 8.61733034e-5
BB = 1./(kB * TT)

''' which gamma(eps) '''
options = ['const','eps**2','sc_eps**2']
style = options[2]

''' gamma(eps) = gamma0 + (eps-mu)**2 * gamma2 '''
mu_arr = []
gamma0 = 1e-4
gamma2 = 5e-5

if style == 'const':
  gamma = gamma0
  for beta in BB[::-1]:
    mubisec = scipy.optimize.bisect(deviation, a=-0.2, b=0.2, args=(ene,gamma,beta))
    print(beta, mubisec)
    mu_arr.append(mubisec)
elif style == 'eps**2': # gamma + eps**2 w.r.t. DFT mu
  gamma = gamma0 + (ene-mudft)**2 * gamma2
  for beta in BB[::-1]:
    mubisec = scipy.optimize.bisect(deviation, a=-0.2, b=0.2, args=(ene,gamma,beta))
    print(beta, mubisec)
    mu_arr.append(mubisec)
elif style == 'sc_eps**2': # self-consisnt solution with gamma + eps**2 to mu + feedback loop
  mixing = 0.5
  for beta in BB[::-1]:
    mu_bisec = [mudft-0.5,mudft] # to enter the first while step
    while (abs(mu_bisec[-1] - mu_bisec[-2]) > 1e-6):
      gamma = gamma0 + (ene-mu_bisec[-1])**2 * gamma2
      muroot = scipy.optimize.bisect(deviation, a=-0.2, b=0.2, args=(ene,gamma,beta))
      mu_bisec.append(muroot*mixing + mu_bisec[-1]*(1-mixing))
    print(beta, mu_bisec[-1], len(mu_bisec)-1)
    mu_arr.append(mu_bisec[-1])

mu_arr = np.array(mu_arr)
''' this format can be used as OldMuText '''
np.savetxt('mu_psi.txt', np.hstack((TT[:,None],mu_arr[::-1,None])))
plt.plot(TT,mu_arr[::-1])
plt.axhline(y=mudft,color='gray', lw=1, zorder=0)
if gapped:
  plt.axhline(y=ene_cband,color='black', lw=2, zorder=0)
  plt.axhline(y=ene_vband,color='black', lw=2, zorder=0)
plt.show()
