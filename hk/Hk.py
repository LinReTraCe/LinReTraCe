from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np
import json
import sympy
from sympy import sympify

from structure.auxiliary import levicivita
from structure.auxiliary import progressBar
from structure.model     import Model
from hk.hk_matrix_elements import *

class hk(Model):
  '''
  Symbolic K-space Hamiltonian class which calculates the matrix elements symbolically using SymPy package

  Made to be in the same form as the tb class so a lrtc energy-file of the same form can be generated for best 
  integration with the rest of the LinReTraCe code 
  '''
  def __init__(self, nkx=1, nky=1, nkz=1, kshift=False,write_to_file=False):
    super(hk, self).__init__(nkx,nky,nkz)
    self.irreducible = False  # generate reducible k-point grid 
    self.kshift      = kshift       # shift by half a k-point to avoid Gamma point
    self.opticfull= True
    self.write_to_file=write_to_file

  def _read_hk_file(self): #reads in the coeffients saved in the Hamiltonian.json file and calulcates matrix elements 
      
    Hamk,b1,b2,b3,name=load_hk(self.hk_file)

    self.n_bands=int(np.sqrt(len(Hamk))) #defines the number of bands in the system 

    self.name=name

    self.kvec=np.zeros((3,3))
    self.rvec=np.zeros((3,3))

    #makes sure that the vectors are read in as numerical values not as Sympy values

    for i in range(3):
      self.kvec[0,i]=N(b1[i])
      self.kvec[1,i]=N(b2[i])
      self.kvec[2,i]=N(b3[i])

    Vk=np.abs(np.dot(np.cross(self.kvec[0,:],self.kvec[1,:]),self.kvec[2,:])) #volume of unit cell in reciprocal space 

    self.rvec[0,:] = np.cross(self.kvec[1,:],self.kvec[2,:]) / Vk 
    self.rvec[1,:] = np.cross(self.kvec[2,:],self.kvec[0,:]) / Vk
    self.rvec[2,:] = np.cross(self.kvec[0,:],self.kvec[1,:]) / Vk
    self.rvec *= 2*np.pi
 
    if self.nkz==1: #if there is only 1 k-point in z-direc treats the system as 2D 
      self.vol = np.abs(np.dot(np.cross(self.rvec[0,:],self.rvec[1,:]),np.array([0,0,1]))) #area of BZ for a 2D lattice - 
    else:

      self.vol= np.abs(np.dot(np.cross(self.rvec[0,:],self.rvec[1,:]),self.rvec[2,:])) #volume of BZ for 3D lattice 


    #calculates the matrix elements on the reducible k-mesh using Sympy package 
    energy_array,diagonal_array,full_array,berry_array=Matrix_Elements_Gen(Hamk,b1,b2,b3,self.nkx,self.nky,self.nkz,self.write_to_file)    

    #load in the energies and matrix elements to hk_class 
    self.energies = energy_array
    self.opticalDiag = diagonal_array
    self.opticalMoments = full_array
    self.berry_curv=berry_array

  def _setup_redkmesh(self): #sets up the reducible k-mesh 
      
      kgrid = np.array([self.nkx,self.nky,self.nkz], dtype=int)
      if self.kshift:
        is_shift = np.array([int(i) for i in self.dims], dtype=np.float64)
      else:
        is_shift = np.array([0,0,0], dtype=int)

      
      self.nkp = self.nkx * self.nky * self.nkz
      self.kpoints = []
      for ikx in np.linspace(0,1,self.nkx,endpoint=False):
        for iky in np.linspace(0,1,self.nky,endpoint=False):
          for ikz in np.linspace(0,1,self.nkz,endpoint=False):
            self.kpoints.append([ikx,iky,ikz])
      self.kpoints = np.array(self.kpoints, dtype=np.float64)
      if self.kshift:
        self.kpoints += is_shift.astype(np.float64)/2. / kgrid[None,:].astype(np.float64)

      self.multiplicity = np.ones((self.nkp), dtype=int)
      self.weights      = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))
      

  def computeData(self, hk_file, charge, mu=None, mushift=False, corronly=False, vector=False):
    
    self.hk_file      = hk_file #file which defines the Hamiltonian in k-space 
    self.charge       = charge #filling 
    self.corronly     = corronly
    self.vector       = vector


    self._read_hk_file()

    self.energyBandMax = self.n_bands 
    self.spins = 1

    if self.charge < 0 or self.charge > self.energyBandMax*2:
      raise ValueError('Provided charge does not match provided bands : charge in [0,2*bands]')

    self._setup_redkmesh()

    ''' setting some important variables '''
    self.opticalBandMin = 0
    self.opticalBandMax = self.energyBandMax

    self._calcFermiLevel(mu)

    if mushift:
      logger.info('Shifting energies -> Chemical Potential === 0.')
      self.energies[0][...] -= self.mu
      # dont forget about possible valence and conduction energies ... np.nan + value = np.nan
      self.ecb[0] -= self.mu
      self.evb[0] -= self.mu
      self.mu = 0.0

    

