#! /usr/bin/env python

import numpy as np
import dft
import sys

try:
  import BoltzTraP2.dft as BTP
  import BoltzTraP2.bandlib as BL
  import BoltzTraP2.io as IO
  from BoltzTraP2 import sphere
  from BoltzTraP2 import fite
  from BoltzTraP2 import serialization
  from BoltzTraP2.misc import ffloat
  from BoltzTraP2 import units
except ImportError:
  print("BoltzTrap2 Library could not be loaded...")
  print("Exiting...")
  sys.exit()


class MetaW2kLoader(BTP.GenericWien2kLoader):
  ''' BoltzTrap Custom Wien2k Loader.
  After setting the class variables one can register the Loader
  und use the provided custom files.
  The usual Wien2kLoader can only access the energy and energyso files.
  We also want to access up and dn files in spin-polarized calculations. '''

  # define class variables
  weightsum = None
  fscf      = None
  fstruct   = None
  fenergy   = None

  # access them here
  def __init__(self, directory):
      BTP.GenericWien2kLoader.__init__(self, "SystemName", \
                                       MetaW2kLoader.weightsum, \
                                       MetaW2kLoader.fscf, \
                                       MetaW2kLoader.fstruct, \
                                       MetaW2kLoader.fenergy)

  @classmethod
  def setfiles(cls, weightsum, fscf, fstruct, fenergy):
    cls.weightsum = weightsum
    cls.fscf      = fscf
    cls.fstruct   = fstruct
    cls.fenergy   = fenergy



class BTPInterpolation(object):
  def __init__(self, dftcalc, niter=30):
    self.dftcalc = dftcalc
    self.niter   = niter # interpolation parameter

    self.energies   = []
    self.velocities = []
    self.curvatures = []

  def interpolate(self):
    if isinstance(self.dftcalc, dft.w2kcalculation):
      for ispin in range(self.dftcalc.spins):

        # set class energy files
        MetaW2kLoader.setfiles(self.dftcalc.weightsum, self.dftcalc.fscf, \
                               self.dftcalc.fstruct, self.dftcalc.fenergyaccess[ispin])
        # and register the new loader
        BTP.register_loader("Custom", MetaW2kLoader)
        # interpolate and save
        self._interp()
        self._save1()
    elif isinstance(self.dftcalc, dft.vaspcalculation):
      BTP.register_loader("VASP", BTP.VASPLoader)
      self._interp()
      if (self.dftcalc.spins == 1):
        self._save1()
      else:
        self._save2() # we have to split the interpolated arrays we get from BTP


  def _interp(self):
    ''' Standard BoltzTrap2 Library interface to interpolate
    the band energies, band velocities and band curvatures '''

    self.data = BTP.DFTData(self.dftcalc.directory, derivatives=False)
    self.equivalences = sphere.get_equivalences(self.data.atoms, self.data.magmom, \
                                                self.niter * len(self.data.kpoints))
    self.coeffs = fite.fitde3D(self.data, self.equivalences)
    self.metadata = serialization.gen_bt2_metadata(self.data, self.data.magmom is not None)

    self.lattvec = self.data.get_lattvec()
    self.interp_energies, self.interp_velocities, self.interp_curvatures = \
        fite.getBands(self.data.kpoints, self.equivalences, self.lattvec, self.coeffs, curvature=True)

  def _save1(self):
    ''' Saving the energies, velocities and curvatures by simply
    appending them to the lists we initialized at init time '''

    # we get the energies on the Hartree scale
    # rescaling to eV!
    self.energies.append(self.interp_energies.transpose(1,0) / units.eV)
    self.velocities.append(self.interp_velocities.transpose(2,1,0) / units.eV)
    self.curvatures.append(self.interp_curvatures.transpose(3,2,0,1) / units.eV)

  def _save2(self):
    ''' In the case of a spin-dependent VASP calculation BTP2
    creates arrays which lists the data in order energyup energydn.
    We want them to be split, which is what this routines does '''

    nbands = self.interp_energies.shape[0] # this is guaranteed to be even here
    self.energies.append(self.interp_energies[:nbands//2,:].transpose(1,0) / units.eV)
    self.energies.append(self.interp_energies[nbands//2:,:].transpose(1,0) / units.eV)
    self.velocities.append(self.interp_velocities[:,:nbands//2,:].transpose(2,1,0) / units.eV)
    self.velocities.append(self.interp_velocities[:,nbands//2:,:].transpose(2,1,0) / units.eV)
    # the last two elements dont matter here, since its symmetric anyways
    self.curvatures.append(self.interp_curvatures[:,:,:nbands//2,:].transpose(3,2,0,1) / units.eV)
    self.curvatures.append(self.interp_curvatures[:,:,nbands//2:,:].transpose(3,2,0,1) / units.eV)

if __name__ == '__main__':
  pass
