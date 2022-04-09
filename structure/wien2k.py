#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import logging
logger = logging.getLogger(__name__)

from structure.aux import progressBar

import numpy as np

try:
  import ase
  import ase.io
  import ase.spacegroup
  from ase import Atoms
  from ase.spacegroup import get_spacegroup
  mymethods = {'get_volume', 'get_global_number_of_atoms'}
  Atommethods = set([method for method in dir(Atoms) if callable(getattr(Atoms, method))])
  ase_exists = True if mymethods.issubset(Atommethods) else False
except ImportError:
  ase_exists = False

from structure.dft import DftCalculation
from structure     import units

class Wien2kCalculation(DftCalculation):
  '''
  Wien2k calculation class which can load all possible Wien2K calculations
  (unpolarized, unpolarized + LS, polarized, polarized + LS)
  and extract the necessary information for a transport calculation.
  We require the files: case.scf, case.klist, case.struct, case.energy*, case.symmat*
  The exact energy and symmat files depend on the type of calculation.
  For this reason we can either detect the calculation (calctype=0)
  or provide the type of calculation:
    0: detect it automatically
    1: unpolarized:
    2: polarized:
    3: unpolarized + SOC + inversion symmetry:
    4: unpolarized + SOC + no inversion symmetry:
    5: polarized + SOC:
  The optical elements (stored in case.symmat*) can either be ignored (optic=False)
  or be loaded (optic=True)
  The file beginning (case) is extracted automatically (case=None)
  or can be provided specifically (e.g. case="abc")
  '''

  def __init__(self, directory, optic=False, calctype=0, case=None, **kwargs):
    logger.info("\nInitializing Wien2k calculation.")
    super(Wien2kCalculation, self).__init__(directory)
    self.optic         = optic   # exising full optical elements
    self.calctype      = calctype
    self.case          = case

    if self.optic: # otherwise already set by parent class
      self.opticdiag = True
      self.opticfull = True

    self._checkDirectory()
    if self.case is None:
      self._defineCase(custom=False)
    else:
      self._defineCase(custom=True)
    self._defineFiles()

    if (self.calctype == 0): # we detect it
      self._detectCalculation()
    self._checkFiles()
    logger.info("Files sucessfully loaded.")

    if ase_exists:
      ''' get all the crystal information from the ase atom object '''
      self._getAuxiliaryInformation()         # self.vol; self.atms; self.latVec; self.latAngle
      logger.info('ASE: unit cell volume: {} [Angstrom^3]'.format(self.vol))

  def __repr__(self):
    return ('w2kcalculation(directory={0.directory}, '
            'optic={0.optic!r}, '
            'calctype={0.calctype!r}, '
            'case={0.case!r})'.format(self))

  # internal routine to get the Path prefix case
  def _defineCase(self, custom=False):
    '''
    Get the Path prefix for a Wien2k calculation.
    case = /path/to/folder/calc/calc
    the files then can be described with case.fileending
    Wien2K usually enforces folder/folder.filending
    However sometimes this is not the case -> custom fileneding
    -> path/custom.filending
    '''

    if not custom:
      self.case = os.path.basename(self.directory) # get the case prefix
    self.case = os.path.join(self.directory,self.case)

  # define the files which we might use
  def _defineFiles(self):
    self.fscf        = self.case + '.scf'
    self.fstruct     = self.case + '.struct'
    self.fklist      = self.case + '.klist'

    self.fenergysoup = self.case + '.energysoup'
    self.fenergysodn = self.case + '.energysodn'
    self.fenergyso   = self.case + '.energyso'
    self.fenergyup   = self.case + '.energyup'
    self.fenergydn   = self.case + '.energydn'
    self.fenergy     = self.case + '.energy'

    self.fmomentup   = self.case + '.symmatup'
    self.fmomentdn   = self.case + '.symmatdn'
    self.fmoment     = self.case + '.symmat'

    if ase_exists:
      self.aseobject = ase.io.read(self.fstruct)

  def _checkFiles(self):
    '''
    Check the existance of files.
    Define the energy and symmat files we will access when reading later on.
    Define the number of different spins which we have to handle
    Define the weightsum (1 spin -> 2; 2 spins -> 1) for the kpointlist
    '''

    # mandatory files
    if not os.path.isfile(self.fstruct):
      raise IOError("Error: case.struct file missing")
    if not os.path.isfile(self.fklist):
      raise IOError("Error: case.klist file missing")

    if self.calctype == 1: # unpolarized
      self.spinorbit = False
      self.weightsum = 2
      self.spins = 1
      self.fenergyaccess = [self.fenergy]
      self.fmomentaccess = [self.fmoment]
    elif self.calctype == 2: # spin-polarized
      self.spinorbit = False
      self.weightsum = 1
      self.spins = 2
      self.fenergyaccess = [self.fenergyup, self.fenergydn]
      self.fmomentaccess = [self.fmomentup, self.fmomentdn]
    elif self.calctype == 3: # spin-orbit with inversion
      self.spinorbit = True
      self.weightsum = 1    # this is 1 because w2k outputs essentially double the number of bands per k-point
                             # in this kind of calculation
      self.spins = 2
      self.fenergyaccess = [self.fenergyso]
      self.fmomentaccess = [self.fmoment]
    elif self.calctype == 4: # spin-orbit without inversion
      self.spinorbit = True
      self.weightsum = 1
      self.spins = 2
      self.fenergyaccess = [self.fenergyso]
      self.fmomentaccess = [self.fmomentup]
    elif self.calctype == 5: # spin-polarized + spin-orbit
      self.spinorbit = True
      self.weightsum = 1
      self.spins = 2
      self.fenergyaccess = [self.fenergysoup]
      self.fmomentaccess = [self.fmomentup]

    # check energy files
    for i in self.fenergyaccess:
      if not os.path.isfile(i):
        raise IOError("Error: case.energy* files missing")
      if os.stat(str(i)).st_size == 0:
        raise IOError("Error: case.energy* files empty")

    # check symmat files
    if self.optic:
      for i in self.fmomentaccess:
        if not os.path.isfile(i):
          raise IOError("Error: case.symmat* files missing")
        if os.stat(str(i)).st_size == 0:
          raise IOError("Error: case.symmat* files empty")

    logger.info("  Number of spins: {}".format(self.spins))


  def _detectCalculation(self):
    '''
    Deduce from the existing energy files the type of the wien2k calculation.
    set the calculation type for the following methods.
    '''

    if os.path.isfile(self.fenergysoup) and os.stat(str(self.fenergysoup)).st_size != 0:
       # and os.path.isfile(self.fenergysodn) and os.stat(str(self.fenergysodn)).st_size != 0:
      self.calctype = 5
      logger.info("Detected spin-polarized calculation with spin-orbit coupling.")
    elif os.path.isfile(self.fenergyso) and os.stat(str(self.fenergyso)).st_size != 0:
      self.calctype = 3
      logger.info("Detected spin-unpolarized calculation with spin-orbit coupling.")
      # further differentiate the type of symmat file we have available
      if self.optic:
        if os.path.isfile(self.fmomentup) and os.stat(str(self.fmomentup)).st_size != 0:
          self.calctype = 4
    elif os.path.isfile(self.fenergyup) and os.stat(str(self.fenergyup)).st_size != 0 \
       and os.path.isfile(self.fenergydn) and os.stat(str(self.fenergydn)).st_size != 0:
      self.calctype = 2
      logger.info("Detected spin-polarized calculation without spin-orbit coupling.")
    elif os.path.isfile(self.fenergy) and os.stat(str(self.fenergy)).st_size != 0:
      self.calctype = 1
      logger.info("Detected spin-unpolarized calculation without spin-orbit coupling.")
    else:
      raise IOError("Error: No matching energy file combinations found")


  def readData(self):
    self._readScf()
    self._readStruct()
    self._readKlist()
    self._readEnergies()
    if (self.optic):
      self._readMoments()
    logger.info("Files successfully read.")
    self._calcFermiLevel()

  def _readScf(self):

    logger.info("Reading: {}".format(self.fscf))
    with open(str(self.fscf), 'r') as scf:
      for line in scf:
        if line.startswith(':LABEL3'):
          try:
            self.version = line[15:26] # this stays a string
          except:
            pass
          break
    logger.info('  Wien2K version: {}'.format(self.version))

    with open(str(self.fscf), 'r') as scf:
      for line in scf:
        if line.startswith(':NOE '): # number of electrons
          self.charge = float(line[38:])
          break
      else:
        raise IOError('Wien2k {}: Did not find Number of electrons (:NOE)'.format(str(self.fscf)))
    logger.info('  Number of electrons: {}'.format(self.charge))

    if not ase_exists:
      with open(str(self.fscf), 'r') as scf:
        for line in scf:
          if line.startswith(':VOL '): # volume of the unit cell
            # self.vol = float(line.split()[-1])
            self.vol = float(line[26:])
            self.vol *= units.bohr2angstrom**3 # from bohr**3 to angstrom**3
            break
        else:
          raise IOError('Wien2k {}: Did not find Crystal Cell Volume (:VOL)'.format(str(self.fscf)))
      logger.info('  Unit cell volume: {} [Angstrom^3]'.format(self.vol))

  def _readStruct(self):
    '''
    Read the number of inequivalent atoms.
    Necessary for the header of the energy files.
        Get the orthogonality of our system (check for 90 degree angles).
    '''

    logger.info("Reading: {}".format(self.fstruct))
    with open(str(self.fstruct), 'r') as struct:
      # get the number of inequivalent atoms
      self.iatms = 0 # inequivalent atoms
      for line in struct:
        if line.startswith('ATOM'):
          self.iatms += 1

      if self.iatms == 0:
        raise IOError('Wien2k {}: Reading number of inequivalent atoms failed.'.format(str(self.fstruct)))
      else:
        logger.info("  Number of inequivalent atoms: {}".format(self.iatms))

      struct.seek(0)
      # get the number of symmetry operations
      for line in struct:
        if 'NUMBER OF SYMMETRY OPERATIONS' in line:
          self.nsym = int(line[:4])
          logger.info("  Number of symmetry operations: {}".format(self.nsym))
          break
      else:
        raise IOError('Wien2k {}: Reading number of symmetry operations failed.'.format(str(self.fstruct)))

      # without resetting we can continue reading the operations
      self.symop    = np.zeros((self.nsym,3,3), dtype=np.float64)
      self.invsymop = np.zeros((self.nsym,3,3), dtype=np.float64)
      for isym in range(self.nsym):
        for i in range(3):
          temp = struct.readline()
          self.symop[isym,i,:] = np.array([temp[:2],temp[2:4],temp[4:6]], dtype=np.float64)
        struct.readline() # skip
        self.invsymop[isym] = np.linalg.inv(self.symop[isym])

      if struct.readline() != "": # exactly at the EOF
        raise IOError('Wien2K {} is not at the EOF after reading'.format(str(self.fstruct)))

      logger.debug('Symmetry operations')
      logger.debug(self.symop)


  def _readKlist(self):
    '''
    Read the kpoint list.
    Read the multiplicity.
    Define the weights.
    '''

    logger.info("Reading: {}".format(self.fklist))
    with open(str(self.fklist), 'r') as klist:
      # first we get the divisor
      firstline = klist.readline()
      pos1 = firstline.find('(')+1
      pos2 = firstline.find(')')
      divisor = []
      divisor.append(int(firstline[pos1:pos1+3]))
      divisor.append(int(firstline[pos1+3:pos1+6]))
      divisor.append(int(firstline[pos1+6:pos1+9]))
      divisor = np.array(divisor, dtype=np.int)

      self.nkx, self.nky, self.nkz = divisor
      # determine dimension whether we find '1's in the divisor
      self.dims = np.logical_not(divisor == np.ones(3, dtype=np.int))
      self.ndim = 3 - np.sum(divisor == np.ones(3, dtype=np.int))


      # now we reset the file
      klist.seek(0)

      myklist = []
      mymult  = []
      for line in klist:
        if ("END" in line):
          break
        else:
          kx   = float(line[10:20])
          ky   = float(line[20:30])
          kz   = float(line[30:40])
          mult = int(float(line[50:55]))
          myklist.append([kx,ky,kz])
          mymult.append(mult)
      else: # if we did not find the break condition
        raise IOError('Wien2k {}: Did not find END statement.'.format(str(self.fklist)))

    self.kpoints = np.array(myklist, dtype=np.float64)
    self.multiplicity = np.array(mymult, dtype=int)
    self.weights = self.weightsum * self.multiplicity / float(np.sum(self.multiplicity))
    self.nkp = self.kpoints.shape[0]
    self.irreducible = not (self.nkx*self.nky*self.nkz == self.nkp)

    for i in range(self.nkp):
      self.kpoints[i,:] /= divisor


    logger.info("  Number of dimensions: {}".format(self.ndim))
    logger.info("  Number of k-points: {}".format(self.nkp))

  def _readEnergies(self):
    '''
    Read all the energies from the files we defined.
    Truncate energies to the least number of existent bands per k-point.
    Rescale read energies to ElectronVolt.
    '''

    for i in self.fenergyaccess:
      logger.info("Reading: {}".format(i))
      with open(str(i), 'r') as elist:
        try:
          myenergy = []
          bands    = []
          for _ in range(2*self.iatms):
            elist.readline()
          for ikp in range(self.nkp):
            # bands.append(int(elist.readline()[57:].split()[2])) # number of bands from line above energy data
            bands.append(int(elist.readline()[73:79])) # number of bands from line above energy data
            myenergy.append([])
            for iband in range(bands[ikp]):
              # myenergy[ikp].append(float(elist.readline().split()[1]))
              myenergy[ikp].append(float(elist.readline()[12:]))
        except Exception as e:
          raise IOError('{}\nWien2K energy file {} contains an illegal data format'.format(str(e),str(i)))

        if elist.readline() != "": # we exactly have the EOF
          raise IOError('Wien2K energy file {} is not at the EOF after reading'.format(str(i)))

      bands = np.array(bands)
      self.energyBandMax = bands.min() # take only the common bands

      energies = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)

      # put the data in one common array
      # and remove the non-common bands
      for ikp in range(self.nkp):
        energies[ikp,:] = np.array(myenergy[ikp][:self.energyBandMax])

      # for the spin unpolarized + spin-orbit case we have to extract the spins from one file ...
      # in the energyso file the energies are listed with up dn up dn up dn ...
      if self.spinorbit:
        self.energies.append(energies[:,::2]*units.rydberg2eV) # spin up
        self.energies.append(energies[:,1::2]*units.rydberg2eV) # spin up
        self.energyBandMax = int(self.energyBandMax / 2.) # we definitely have 2n bands from wien2k
      else:
        self.energies.append(energies*units.rydberg2eV)

    # check for spin dependent maximum number of bands
    # if there is a difference truncate the bigger one
    if self.spins == 2:
      bmax1, bmax2 = self.energies[0].shape[-1], self.energies[1].shape[-1]
      if bmax1 != bmax2:
        bmax = min(bmax1,bmax2)
        self.energies[0] = self.energies[0][:,:,:bmax]
        self.energies[1] = self.energies[1][:,:,:bmax]
        self.energyBandMax = bmax

    logger.info("  Number of common bands: {}".format(self.energyBandMax))

  def _readMoments(self):
    '''
    Read all the optical elements from the files we defined.
    There have to be 3, 6 or 9 entries.
    If there are elements above the maximum number of energy bands, truncate them.
    '''

    for i in self.fmomentaccess:
      logger.info("Reading: {}".format(i))
      with open(str(i), 'r') as mlist:
        try:
          # number of optical entries
          line = mlist.readline()
          entries = int(line[10:11])
          bentries = np.zeros((entries,), dtype=np.float64)
          if (entries % 3 != 0):
            raise CustomError("Wien2k symmat file {}: Number of entries must be divisable by 3".format(str(i)))

          symm = np.ones((entries,), dtype=np.float64)
          if entries == 9:
            symm[6:] = -1.0

          # we require Re[xx] Re[yy] Re[zz] (ortho, with and without spin-orbit)
          # or         Re[xx] Re[yy] Re[zz] Re[xy] Re[xz] Re[yz] (non-ortho, no spin-orbit)
          # or         Re[xx] Re[yy] Re[zz] Re[xy] Re[xz] Re[yz] Im[xy] Im[xz] Im[yz] (non-ortho, spin-orbit)

          mlist.readline() # skip the next line

          mymoments = []
          bandranges = []
          # moments = np.zeros((self.nkp,bmax,bmax,entries), dtype=np.float64)

          for ikp in range(self.nkp):
            if logger.getEffectiveLevel() in [logging.DEBUG,logging.INFO]:
              progressBar(ikp+1,self.nkp, status='k-points')
            # temp = mlist.readline().split() # read the KP: line
            temp = mlist.readline()

            # local (per k-point) band limits
            bmin = int(temp[26:32]) - 1
            bmax = int(temp[32:37])
            # bmin and bmax are now properly iterable in python

            # bandranges contain 1-indexed values
            bandranges.append([bmin,bmax])

            # we define the 'load-array':
            # here we ignore the lower bound
            mymoments.append(np.zeros((bmax,bmax,entries), dtype=np.float64))

            mlist.readline() # skip the next empty line

            for b1 in range(bmin, bmax):
              for b2 in range(b1,bmax):

                line = mlist.readline()
                for ientry in range(entries):
                  bentries[ientry] = float(line[11+13*ientry:24+13*ientry]) # because we have possible minus signs
                mymoments[ikp][b1,b2,:] = bentries

                if b1 != b2:
                  mymoments[ikp][b2,b1,:] = bentries*symm

            # if we are not at the last k-points skip 3 lines
            if (ikp+1 != self.nkp):
              mlist.readline()
        except CustomError as e:
          raise e
        except Exception as e: # if any error occurs
          raise IOError('{}\nWien2k symmat file {} contains an illegal data format'.format(str(e),str(i)))

        if mlist.readline() != "": # we exactly have the EOF
          raise IOError('Wien2K symmat file {} is not at the EOF after reading'.format(str(i)))

        # we are done reading the current file
        # Cleaning up:

        # we get the common energy band interval
        # bmin bmax are 1-indexed
        bmin, bmax = bandranges[0][0], bandranges[0][1]
        for ikp in range(1,self.nkp):
          localbmin, localbmax = bandranges[ikp][0], bandranges[ikp][1]
          if localbmax < bmax:
            bmax = localbmax
          if localbmin > bmin:
            bmin = localbmin

        # if in this common interval the bands are for some reason
        # above the band maximum of the energies, truncate it further

        # again : differentiate between the spin-unpolarized + spin-orbit and the rest
        if self.spinorbit:
          if bmax > 2*self.energyBandMax:
            bmax = 2*self.energyBandMax
        else:
          if bmax > self.energyBandMax:
            bmax = self.energyBandMax

        # save this interval as variables
        self.opticalBandMin = bmin
        self.opticalBandMax = bmax

        binterval = bmax-bmin
        moments = np.zeros((self.nkp,binterval,binterval,entries), dtype=np.float64)

        for ikp in range(self.nkp):
          moments[ikp,:,:,:] = mymoments[ikp][bmin:bmax,bmin:bmax,:]

        # unit conversion to eV**2 * angstrom**2
        moments *= units.w2kmom

        if self.spinorbit:
          self.opticalMoments.append(moments[:,::2,::2,:])
          self.opticalMoments.append(moments[:,1::2,1::2,:])
          self.opticalBandMax = self.opticalBandMin + int((self.opticalBandMax - self.opticalBandMin) / 2.)
        else:
          self.opticalMoments.append(moments)


        nelements = self.opticalMoments[0].shape[1]
        if self.spinorbit:
          self.opticalDiag.append(self.opticalMoments[0][:,np.arange(nelements),np.arange(nelements),:])
          self.opticalDiag.append(self.opticalMoments[1][:,np.arange(nelements),np.arange(nelements),:])
        else:
          self.opticalDiag.append(self.opticalMoments[0][:,np.arange(nelements),np.arange(nelements),:])

    logger.info("  Optical band range: {} - {}".format(self.opticalBandMin+1,self.opticalBandMax))
