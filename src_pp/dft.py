#! /usr/bin/env python

import os
import sys
import numpy as np

import ase.io
import ase.spacegroup
import ase.geometry
# for ase documentation visit
# https://wiki.fysik.dtu.dk/ase/ase/ase.html

from pathlib import Path
import xml.etree.ElementTree as et


# rather bare-bone Exception classes
class MissingFileError(IOError):
  def __init__(self, message):
    super().__init__(message)
    self.message = message

class EmptyFileError(IOError):
  def __init__(self, message):
    super().__init__(message)
    self.message = message

# TODO: energy cutoff for extracted bands
# TODO: energy cutoff for extracted optical elements
# TODO: detect orthogonal lattices
# and truncate moment to length 3 xx yy zz

class DFTcalculation(object):
  def __init__(self, directory, optic):
    self.directory      = directory
    self.optic          = optic
    self.nkp            = 0
    self.charge         = 0
    self.weights        = None
    self.weightsum      = 0
    self.energies       = []
    self.energyBandMax  = 0
    self.opticalBandMin = 0
    self.opticalBandMax = 0
    self.spins          = 1

  def readData(self):
    raise NotImplementedError("DFTcalculation is a non-initializable base class")

  def _checkDirectory(self):
    if Path(self.directory).is_file():
      raise IOError("Supplied directory is not a file: " + self.directory)
    if not Path(self.directory).is_dir():
      raise IOError("Supplied Directory does not exist: " + self.directory)


class w2kcalculation(DFTcalculation):
  def __init__(self, directory, optic=False, calctype=0):
    print("\nInitializing Wien2k calculation.")
    super(w2kcalculation, self).__init__(directory, optic)
    self.calctype      = calctype
    self.moments       = []

    self._checkDirectory()
    # 0 -> detect; 1 -> spin-unpolarized; 2 -> spin-polarized; 3 -> spin-orbit; 4 -> so+sp
    self._defineCase()                      # self.case
    self._defineFiles()

    try:
      if (self.calctype == 0): # we detect it
        self._detectCalculation()
      self._checkFiles()
    except MissingFileError as e:
      print(e.message)
      print("Exiting...")
      sys.exit()
    except EmptyFileError as e:
      print(e.message)
      print("Exiting...")
      sys.exit()
    else:
      print("Files sucessfully loaded.")

    self._getAuxiliaryInformation()         # self.vol; self.atms; self.latVec; self.latAngle

  # internal routine to get the Path prefix case
  def _defineCase(self):
    ''' Get the Path prefix for a Wien2k calculation.
    case = /path/to/folder/calc/calc
    the files then can be described with case.fileending '''

    self.case = self.directory.split('/')
    self.case = list(filter(None,self.case))[-1]

    if self.directory[-1] == '/':
      self.case = self.directory + self.case
    else:
      self.case = self.directory + '/' + self.case

  # define the files which we might use
  def _defineFiles(self):
    self.fscf        = Path(self.case + '.scf')
    self.fstruct     = Path(self.case + '.struct')
    self.fklist      = Path(self.case + '.klist')

    self.fenergysoup = Path(self.case + '.energysoup')
    self.fenergysodn = Path(self.case + '.energysodn')
    self.fenergyso   = Path(self.case + '.energyso')
    self.fenergyup   = Path(self.case + '.energyup')
    self.fenergydn   = Path(self.case + '.energydn')
    self.fenergy     = Path(self.case + '.energy')

    self.fmomentup   = Path(self.case + '.symmatup')
    self.fmomentdn   = Path(self.case + '.symmatdn')
    self.fmoment     = Path(self.case + '.symmat')

  def _checkFiles(self):
    ''' Check the existance of files.
    Define the energy and symmat files we will access when reading later on.
    Define the number of different spins which we have to handle
    Define the weightsum (1 spin -> 2; 2 spins -> 1) for the kpointlist '''

    # mandatory files
    if not self.fstruct.is_file():
      raise MissingFileError("Error: case.struct file missing")
    if not self.fklist.is_file():
      raise MissingFileError("Error: case.klist file missing")

    if self.calctype == 1: # unpolarized
      self.weightsum = 2.
      self.spins = 1
      self.fenergyaccess = [self.fenergy]
      self.fmomentaccess = [self.fmoment]
    elif self.calctype == 2: # spin-polarized
      self.weightsum = 1.
      self.spins = 2
      self.fenergyaccess = [self.fenergyup, self.fenergydn]
      self.fmomentaccess = [self.fmomentup, self.fmomentdn]
    elif self.calctype == 3: # spin-orbit
      self.weightsum = 2.
      self.spins = 1
      self.fenergyaccess = [self.fenergyso]
      self.fmomentaccess = [self.fmoment]
    elif self.calctype == 4: # spin-polarized + spin-orbit
      self.weightsum = 1.
      self.spins = 2
      self.fenergyaccess = [self.fenergysoup, self.fenergysodn]
      self.fmomentaccess = [self.fmomentup, self.fmomentdn]

    # check energy files
    for i in self.fenergyaccess:
      if not i.is_file():
        raise MissingFileError("Error: case.energy* files missing")
      if os.stat(i).st_size == 0:
        raise EmptyFileError("Error: case.energy* files empty")

    # check symmat files
    if self.optic:
      for i in self.fmomentaccess:
        if not i.is_file():
          raise MissingFileError("Error case.symmat* files missing")
        if os.stat(i).st_size == 0:
          raise EmptyFileError("Error: case.symmat* files empty")


  def _detectCalculation(self):
    ''' Deduce from the existing energy files the type of the wien2k calculation.
    set the calculation type for the following methods.'''

    if self.fenergysoup.is_file() and os.stat(self.fenergysoup).st_size != 0 \
       and self.fenergysodn.is_file() and os.stat(self.fenergysodn).st_size != 0:
      self.calctype = 4
      print("Detected spin-polarized calculation with spin-orbit coupling.")
    elif self.fenergyso.is_file() and os.stat(self.fenergyso).st_size != 0:
      self.calctype = 3
      print("Detected spin-unpolarized calculation with spin-orbit coupling.")
    elif self.fenergyup.is_file() and os.stat(self.fenergyup).st_size != 0 \
       and self.fenergydn.is_file() and os.stat(self.fenergydn).st_size != 0:
      self.calctype = 2
      print("Detected spin-polarized calculation without spin-orbit coupling.")
    elif self.fenergy.is_file() and os.stat(self.fenergy).st_size != 0:
      self.calctype = 1
      print("Detected spin-unpolarized calculation without spin-orbit coupling.")
    else:
      raise MissingFileError("Error: No matching energy file combinations found")

  def _getAuxiliaryInformation(self):
    ''' construct an ase Atoms object and get the information
    which is relevant to us '''

    asestruct = ase.io.read(self.fstruct)
    self.vol = asestruct.get_volume() # volume of unit cell
    self.atms = asestruct.get_number_of_atoms() # total number of atoms
    laa = asestruct.get_cell_lengths_and_angles()
    self.latLength = np.array(laa[:3])
    self.latAngle = np.array(laa[3:])
    self.spacegroup = ase.spacegroup.get_spacegroup(asestruct)

  def readData(self):
    self.readScf()
    self.readStruct()
    self.readKlist()
    self.readEnergies()
    if (self.optic):
      self.readMoments()
    print("Files successfully read.")

  def readScf(self):
    print("Reading: ", self.fscf)
    with open(self.fscf, 'r') as scf:
      for line in scf:
        if(':CHA ' in line):
          self.charge = float(line.split()[-1])
          break
        # if(':VOL' in line):
        #   self.vol = float(line.split()[-1])
        #   print(self.vol)

  def readStruct(self):
    ''' Read the number of inequivalent atoms.
    Necessary for the header of the energy files '''

    print("Reading: ", self.fstruct)
    with open(self.fstruct, 'r') as struct:
      self.iatms = 0 # inequivalent atoms
      for line in struct:
        if 'ATOM  ' in line:
          self.iatms = self.iatms + 1

  def readKlist(self):
    ''' Read the kpoint list.
    Read the multiplicity.
    Define the weights. '''

    print("Reading: ", self.fklist)
    with open(self.fklist, 'r') as klist:
      # first we get the divisor
      firstline = klist.readline()
      pos1 = firstline.find('(')+1
      pos2 = firstline.find(')')
      divisor = np.array(firstline[pos1:pos2].split()).astype(int)

      # now we reset the file
      klist.seek(0)

      myklist = []
      mymult  = []
      for line in klist:
        if ("END" in line):
          break
        else:
          myklist.append(line.split()[1:4])
          mymult.append(line.split()[5])

    self.kpoints = np.array(myklist).astype(np.float64)
    self.multiplicity = np.array(mymult).astype(np.float64)
    self.weights = self.multiplicity / sum(self.multiplicity) * self.weightsum
    self.nkp = self.kpoints.shape[0]

    for i in range(self.nkp):
      self.kpoints[i,:] /= divisor

  def readEnergies(self):
    ''' Read all the energies from the files we defined.
    Truncate energies to the least number of existent bands per k-point.
    Rescale read energies to ElectronVolt. '''

    for i in self.fenergyaccess:
      print("Reading: ", i)
      with open(i, 'r') as elist:
        myenergy = []
        bands    = []
        for _ in range(2*self.iatms):
          elist.readline()
        for ikp in range(self.nkp):
          bands.append(int(elist.readline().split()[5])) # number of bands from line above energy data
          myenergy.append([])
          for iband in range(bands[ikp]):
            myenergy[ikp].append(float(elist.readline().split()[1]))

      bands = np.array(bands)
      self.energyBandMax = bands.min() # take only the common bands

      energies = np.empty((self.nkp, self.energyBandMax), dtype=np.float64)

      # put the data in one common array
      # and remove the non-common bands
      for ikp in range(self.nkp):
        energies[ikp,:] = np.array(myenergy[ikp][:self.energyBandMax])

      self.energies.append(energies*13.6056980659) # Rydberg 2 eV

  def readMoments(self):
    ''' Read all the optical elements from the files we defined.
    There have to be 3, 6 or 9 entries.
    If there are elements above the maximum number of energy bands, truncate them.
    '''

    for i in self.fmomentaccess:
      print("Reading: ", i)
      with open(i, 'r') as mlist:

        # number of optical entries
        entries = int(mlist.readline().split()[0]) # read the very first line
        bentries = np.zeros((entries,), dtype=np.float64)
        if (entries % 3 != 0):
          raise ValueError("Number of entries in case.symmat* must be divisable by 3")

        # we require Re[xx] Re[yy] Re[zz]
        # or         Re[xx] Re[yy] Re[zz] Re[xy] Re[xz] Re[yz]
        # or         Re[xx] Re[yy] Re[zz] Re[xy] Re[xz] Re[yz] Im[xy] Im[xz] Im[yz]

        mlist.readline() # skip the next line

        mymoments = []
        bandranges = []
        # moments = np.zeros((self.nkp,bmax,bmax,entries), dtype=np.float64)

        for ikp in range(self.nkp):
          temp = mlist.readline().split() # read the KP: line
          # local (per k-point) band limits
          bmin = int(temp[5])
          bmax = int(temp[6])

          bandranges.append([bmin,bmax])
          mymoments.append(np.zeros((bmax,bmax,entries), dtype=np.float64))

          brange = bmax-bmin+1
          ndatapoints = (brange*(brange+1))//2

          mlist.readline() # skip the next empty line

          for _ in range(ndatapoints): # start reading in the datarange
            line = mlist.readline()
            b1 = int(line[:7])   - 1 # index starts at 0
            b2 = int(line[7:11]) - 1 # same here
            for ientry in range(entries):
              bentries[ientry] = float(line[11+13*ientry:24+13*ientry]) # because we have possible minus signs
            mymoments[ikp][b1,b2,:] = bentries

          # if we are not at the last k-points skip 3 lines
          if (ikp+1 != self.nkp):
            mlist.readline()

        # we are done reading the current file
        # Cleaning up:

        # we get the common energy band interval
        bmin, bmax = bandranges[0][0], bandranges[0][1]
        for ikp in range(1,self.nkp):
          localbmin, localbmax = bandranges[ikp][0], bandranges[ikp][1]
          if localbmax < bmax:
            bmax = localbmax
          if localbmin > bmin:
            bmin = localbmin

        # if in this common interval the bands are for some reason
        # above the band maximum of the energies, truncate it further
        if bmax > self.energyBandMax:
          bmax = self.energyBandMax

        # save this interval as variables
        self.opticalBandMin = bmin
        self.opticalBandMax = bmax

        binterval = bmax-bmin+1
        # this has always 6 entries and is complex
        moments = np.zeros((self.nkp,binterval,binterval,6), dtype=np.complex128)

        for ikp in range(self.nkp):
          if (entries == 3): # no off-diagonal elements
            moments[ikp,:,:,:3] = mymoments[ikp][bmin-1:bmax,bmin-1:bmax,:]
          elif (entries == 6): # real off-diagonal elements
            moments[ikp,:,:,:] = mymoments[ikp][bmin-1:bmax,bmin-1:bmax,:]  # real off-diagonal elements
          elif (entries == 9): # complex off-diagonal elements
            moments[ikp,:,:,:]   = mymoments[ikp][bmin-1:bmax,bmin-1:bmax,:6]
            moments[ikp,:,:,3:6] += 1j * mymoments[ikp][bmin-1:bmax,bmin-1:bmax,6:9]


        # now we symmetrize the bands
        # because we only got the upper triangle

        # in order to get the 'lower half' we use the relation
        # mopt(m,n)^ij = mopt(n,m)^ji = [mopt(n,m)^ij]*
        # which means:
        # the lower band triangle = upper band triangle conjugated
        for ikp in range(self.nkp):
          for i in range(6):
            view = moments[ikp,:,:,i] # this is a memory view
            # we need this in order to do this fancy accessing
            il1 = np.tril_indices(binterval,-1)
            iu1 = np.triu_indices(binterval,1)
            view[il1] = view[iu1].conjugate() # for i = 0,1,2 the conjugate doesnt do anything

        # a list of spin-dependent arrays
        self.moments.append(moments)

  def getData(self):
    pass


class vaspcalculation(DFTcalculation):
  def __init__(self, directory):
    print("\nInitializing VASP calculation.")
    super(vaspcalculation, self).__init__(directory, optic=False)

    self._checkDirectory() # from base class
    self._defineCase()     # file prefix
    self._defineFiles()
    self._getAuxiliaryInformation()

    try:
      self._checkFiles()
    except MissingFileError as e:
      print(e.message)
      print("Exiting...")
      sys.exit()
    except EmptyFileError as e:
      print(e.message)
      print("Exiting...")
      sys.exit()
    else:
      print("Files sucessfully loaded.")

  # internal routine to get the Path prefix case
  def _defineCase(self):
    if self.directory[-1] == '/':
      self.case = self.directory
    else:
      self.case = self.directory + '/'

  # we only care about the vasprun file
  def _defineFiles(self):
    self.fvasprun    = Path(self.case + 'vasprun.xml')

  def _checkFiles(self):
    if not self.fvasprun.is_file():
      raise MissingFileError("Error: vasprun.xml file missing.")
    if os.stat(self.fvasprun).st_size == 0:
      raise EmptyFileError("Error: vasprun.xml is empty.")

  def _getAuxiliaryInformation(self):
    asestruct = ase.io.read(self.fvasprun)
    self.vol = asestruct.get_volume() # volume of unit cell
    self.atms = asestruct.get_number_of_atoms() # total number of atoms
    laa = asestruct.get_cell_lengths_and_angles()
    self.latLength = np.array(laa[:3])
    self.latAngle = np.array(laa[3:])
    self.spacegroup = ase.spacegroup.get_spacegroup(asestruct)

  def readData(self):
    print("Reading: ", self.fvasprun)
    self.xmlroot = et.parse(self.fvasprun).getroot()
    self._read_ispin()
    self._read_nelect()
    self._read_vol()
    self._read_kpointlist()
    self._read_weights()
    self._read_energies()
    print("Files successfully read.")

  def _read_ispin(self):
    self.spins = int(self.xmlroot.find('incar/i[@name="ISPIN"]').text)
    if (self.spins == 1):
      self.weightsum = 2
    elif (self.spins == 2):
      self.weightsum = 1

  def _read_nelect(self):
    self.charge = float(self.xmlroot.find('parameters/separator[@name="electronic"]/i[@name="NELECT"]').text)

  def _read_vol(self):
    self.vol = float(self.xmlroot.find('structure[@name="primitive_cell"]/crystal/i[@name="volume"]').text)

  def _read_kpointlist(self):
    arr = self.xmlroot.find('kpoints/varray[@name="kpointlist"]')
    kpointlist = []
    for ikp in arr:
      kpointlist.append(ikp.text.split())

    self.nkp = len(kpointlist)
    self.kpoints = np.array(kpointlist, dtype=np.float64)

  def _read_weights(self):
    arr = self.xmlroot.find('kpoints/varray[@name="weights"]')
    weightlist = []
    for ikp in arr:
      weightlist.append(ikp.text)

    self.weights = np.array(weightlist, dtype=np.float64) * self.weightsum

  def _read_energies(self):
    for ispin in range(self.spins):
      energylist = []
      spinarr = self.xmlroot.find('calculation/eigenvalues/array/set/set[@comment="spin {}"]'.format(ispin+1))
      for ikp in range(self.nkp):
        energylist.append([])
        kpointarr = spinarr.find('set[@comment="kpoint {}"]'.format(ikp+1))
        for ene in kpointarr:
          energylist[ikp].append(ene.text.split()[0])
      self.energies.append(np.array(energylist, dtype=np.float64))


if __name__ == '__main__':

  calc = w2kcalculation('../aFe-sp', optic=True)
  calc.readData()

  # print("\n\n\n")
  # calc = vaspcalculation('../test_vasp/')
  # calc.readData()

  print(calc.kpoints)
  print(calc.weights)
  print(np.sum(calc.weights))
  print(calc.energies)
  print(calc.moments[1][28])
  print(calc.energyBandMax)
  print(calc.opticalBandMin)
  print(calc.opticalBandMax)
