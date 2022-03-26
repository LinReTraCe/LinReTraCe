#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import argparse
import warnings
from collections import OrderedDict

import numpy as np
with warnings.catch_warnings():
  warnings.filterwarnings("ignore",category=FutureWarning)
  import h5py

from structure import es
from structure.dos import calcDOS

class LRTCoutput(object):
  '''
  Output class for the main output of LRTC
  We initialize the file with either a path to a hdf5 file
  or with "latest". Latest will iterate through all available hdf5 files
  in the current folder and check for the LRTC "identifier"

  The two main output routines are "saveData" and "outputData"
  saveData purely saves the requested data in the data dictionary
  while outputData additionally outputs said data either to stdout or plots it via matplotlib
  '''
  def __init__(self, fname, altaxis=False, dp=None):
    self.fname    = fname.strip()
    self.datasets = OrderedDict({})
    self.owned    = OrderedDict({})
    self.data     = None
    self.dataspinsum = None

    self._parse()        # runmode, quad, dimensions
    if dp:
      self.quad = not dp # overwrite quad if wanted
      print('#\n#   DEBUG: USING DOUBLE PRECISION DATA')

    self._get_axis(altaxis) # get T / beta / mu / carrier axis, set the wanted one to self.axis
    self._defineDicts()     # define all possible response datasets internally
    self._retrieve_groups() # define available response datasets + spin

    if sys.version_info >= (3, 0): #  for numpy
      self.textpipe = sys.stdout.buffer
    else:
      self.textpipe = sys.stdout
    sys.stdout.flush()      # so all the output afterwards comes after ther information

  def __repr__(self):
    return ('LRTCoutput(fname={0.fname!r})'.format(self))

  def __getitem__(self,key):
    if self.data is not None:
      return self.data[key]
    else:
      return None

  def __len__(self):
    if self.data is not None:
      return len(self.data)
    else:
      return 0

  def __iter__(self):
    if self.data is not None:
      return iter(self.data)
    else:
      return None

  def _defineDicts(self):
    '''
    Define all possible available dataset.
    Raw datasets contain a path in the hdf5 file
    while derived datasets contain a list of required raw datasets
    '''

    # for raw quantities
    # key: userinput; value: (raw dset, (internal path(s)), description, response, magnetic)
    #                        (True    ,  ...         ,  ...       , True/False, True/False)
    # for derived quantities
    # key: userinput; value: (raw dset, requirements,  description, response, magnetic)
    #                        (False   ,  ...         ,  ...       , True      , True/False)

    # key : what the user inputs
    # raw dset : direct datasets or onsager coefficients // conductivity - resistivity - etc are ' derived '
    # internal path : hdf5 path
    # requirements : refers to keys of raw quantities
    # response : false for direct datasets, true for onsager and derived thereof
    # magnetic : magnetic onsager or derived quantities that require magnetic onsager (hall, nernst, mobilities)


    # 'raw' direct quantities
    # dos is only listed here so the 'list' command shows is, we use a different method to calculate it

    self.datasets.update({'dos':        (True, '.energies',                    'Density of States',                          False, False)})
    self.datasets.update({'energy':     (True, '.quantities/energy',           'Total of energy of the system [eV]',         False, False)})
    self.datasets.update({'mu':         (True, '.quantities/mu',               'Chemical potential [eV]',                    False, False)})
    self.datasets.update({'occupation': (True, '.quantities/occupation',       'Total occupation in the system',             False, False)})
    self.datasets.update({'carrier':    (True, '.quantities/carrier',          'Carrier concentration in the system w.r.t. neutral charge [cm^-3]', False, False)})
    self.datasets.update({'electrons':  (True, '.quantities/electrons',        'Thermally activated electrons',              False, False)})
    self.datasets.update({'holes':      (True, '.quantities/holes',            'Thermally activated holes',                  False, False)})
    self.datasets.update({'imp':        (True, '.quantities/imp_contribution', 'Thermally activated impurity electrons',     False, False)})


    quadprefix = 'Quad' if self.quad else ''

    # 'raw' Onsager coefficients
    for iL, unit in zip(['L11','L12','L22'],['V/(A*m)','A/m', 'V*A/m']):
      for ii in ['inter','intra']:
        for iM, iMdescr, iMflag in zip(['','M'],['',' in magnetic field'],[False,True]):
          for iB, iBdescr in zip(['','Boltz'],['','Boltzmann']):
            key = '{}{}-{}{}'.format(iL,iM,ii,iB) # L11M-intraBoltz
            if key[-1] == '-': key = key[:-1]
            internalpath = '{}{}/{}{}{}/sum'.format(iL,iM,ii,quadprefix,iBdescr.strip())
            description  = '{} {} {}{} [{}{}]'.format(iL,ii,iBdescr,iMdescr, unit, ' (m*m)/(V*s)' if iM else '') # m^2/(Vs) = 1/T
            self.datasets.update({key : (True, internalpath, description, True, iMflag)})

    # 'derived' quantities ... that are constructed from the Onsager coefficients
    for iL, iLreq, iLdescr, unit, magnetic in zip(['r','c','p','s','tc','tr','rh','n','muh','mut'], \
            [('L11',),('L11',),('L11','L12'),('L11','L12'),('L11','L12','L22'),('L11','L12','L22'),('L11M','L11'),('L11M','L12M','L11','L12'),('L11','L11M'),('L12','L12M')], \
            ['resistivity', 'conductivity','peltier', 'seebeck', 'thermal conductivity', 'thermal resistivity', 'hall', 'nernst', 'hall mobility', 'thermal mobility'], \
            ['[Ohm*m]','[1/(Ohm*m)]','[V]','[V/K]','[W/(m*K)]','[m*K/W]', '[m^3/C]', '[V/(K*T)]', '[1/T]', '[1/T]'], \
            [False,False,False,False,False,False,True,True,True,True]):
      for ii, iireq in zip(['inter','intra','total'], [('inter',), ('intra',), ('inter','intra')]):
        for iB, iBdescr in zip(['','Boltz'],['','Boltzmann']):

          key = '{}-{}{}'.format(iL,ii,iB)
          requirement = []
          for i in iLreq:
            for j in iireq:
              requirement.append(i+'-'+j+iB)

          description = '{} {} {} {}'.format(iLdescr,ii,iBdescr,unit)
          self.datasets.update({key : (False, requirement, description, True, magnetic)})


  def saveData(self, command, diag=False, *args):
    '''
    Check if the provided command is valid for the given file.
    Save the collected data in self.data in form of a OrderedDict
    '''

    command = command.strip()

    if command in self.owned:
      derived  = not self.owned[command][0]
      response = self.owned[command][3]
      magnetic = self.owned[command][4]
      if derived:
        commands = self.owned[command][1]
      else:
        commands = [command]
    else:
      raise IOError('Provided dataset does not exist.')

    if self.data is None:
      # we define the dictionary and the temperatures
      # the first time we save something
      self.data        = OrderedDict({})
      self.dataspinsum = OrderedDict({})
      self.data.update({'temp':self.temp})
      self.data.update({'invtemp':self.invtemp})
      self.data.update({'mu':self.mu})
      self.data.update({'carrier':self.carrier})

    if response:
      for icmd in commands: # iterate through all the required items
        key = self.owned[icmd][1] # path
        # get the spin-resolved onsager coefficients
        out = self._getResponseCombination(key, spinsum=False)
        self.data.update({icmd:out})
        # get the spin-summed onsager coefficients
        out = self._getResponseCombination(key, spinsum=True)
        self.dataspinsum.update({icmd:out})
    else:
      if len(args) != 0:
        print('#   Warning: This group does not take additional arguments')
      key = self.owned[command][1] # internal hdf5 key
      out = self._getQuantity(key)
      self.data.update({command:out})


    # transform the data from onsager coefficients into physical observables
    if response and derived: # combine the saved data
      requirements = sorted(self.owned[command][1]) # so we have L11 L11M L12 L12M L22 L22M
      # print(requirements)
      # this sorting is vital for the array indexing below

      combined = []
      combinedspinsum = []
      for ireq in requirements:
        combined.append(self.data[ireq])
        combinedspinsum.append(self.dataspinsum[ireq])

      if command.find('total') != -1: # found
        ''' we add inter + intra '''
        total = []
        totalspinsum = []
        ''' we appended the data in a sorted fashion '''
        for i in range(len(combined)//2):
          total.append(combined[2*i]+combined[2*i+1])
          totalspinsum.append(combinedspinsum[2*i]+combinedspinsum[2*i+1])
      else:
        total        = combined
        totalspinsum = combinedspinsum

      del combined
      del combinedspinsum

      for i, itotal in enumerate([total,totalspinsum]):

        if i==0:
          if magnetic:
            temp = self.temp[:,None,None,None,None] # steps, spins, dir, dir, dir
          else:
            temp = self.temp[:,None,None,None] # steps, spins, dir, dir
        else:
          if magnetic:
            temp = self.temp[:,None,None,None] # steps, dir, dir, dir
          else:
            temp = self.temp[:,None,None] # steps, dir, dir

        if command.startswith('c-'): # conductivity
          tosave = itotal[0]
        elif command.startswith('r-'): # resistivity
          tosave = self.invert(itotal[0])
        elif command.startswith('p-'): # peltier
          tosave = -np.einsum('...ij,...jk->...ik', self.invert(itotal[0]), itotal[1])
        elif command.startswith('s-'): # seebeck
          tosave = -np.einsum('...ij,...jk->...ik', self.invert(itotal[0]), itotal[1]) / temp
        elif command.startswith('tc-'): # thermal conductivity
          tosave = itotal[2] - np.einsum('...ij,...jk,...kl->...il', itotal[1], self.invert(itotal[0]), itotal[1])
          tosave /= temp
        elif command.startswith('tr-'): # thermal resistivity
          tosave = itotal[2] - np.einsum('...ij,...jk,...kl->...il', itotal[1], self.invert(itotal[0]), itotal[1])
          tosave /= temp
          tosave = self.invert(tosave)
        elif command.startswith('rh-'): # Hall coefficient
          tosave = np.einsum('...ij,...jkz,...kl->...ilz', self.invert(itotal[0]), itotal[1], self.invert(itotal[0]))
        elif command.startswith('n-'): # Nernst coefficient
          tosave  = np.einsum('...ij,...jkz,...kl,...lm->...imz', self.invert(itotal[0]), itotal[1], itotal[2], self.invert(itotal[0]))
          tosave -= np.einsum('...ij,...jkz,...kl,...lm->...imz', self.invert(itotal[0]), itotal[3], itotal[0], self.invert(itotal[0]))
          tosave /= temp
          tosave *= (-1.)
        elif command.startswith('muh-'): # Hall mobility
          tosave  = np.einsum('...ij,...jkz->...ikz', self.invert(itotal[0]), itotal[1])
        elif command.startswith('mut-'): # Thermal mobility
          tosave  = np.einsum('...ij,...jkz->...ikz', self.invert(itotal[0]), itotal[1])
        else:
          raise IOError('Cannot recognize command')

        if i==0:
          self.data.update({command:tosave})
        else:
          self.dataspinsum.update({command:tosave})


  def plotBandgap(self):
    '''
    When plotting the chemical potential:
    Plot the maximum of the valence and the minimum of the conduction band
    Plot the impurity states.
    '''

    import matplotlib.pyplot as plt

    with h5py.File(self.fname,'r') as h5:
      enev = []
      enec = []
      fullgap = True
      for ispin in range(self.spins):
        if self.spins == 1:
          prefix = '/'
        else:
          if ispin == 0:
            prefix = '/up'
          else:
            prefix = '/dn'

        gapped = h5['.quantities/bandgap'+prefix+'/gapped'][()]

        if gapped:
          enev.append(h5['.quantities/bandgap'+prefix+'/ene_vband'][()])
          enec.append(h5['.quantities/bandgap'+prefix+'/ene_cband'][()])
        else:
          fullgap = False
          enev.append(np.nan)
          enec.append(np.nan)

      if fullgap:
        plt.axhline(y=np.min(enev), color='black', lw=2)
        plt.axhline(y=np.max(enec), color='black', lw=2)

      nimp = h5['.quantities/impurities/nimp'][()]
      if nimp > 0:
        for iimp in range(nimp):
          eimp = h5['.quantities/impurities/imp-{:03}/energy'.format(iimp+1)][()]
          dop  = h5['.quantities/impurities/imp-{:03}/dopant'.format(iimp+1)][()]
          wid  = h5['.quantities/impurities/imp-{:03}/width'.format(iimp+1)][()]
          if wid < 1e-7:
            plt.axhline(y=eimp, color='red' if dop==1. else 'blue', lw=2)
          else:
            plt.axhspan(eimp-wid/2., eimp+wid/2., color='red' if dop==1. else 'blue', alpha=0.5)

      if nimp == 1 and fullgap:
        eimp = h5['.quantities/impurities/imp-001/energy'][()]
        dop  = h5['.quantities/impurities/imp-001/dopant'][()]
        wid  = h5['.quantities/impurities/imp-001/width'][()]

        enec = np.min(enec)
        enev = np.min(enev)

        if dop == 1:
          elvl = ( enec + (eimp + wid/2.) ) /2.
        else:
          elvl = ( enev + (eimp - wid/2.) ) /2.
        plt.axhline(y=elvl,  color='gray', lw=1, ls='-.')

      # # some matplotlib stuff to create a secondary axis, fixed to the primary one
      # # used only for labellying
      # # does not work properly
      # ax1 = plt.gca()
      # ax2 = ax1.twinx()
      # ax2.set_ylim(ax1.get_ylim())
      # def fix_secondaxis(ax1):
      #   ax2.set_ylim(ax1.get_ylim())
      #   ax2.figure.canvas.draw()
      # ax1.callbacks.connect("ylim_changed", fix_secondaxis)
      # plt.sca(ax1) # return focus



  def outputData(self, command, imag=False, plot=False, diag=False, *args):
    '''
    User interface for lprint.
    Save the data via saveData
    Output the collected data to stdout
    or plot it with matplotlib.
    '''

    if plot:
      import matplotlib.pyplot as plt

    self.saveData(command, diag, *args)
    self.headerwritten = False

    response = self.owned[command][3]
    magnetic = self.owned[command][4]

    '''
    check arguments in more detail
    raise valuerrors if inconsistencies are detected
    '''
    if len(args) > 0:
      if magnetic:
        for icomb in args:
          if len(icomb) == 3 and str(icomb)[0] == '0': # gets added in the main file automatically
            raise ValueError("Invalid directional argument: -- incorrect argument length [use e.g. xyz, uxxz'")
          if len(icomb) < 3 or len(icomb) > 4:
            raise ValueError("Invalid directional argument: -- incorrect argument length [use e.g. xyz, uxxz'")
      else:
        for icomb in args:
          if len(icomb) < 2 or len(icomb) > 3:
            raise ValueError("Invalid directional argument: -- incorrect argument length [use e.g. xx, uxy]")

    if response:
      outfull = self.data[command]
      outspinsum = self.dataspinsum[command]

      for ispin, ispindescr in zip(range(-1,2), ['', 'up ', 'dn ']): # -1 refers to the spin-summed quantity
        for idir1, idir1descr in zip(range(3), ['x','y','z']):
          for idir2, idir2descr in zip(range(3), ['x','y','z']):
            for idir3, idir3descr in zip(range(3), ['x','y','z']):

              if len(outfull.shape) == 4: # steps, spins, dir1, dir2
                idir3 = None

              if len(args)==0: # iterate through all possible combinations

                if ispin == -1: # skip spin-summed elements if we plot all combinations
                  continue

                if diag: # skip non-diagonal elements
                  if idir3 is None:
                    if idir1!=idir2: continue
                  else:
                    if idir1==idir2 or idir2==idir3 or idir1==idir3: continue

                if idir3 is None:
                  icomb =  str(ispin+1)+str(idir1+1)+str(idir2+1)
                  icombdescr = ispindescr+idir1descr+idir2descr
                else:
                  icomb =  str(ispin+1)+str(idir1+1)+str(idir2+1)+str(idir3+1)
                  icombdescr = ispindescr+idir1descr+idir2descr+idir3descr

              else: # check the input combinations
                if idir3 is None:
                  icomb =  str(ispin+1)+str(idir1+1)+str(idir2+1)
                  icombdescr = ispindescr+idir1descr+idir2descr
                else:
                  icomb =  str(ispin+1)+str(idir1+1)+str(idir2+1)+str(idir3+1)
                  icombdescr = ispindescr+idir1descr+idir2descr+idir3descr

                if icomb not in args:
                  continue

              if ispin >= 0:
                if idir3 is None:
                  outarray = outfull[:,ispin,idir1,idir2]
                else:
                  outarray = outfull[:,ispin,idir1,idir2,idir3]
              else:
                if idir3 is None:
                  outarray = outspinsum[:,idir1,idir2]
                else:
                  outarray = outspinsum[:,idir1,idir2,idir3]


              if plot:
                plt.plot(self.axis, outarray.real, label='{}.real [{}]'.format(command, icombdescr))
                if imag: plt.plot(self.axis, outarray.imag, label='{}.imag [{}]'.format(command, icombdescr))
              else:
                if idir3 is None:
                  auxarray = np.zeros((self.nT,3), dtype=np.int)
                  auxarray[None,:] = np.array([ispin+1,idir1+1,idir2+1], dtype=np.int)

                  if not self.headerwritten:
                    np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None], auxarray)), \
                               fmt='%25.15e %30.18e %30.18e %5i %2i %2i', \
                               header='  {0}{1}, {2:>31}.real, {2:>24}.imag,           is id1 id2'.format \
                               (self.axisname,self.axisunit,command))
                    self.headerwritten = True

                  else:
                    np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None], auxarray)), \
                               fmt='%25.15e %30.18e %30.18e %5i %2i %2i', comments='', header='\n')
                else:
                  auxarray = np.zeros((self.nT,4), dtype=np.int)
                  auxarray[None,:] = np.array([ispin+1,idir1+1,idir2+1,idir3+1], dtype=np.int)

                  if not self.headerwritten:
                    np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None], auxarray)), \
                               fmt='%25.15e %30.18e %30.18e %5i %2i %2i %2i', \
                               header='  {0}{1}, {2:>31}.real, {2:>24}.imag,           is id1 id2 id3'.format \
                               (self.axisname,self.axisunit,command))
                    self.headerwritten = True

                  else:
                    np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None], auxarray)), \
                               fmt='%25.15e %30.18e %30.18e %5i %2i %2i %2i', comments='', header='\n')

              # we have plotted it now, now break the idir3 loop
              # if this is not done we do it twice more
              if idir3 is None:
                break
    else:
      outarray = self.data[command]
      if plot:
        plt.plot(self.axis, outarray, label='{}'.format(command))
      else:
        np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray[:,None])), header='{}{}, {}'.format \
        (self.axisname,self.axisunit, self.owned[command][1]))


  def outputList(self, full=False):
    '''
    List the internally existing data sets.
    full=False does not list the raw responses (L0 ...)
    full=True  lists all datasets
    '''

    barlength = 55

    print('\n{:<12}  {}'.format('Key', 'Description [unit]'))
    print(barlength*u'\u2500')

    # quantities
    for (key, value) in self.owned.items():
      raw_dset, path, description, response, magnetic = value
      if not response:
        print('{:<12}  {}'.format(key, description))
    print(barlength*u'\u2500')

    if full:
      # raw responses
      for (key, value) in self.owned.items():
        raw_dset, path, description, response, magnetic = value
        if response and raw_dset:
          print('{:<18}  {}'.format(key, description))
      print(barlength*u'\u2500')
    else:
      # derived responses
      for (key, value) in self.owned.items():
        raw_dset, requirements, description, response, magnetic = value
        if response and not raw_dset:
          print('{:<18}  {}'.format(key, description))
      print(barlength*u'\u2500')

  def outputConfig(self):
    '''
    List all the saved config parameters
    '''

    with h5py.File(self.fname, 'r') as h5:
      self.config = h5['.config'].attrs

      barlength = 55

      print('\n{:<22}  {}'.format('Parameter', 'Value'))
      print(barlength*u'\u2500')

      for i in self.config:
        try:
          print('{:<22}  {}'.format(i, self.config[i].decode('utf-8'))) # byte string to string
        except:
          print('{:<22}  {}'.format(i, self.config[i]))

      print(barlength*u'\u2500')


      if self.spins == 1:
        try:
          gapped = h5['.quantities/bandgap/gapped'][()]
          gap    = h5['.quantities/bandgap/gapsize'][()]
          print('gap [eV]:',gap)
        except:
          print('no gap')
      else:
        try:
          gappedup = h5['.quantities/bandgap/up/gapped'][()]
          gapup    = h5['.quantities/bandgap/up/gapsize'][()]
          print('up: gap [eV]:',gapup)
        except:
          print('up: no gap')
        try:
          gappeddn = h5['.quantities/bandgap/up/gapped'][()]
          gapdn    = h5['.quantities/bandgap/dn/gapsize'][()]
          print('dn: gap [eV]:',gapdn)
        except:
          print('dn: no gap')

      if (self.config['doping']):
        print('doping: ', h5['.quantities/doping'][()])

      print(barlength*u'\u2500')

      print('{}-mode steps: {}'.format(self.mode, self.temp.shape[0]))
      if self.mode == 'temp':
        print('tmin [K]: {}\ntmax [K]: {}'.format(self.temp[0],self.temp[-1]))
      elif self.mode == 'mu':
        print('mumin [eV]: {}\nmumax [eV]: {}'.format(self.mu[0],self.mu[-1]))

      print(barlength*u'\u2500')

      for i in h5['.scattering'].keys():
        print('{:<22}  {}'.format(i, h5['.scattering'][i][()]))

      print(barlength*u'\u2500')

      try:
        if self.config['impurities']:
          print('----- : dopant, density, energy, degeneracy, width')
          nimp = h5['/.quantities/impurities/nimp'][()]
          for i in range(1,nimp+1):
            deg = h5['/.quantities/impurities/imp-001/degeneracy'][()]
            dens= h5['/.quantities/impurities/imp-001/density'][()]
            dop = h5['/.quantities/impurities/imp-001/dopant'][()]
            ene = h5['/.quantities/impurities/imp-001/energy'][()]
            wid = h5['/.quantities/impurities/imp-001/width'][()]
            print('imp {} : {} {} {} {} {}'.format(i,dop,dens,ene,deg,wid))
          print(barlength*u'\u2500')
      except:
        pass


  def _getResponseCombination(self, key, spinsum):
    '''
    Get the key numpy array from the file
    '''

    with h5py.File(self.fname,'r') as h5:
      outputarray = h5['{}'.format(key)][()]

    if spinsum:
      data = np.sum(outputarray, axis=1)
      return data

    else:
      # artificially introduce spins
      if self.spins == 1:
        shape = list(outputarray.shape)
        shape[1] = 2
        spinshape = np.zeros(shape, dtype=np.complex128)
        spinshape[...] = outputarray / 2.
        return spinshape

      else:
        return outputarray

  def _getQuantity(self, key):
    '''
    Retrieve the given key (for a quantities (mu, energy)
    and return the array.
    '''

    outputarray = np.zeros((self.nT,), dtype=np.float64)
    with h5py.File(self.fname,'r') as h5:
      outputarray[:] = h5[format(key)][()]
    return outputarray


  def _retrieve_groups(self):
    '''
    Iterate through all the allowed entries and check
    for existance in our file.
    If the entry exist, add it to a dict.
    '''

    # first we iterate through the raw_datasets
    with h5py.File(self.fname, 'r') as h5:
      for (key, value) in self.datasets.items(): # for compatibility reasons
        raw_dset, path, description, response, magnetic = value

        if not raw_dset:
          continue
        if response:
          exist = '{}'.format(path) in h5
        else:
          exist = path in h5
        if exist:
          self.owned.update({key:value})

    # now we iterate through the derived datasets
    for (key, value) in self.datasets.items():
      raw_dset, requirements, desdcription, response, magnetic = value

      if raw_dset:
        continue
      allcontained = True
      for ireq in requirements:
        if ireq not in self.owned:
          allcontained = False
          break
      if allcontained:
        self.owned.update({key:value})

  def _parse(self):
    '''
    user method to parse the file and
    retrieve the available groups.
    '''

    if self.fname == 'latest':
      self._parse_latest()
    else:
      self._parse_file()


  def _parse_latest(self):
    '''
    Get all the hdf5 files sorted by timestamp.
    In reverse order: check if we have one of 'our' output files
    we can parse.
    '''

    directoryfiles = [f for f in os.listdir('.') if os.path.isfile(f)]
    filteredfiles  = filter(lambda x: x.endswith('.hdf5'), directoryfiles)
    sortedfiles    = sorted(filteredfiles, key=lambda x: os.path.getmtime(x), reverse=True)

    success = False
    for fi in sortedfiles:
      try:
        with h5py.File(fi,'r') as hfi:
          if hfi['.quantities'].attrs['identifier'].decode("utf-8") == 'LRTC':
            # problems here will be caught by the outer try - except
            self.fname = fi
            break
          else:
            continue
      except:
        continue
    else:
      raise IOError('No LRTC output files detected in current folder.')

    self._parse_file()


  def _parse_file(self):
    '''
    Check if the provided file is one of 'our' output files.
    Detect the 'run-mode' of the calculation
    Detect if we have quad precision response
    Detect the number of dimensions and which dimensions are valid
    (necessary for quantities that require an inversion)
    '''

    try:
      with h5py.File(self.fname,'r') as hfi:
        # we get a byte string here
        # hence we have to decode it to utf-8
        if hfi['.quantities'].attrs['identifier'].decode("utf-8") != 'LRTC':
          raise IOError('Provided file is not an LRTC output file.')
        if hfi['.quantities'].attrs['mode'].decode("utf-8") == 'temp':
          self.mode = 'temp'
        elif hfi['.quantities'].attrs['mode'].decode("utf-8") == 'mu':
          self.mode = 'mu'
        print('#   Using file:', self.fname)
        print('#   Detected run mode:', self.mode)

        self.quad = hfi['.config'].attrs['quad']
        if self.quad:
          print('#   Detected Quad Precision')
        self.spins = hfi['.quantities/ispin'][()]
        self.ndim = hfi['.unitcell/ndim'][()]
        self.dims = hfi['.unitcell/dims'][()]
        self.dimmask2 = np.logical_and(self.dims[:,None], self.dims[None,:])
        self.dimmask3 = np.logical_and(np.logical_and(self.dims[:,None], self.dims[None,:])[:,:,None], self.dims[None,None,:])
        print('#   Detected {} dimensions: {}'.format(self.ndim, np.array(["x","y","z"])[self.dims]))
    except:
      raise IOError('Provided file is not an LRTC output file.')

  def _get_axis(self, altaxis):
    '''
    Get the temperature and inverse temperature axis
    Also save the number of temperature steps
    '''
    with h5py.File(self.fname, 'r') as h5:
      self.temp    = h5['.quantities/tempAxis'][()]
      self.invtemp = h5['.quantities/betaAxis'][()]
      self.carrier = h5['.quantities/carrier'][()]
      self.mu      = h5['.quantities/mu'][()]
      self.mudft   = h5['.quantities/mudft'][()]
      self.nT      = self.temp.shape[0]

    if self.mode == 'temp':
      if altaxis:
        self.axis = self.invtemp
        self.axisname = 'beta'
        self.axisunit = '[eV^{-1}]'
        self.axislatex = r'$\beta$ [eV$^{-1}$]'
      else:
        self.axis = self.temp
        self.axisname = 'T'
        self.axisunit = '[K]'
        self.axislatex = r'$T$ [K]'
    elif self.mode == 'mu':
      if altaxis:
        self.axis = self.carrier
        self.axisname = 'n'
        self.axisunit = '[cm^{-3}]'
        self.axislatex = r'$n$ [cm$^{-3}$]'
      else:
        self.axis = self.mu
        self.axisname = 'mu'
        self.axisunit = '[eV]'
        self.axislatex = r'$\mu$ [eV]'
    else:
      raise ValueError('no properly defined x-axis to plot')

  def outputDOS(self, plot, broadening=0.02):
    '''
    Print DOS/NOS
    '''

    import matplotlib.pyplot as plt

    with h5py.File(self.fname,'r') as h5:
      mu  = h5['.quantities/mudft'][()]
      spins = h5['.quantities/ispin'][()]
      weights = h5['.quantities/weights'][()]

      if spins==1:
        ene = h5['.energies'][()]
        dosaxis, dos, nos = calcDOS(ene, weights, gamma=broadening, windowsize=1.1)
        del ene
      else:
        eneup = h5['.energies/up'][()]
        dosaxisup, dosup, nosup = calcDOS(eneup, weights, gamma=broadening, windowsize=1.1)
        del eneup
        enedn = h5['.energies/dn'][()]
        dosaxisdn, dosdn, nosdn = calcDOS(enedn, weights, gamma=broadening, windowsize=1.1)
        del enedn

      if plot:
        if spins==1:
          plt.plot(dosaxis, dos, color='black', lw=2, label='DOS')
          plt.legend(loc='upper left')
          plt.ylabel(r'DOS [eV$^{-1}$]')
          plt.xlabel(r'$\varepsilon$ [eV]')

          plt.twinx()
          plt.ylabel('NOS')
          plt.plot(dosaxis, nos, color='gray', lw=2, label='NOS')
          plt.legend(loc='upper right')
        else:
          plt.plot(dosaxisup, dosup, color='blue', lw=2, label='DOS up')
          plt.plot(dosaxisdn, -dosdn, color='red', lw=2, label='DOS dn')
          ylim = max(np.max(dosup),np.max(dosdn)) * 1.1
          plt.ylim(-ylim,ylim)
          plt.legend(loc='upper left')
          plt.ylabel(r'DOS [eV$^{-1}$]')
          plt.xlabel(r'$\varepsilon$ [eV]')

          plt.twinx()
          plt.plot(dosaxisup, nosup, color='deepskyblue', lw=2, label='NOS up')
          plt.plot(dosaxisup, -nosdn, color='indigo', lw=2, label='NOS dn')
          ylim = max(np.max(nosup),np.max(nosdn)) * 1.1
          plt.ylim(-ylim,ylim)
          plt.legend(loc='upper right')
          plt.ylabel('NOS')
        plt.axvline(x=mu, ls='--', color='gray', lw=2)
        plt.xlabel('energy')
      else:
        if spins==1:
          np.savetxt(self.textpipe, np.hstack((dosaxis[:,None], dos[:,None], nos[:,None])), \
                     fmt='%25.15e %25.15e %25.15e', comments='', \
                     header='# energy [eV], DOS [eV^-1], NOS]')
        else:
          np.savetxt(self.textpipe, np.hstack((dosaxis[:,None], dosup[:,None], dosdn[:,None], nosup[:,None], nosdn[:,None])), \
                     fmt='%25.15e %25.15e %25.15e %25.15e %25.15e', comments='', \
                     header='#  energy [eV], DOSup [eV^-1], DOSdn [eV^-1], NOSup, NOSdn')

    if plot:
      plt.title(self.fname)
      plt.show()
    print('') # empty line before next CLI input

  def invert(self, data):
    '''
    Invert given data
    Provided data has either shape of [nT,spins,3,3]
    or [nT,3,3] ... need to do checks to differentiate

    Given full dimensionality: straight forward inversion and output
    Reduced dimensionality: reduce data to the reduced dimensions, perform inversion
    and 'blow' up to full 3x3 afterwards

    in this sense the resistivity on an invalid axis is === 0
    necessary to avoid inversion attemps on singular matrices.
    '''

    spins = True if len(data.shape)==4 else False

    if self.ndim < 3 and self.ndim > 0: # treat 0D identical to 3D
      # select direction combinations
      if spins:
        outputmasked = data[:,:,self.dimmask2].copy()
        outputmasked = outputmasked.reshape((data.shape[0],data.shape[1],self.ndim,self.ndim))
      else:
        outputmasked = data[:,self.dimmask2].copy()
        outputmasked = outputmasked.reshape((data.shape[0],self.ndim,self.ndim))
    elif self.ndim == 3:
      outputmasked = data
    elif self.ndim == 0:
      outputmasked = data
      # there cannot be mixed directions in 0D per definition
      outputmasked[...,[0,0,1,1,2,2],[1,2,0,2,0,1]] = 0.0


    # always applied on the last two elements
    # inverted = np.linalg.inv(outputmasked)

    inverted = np.full_like(outputmasked, fill_value=np.nan, dtype=np.complex128)
    if spins:
      for i in range(outputmasked.shape[0]):
        for j in range(outputmasked.shape[1]):
          try:
            inverted[i,j,:,:] = np.linalg.inv(outputmasked[i,j,:,:])
          except:
            pass
    else:
      for i in range(outputmasked.shape[0]):
        try:
          inverted[i,:,:] = np.linalg.inv(outputmasked[i,:,:])
        except:
          pass


    if self.ndim < 3 and self.ndim > 0:
      returned = np.zeros_like(data, dtype=np.complex128)
      ii = -1
      for i in range(3):
        if self.dims[i]:
          ii += 1
        else:
          continue

        jj = -1
        for j in range(3):
          if self.dims[j]:
            jj += 1
          else:
            continue

          if spins:
            returned[:,:,i,j] = inverted[:,:,ii,jj]
          else:
            returned[:,i,j] = inverted[:,ii,jj]
    else:
      returned = inverted

    return returned
