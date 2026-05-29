#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import argparse
import warnings
import logging
logger = logging.getLogger(__name__)

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

  Standard-deviation support (disorder averaging workflow)
  ---------------------------------------------------------
  When the HDF5 file was produced by the disorder averaging script
  (workflow_step7_transport_avg.py), each Onsager group contains companion
  datasets ``sum_std`` and ``sum_std_imag`` alongside ``sum``.
  This class detects them via ``self.has_std`` and exposes them through
  ``self.datastd`` / ``self.datastd_spinsum``.

  For raw Onsager coefficients and conductivity (c-) the std is read
  directly.  For derived ratio quantities (Seebeck s-, Hall rh-, Nernst n-,
  etc.) an extremal-envelope approach is used: the formula is evaluated at
  all 2^M corners of the +/-sigma hypercube (M = number of participating
  Onsager inputs) and the pointwise min/max across those corners gives a
  conservative uncertainty band.  See saveData for details.

  Notes on dtypes
  ---------------
  The Fortran code writes Onsager arrays as complex(8); h5py reads them as
  numpy complex128.  L11 is purely real; L12/L22 have small imaginary parts.
  The ``sum_std`` datasets are float64 (real-valued), written by
  np.std(arr.real, ...) in the averaging script.
  '''
  def __init__(self, fname, altaxis=False):
    self.fname    = fname.strip()
    self.datasets = {}
    self.owned    = {}
    self.data        = None
    self.dataspinsum = None
    self.datastd         = None   # std of real parts (float64); populated when has_std is True
    self.datastd_spinsum = None
    self.datastd_lo      = None   # extremal-envelope lower bound for derived quantities
    self.datastd_lo_ss   = None
    self.datastd_hi      = None   # extremal-envelope upper bound for derived quantities
    self.datastd_hi_ss   = None

    self._parse()        # runmode, quad, dimensions

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
    # dos is only listed here so the 'list' command shows it, we use a different method to calculate it

    self.datasets.update({'dos':        (True, '.structure/energies',             'Density of States',                                    False, False)})
    self.datasets.update({'energy':     (True, '.quantities/energy',              'Total of energy of the system [eV]',                   False, False)})
    self.datasets.update({'mu':         (True, '.quantities/mu',                  'Chemical potential [eV]',                              False, False)})
    self.datasets.update({'occupation': (True, '.quantities/occupation',          'Total occupation in the system',                       False, False)})
    self.datasets.update({'carrier':    (True, '.quantities/carrier',             'Carrier concentration w.r.t. neutral charge [cm^-3]',  False, False)})
    self.datasets.update({'electrons':  (True, '.quantities/electrons',           'Thermally activated electrons',                        False, False)})
    self.datasets.update({'holes':      (True, '.quantities/holes',               'Thermally activated holes',                            False, False)})
    self.datasets.update({'impurity':   (True, '.quantities/imp_contribution',    'Thermally activated impurity electrons',               False, False)})
    self.datasets.update({'doping':     (True, '.quantities/doping_contribution', 'Additional system doping',                             False, False)})


    # 'raw' Onsager coefficients
    for iL, unit in zip(['L11','L12','L22'],['V/(A*m)','A/m', 'V*A/m']):
      for ii in ['inter','intra']:
        for iM, iMdescr, iMflag in zip(['','B'],['',' in magnetic field'],[False,True]):
          for iB, iBdescr in zip(['','Boltz'],['','Boltzmann']):
            key = '{}{}-{}{}'.format(iL,iM,ii,iB) # L11B-intraBoltz
            if key[-1] == '-': key = key[:-1]
            internalpath = '{}{}/{}{}/sum'.format(iL,iM,ii,iBdescr.strip()) # L11B/intraBoltzmann/sum

            quantity_description = '{}{} Onsager {}'.format(iL,iM if iM else ' ',iBdescr) #  resistivitiy Boltzmann
            type_description = '({})'.format(ii)
            unitplus = '[' + unit + '{}]'.format(' 1/T' if iM else '')
            description = '{0:<35} {1:<8} {2:>15}'.format(quantity_description, type_description, unitplus)
            # description  = '{} {} {}{} [{}{}]'.format(iL,ii,iBdescr,iMdescr, unit, ' (m*m)/(V*s)' if iM else '') # m^2/(Vs) = 1/T
            self.datasets.update({key : (True, internalpath, description, True, iMflag)})

    # 'derived' quantities ... that are constructed from the Onsager coefficients
    for iL, iLreq, iLdescr, unit, magnetic in zip(['r','c','p','s','pf','tc','tr','cb','rh','n','muh','mut'], \
            [('L11',),('L11',),('L11','L12'),('L11','L12'),('L11','L12'),('L11','L12','L22'),('L11','L12','L22'),('L11B',),('L11B','L11'),('L11B','L12B','L11','L12'),('L11','L11B'),('L12','L12B')], \
            ['Resistivity', 'Conductivity','Peltier coeff', 'Seebeck coeff', 'Power factor', 'Thermal conductivity', 'Thermal resistivity', 'Hall conductivity', 'Hall coeff', 'Nernst coeff', 'Hall mobility', 'Thermal mobility'], \
            ['[Ohm*m]','[1/(Ohm*m)]','[V]','[V/K]','[W/(K^2*m)]','[W/(K*m)]','[K*m/W]', '[A*m^2/(V^2*s)]', '[m^3/C]', '[V/(K*T)]', '[1/T]', '[1/T]'], \
            [False,False,False,False,False,False,False,True,True,True,True,True]):
      for ii, iireq in zip(['inter','intra','total'], [('inter',), ('intra',), ('inter','intra','inter_anti',)]):
        for iB, iBdescr in zip(['','Boltz'],['','Boltzmann']):

          key = '{}-{}{}'.format(iL,ii,iB)
          requirement = []
          for i in iLreq:
            for j in iireq:
              requirement.append(i+'-'+j+iB)

          quantity_description = '{} {}'.format(iLdescr,iBdescr) #  resistivitiy Boltzmann
          type_description = '({})'.format(ii)
          description = '{0:<35} {1:<8} {2:>15}'.format(quantity_description, type_description, unit)
          self.datasets.update({key : (False, requirement, description, True, magnetic)})


  def saveData(self, command, *args):
    '''
    Check if the provided command is valid for the given file.
    Save the collected data in self.data in form of a dictionary.

    Standard-deviation handling (disorder averaging workflow)
    ---------------------------------------------------------
    When self.has_std is True the method also populates parallel dicts
    self.datastd / self.datastd_spinsum for raw Onsager coefficients,
    and self.datastd_lo / self.datastd_hi (and spinsum variants) for
    derived quantities.

    *Raw Onsager coefficients and conductivity (c-)*:
      The std is read directly from sum_std (float64, real part only).
      For the band-summed quantity (c-) which is just L11, this is exact.

    *Derived ratio quantities (Seebeck s-, Hall rh-, Nernst n-, ...)*:
      Exact error propagation would require the full covariance matrix
      between Onsager inputs, which is not stored.  Instead we use an
      extremal-envelope approach:

        For M participating Onsager inputs each with element-wise std
        sigma_k, evaluate the formula f at all 2^M corners of the
        +/-sigma hypercube:

          f_corner = f(L1 +/- sigma1, L2 +/- sigma2, ...)

        The uncertainty band is:
          lower = pointwise min over all corners
          upper = pointwise max over all corners

        This is a rigorous conservative bound when correlations between
        Onsager inputs are unknown, and is exact when all inputs are
        monotone in the formula (which is true for S = L12/L11 etc.).
        It costs 2^M formula evaluations; since M <= 4 (L11, L12, L22
        and one B-field coefficient) this is at most 16 evaluations.

    Limitations:
      * The 2^M bound is conservative: if disorder fluctuations in L11
        and L12 are positively correlated (they usually are), the true
        Seebeck band is narrower than what we show.  The band is an
        outer bound, not a 1-sigma interval.
      * Near a metal-insulator transition where L11 is small the inversion
        is ill-conditioned and even the outer bound may be misleading.
      * Propagated bands for derived quantities are only shown when the
        user explicitly passes --fullstd to lprint.
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
      raise IOError('Provided dataset does not exist. "lprint <file> list" to list all output containers. ')

    if self.data is None:
      # we define the dictionary and the temperatures
      # the first time we save something
      self.data        = {}
      self.dataspinsum = {}
      self.data.update({'temp':self.temp})
      self.data.update({'invtemp':self.invtemp})
      self.data.update({'mu':self.mu})
      self.data.update({'carrier':self.carrier})

    if self.datastd is None and self.has_std:
      self.datastd         = {}
      self.datastd_spinsum = {}
      self.datastd_lo      = {}
      self.datastd_lo_ss   = {}
      self.datastd_hi      = {}
      self.datastd_hi_ss   = {}

    if response:
      for icmd in commands: # iterate through all the required items
        key = self.owned[icmd][1] # path
        # get the spin-resolved onsager coefficients
        out = self._getResponseCombination(key, spinsum=False)
        self.data.update({icmd:out})
        # get the spin-summed onsager coefficients
        out = self._getResponseCombination(key, spinsum=True)
        self.dataspinsum.update({icmd:out})
        # load std counterpart when available
        if self.has_std:
          out_std = self._getResponseCombination(key, spinsum=False, std=True)
          self.datastd.update({icmd: out_std})
          out_std_ss = self._getResponseCombination(key, spinsum=True, std=True)
          self.datastd_spinsum.update({icmd: out_std_ss})
    else:
      if len(args) != 0:
        print('#   Warning: This group does not take additional arguments')
      key = self.owned[command][1] # internal hdf5 key
      out = self._getQuantity(key)
      self.data.update({command:out})

    ''' we massage the description to produce nice looking ylabels and units '''
    if response and not derived:
      description = self.owned[commands[0]][2]
      ylabel = description.split()[0]
      unit = description[description.find("["):]
      unit  = unit[1:-1].replace(' ','\cdot').replace('*',' ')
      unit = r'$[' + unit + r']$'

    # transform the data from onsager coefficients into physical observables
    if response and derived: # combine the saved data
      requirements = sorted(self.owned[command][1]) # so we have L11 L11B L12 L12B L22 L22B
      # print(requirements)
      # this sorting is vital for the array indexing below

      unit  = self.owned[command][2].split()[-1]
      ''' make the unit look nicer '''
      unit  = unit[1:-1].replace('Ohm', r'\Omega').replace('*',' ')
      unit = r'$[' + unit + r']$'

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

        # --- helper: evaluate the transport formula from a list of Onsager arrays ---
        # Factored out so the same code path serves both the nominal evaluation
        # and the extremal-envelope std computation (see saveData docstring).
        def _eval_derived(cmd, L, tmp):
          '''
          Evaluate the transport formula for cmd given Onsager list L and
          temperature broadcast array tmp.
          Returns (result_array, ylabel_string).
          '''
          if cmd.startswith('c-'):
            return L[0], r'$\sigma$'
          elif cmd.startswith('r-'):
            return self.invert(L[0]), r'$\rho$'
          elif cmd.startswith('p-'):
            return -np.einsum('...ij,...jk->...ik', self.invert(L[0]), L[1]), r'$\Pi$'
          elif cmd.startswith('s-'):
            return -np.einsum('...ij,...jk->...ik', self.invert(L[0]), L[1]) / tmp, r'$S$'
          elif cmd.startswith('pf-'):
            seeb = -np.einsum('...ij,...jk->...ik', self.invert(L[0]), L[1]) / tmp
            return np.einsum('...ij,...jk,...kl->...il', seeb, seeb, L[0]), r'$PF$'
          elif cmd.startswith('tc-'):
            val = L[2] - np.einsum('...ij,...jk,...kl->...il', L[1], self.invert(L[0]), L[1])
            return val / tmp, r'$\kappa_e$'
          elif cmd.startswith('tr-'):
            val = L[2] - np.einsum('...ij,...jk,...kl->...il', L[1], self.invert(L[0]), L[1])
            return self.invert(val / tmp), r'$r$'
          elif cmd.startswith('cb-'):
            return L[0], r'$\sigma_B$'
          elif cmd.startswith('rh-'):
            return np.einsum('...ij,...jkz,...kl->...ilz', self.invert(L[0]), L[1], self.invert(L[0])), r'$R_H$'
          elif cmd.startswith('n-'):
            val  = np.einsum('...ij,...jkz,...kl,...lm->...imz', self.invert(L[0]), L[1], L[2], self.invert(L[0]))
            val -= np.einsum('...ij,...jkz,...kl,...lm->...imz', self.invert(L[0]), L[3], L[0], self.invert(L[0]))
            return -val / tmp, r'$\nu$'
          elif cmd.startswith('muh-'):
            return np.einsum('...ij,...jkz->...ikz', self.invert(L[0]), L[1]), r'$\mu_H$'
          elif cmd.startswith('mut-'):
            return np.einsum('...ij,...jkz->...ikz', self.invert(L[0]), L[1]), r'$\mu_T$'
          else:
            raise IOError('Cannot recognize command')

        # nominal value
        tosave, ylabel = _eval_derived(command, itotal, temp)

        # Zero out T=0 entries for quantities that carry a 1/T prefactor
        # (Seebeck, power factor, thermal conductivity/resistivity, Nernst)
        # Physically these vanish as T->0 (Onsager coefficients decay with at least T^2), but numerically we get NaN from 0/0.
        if command.startswith(('s-', 'pf-', 'tc-', 'tr-', 'n-')):
          t0_mask = (self.temp == 0.0)          # shape: (nT,)
          # Broadcast the mask to the shape of tosave
          slices = (slice(None),) + (np.newaxis,) * (tosave.ndim - 1)
          tosave = np.where(t0_mask[slices], 0.0, tosave)

        # --- std / extremal-envelope for derived quantities ---
        #
        # Two cases:
        #
        # (A) Conductivity (c-): f = L11 (or L11-inter + L11-intra for c-total).
#     The std is simply the std of L11, possibly summed over inter+intra.
#     We propagate std(L11-inter + L11-intra) = std(L11-inter) + std(L11-intra)
#     (linear sum: conservative, appropriate since the two channels are driven
#     by the same disorder realisation and hence positively correlated).
#     Result stored in datastd[command] so the 'raw' display path picks it up.
#
        # (B) All other derived quantities: extremal-envelope over 2^M corners
#     of the +/-sigma hypercube, where M = len(itotal) (the number of SUMMED
#     Onsager arrays after the inter+intra addition for 'total' commands).
#     We build std_itotal to mirror exactly how itotal was built from combined.
#     Result stored in datastd_lo / datastd_hi.
        if self.has_std and not self.owned[command][0]:  # derived quantities only
          if i == 0:
            std_src = self.datastd
          else:
            std_src = self.datastd_spinsum

          reqs_ok = (
            std_src is not None
            and all(ireq in std_src and std_src[ireq] is not None
                    for ireq in requirements)
          )

          if reqs_ok:
            # Build std_itotal in exactly the same way itotal was built:
            # collect per-requirement stds, then sum pairs for 'total' commands.
            std_combined = [np.abs(std_src[ireq]) for ireq in requirements]
            if command.find('total') != -1:
              std_itotal = [
                std_combined[2*k] + std_combined[2*k+1]
                for k in range(len(std_combined)//2)
              ]
            else:
              std_itotal = std_combined

            if command.startswith('c-'):
              # Case A: conductivity -- std is std(L11[total])
              tosave_std = std_itotal[0]
              env_lo = env_hi = None
            else:
              # Case B: extremal envelope over 2^M corners
              import itertools
              M = len(itotal)   # matches len(std_itotal)
              env_lo = None
              env_hi = None
              for signs in itertools.product([-1, 1], repeat=M):
                L_corner = [
                  itotal[k] + signs[k] * std_itotal[k]
                  for k in range(M)
                ]
                val_corner, _ = _eval_derived(command, L_corner, temp)
                val_corner = val_corner.real
                if env_lo is None:
                  env_lo = val_corner.copy()
                  env_hi = val_corner.copy()
                else:
                  env_lo = np.minimum(env_lo, val_corner)
                  env_hi = np.maximum(env_hi, val_corner)
              # Apply same T=0 mask to the envelope bounds
              if command.startswith(('s-', 'pf-', 'tc-', 'tr-', 'n-')):
                slices = (slice(None),) + (np.newaxis,) * (env_lo.ndim - 1)
                env_lo = np.where(t0_mask[slices], 0.0, env_lo)
                env_hi = np.where(t0_mask[slices], 0.0, env_hi)
              tosave_std = None
          else:
            tosave_std = env_lo = env_hi = None
        else:
          tosave_std = env_lo = env_hi = None

        if i==0:
          self.data.update({command:tosave})
          if tosave_std is not None:
            self.datastd[command]    = tosave_std
          if env_lo is not None:
            self.datastd_lo[command] = env_lo
            self.datastd_hi[command] = env_hi
        else:
          self.dataspinsum.update({command:tosave})
          if tosave_std is not None:
            self.datastd_spinsum[command] = tosave_std
          if env_lo is not None:
            self.datastd_lo_ss[command] = env_lo
            self.datastd_hi_ss[command] = env_hi

    if response and self.ndim == 3:
      return ylabel + ' ' + unit
    elif response:
      ''' for 2D and below we have no way of knowing what the unit is exactly '''
      return ylabel
    else:
      ''' no label for non-response quantities '''
      return ''


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

        gapped = h5['.structure/bandgap'+prefix+'/gapped'][()]

        if gapped:
          enev.append(h5['.structure/bandgap'+prefix+'/ene_vband'][()])
          enec.append(h5['.structure/bandgap'+prefix+'/ene_cband'][()])
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



  def outputData(self, command, settings, *args):
    '''
    User interface for lprint.
    Save the data via saveData
    Output the collected data to stdout
    or plot it with matplotlib.

    Standard-deviation output (disorder averaging workflow)
    -------------------------------------------------------
    When self.has_std is True and --nostd is not set:

    * Text mode: two extra columns <command>_std.real and <command>_std.imag
      are appended after the usual real/imag columns.  For derived quantities
      the lower/upper envelope bounds are printed instead (as _lo and _hi).
      Only shown with --fullstd for non-raw quantities.

    * Plot mode: a semi-transparent grey fill_between band is drawn.
      For raw Onsager / conductivity: mean +/- std.
      For derived quantities (--fullstd only): the extremal-envelope
      [datastd_lo, datastd_hi] band is shaded.  This is a rigorous
      outer bound (conservative) rather than a 1-sigma band.
    '''

    if settings.plot:
      import matplotlib.pyplot as plt

    if self.mode == 'mu' and settings.convolve:
      from scipy import signal
      muaxis = self.mu
      nmu = muaxis.shape[0]
      murange = np.max(muaxis) - np.min(muaxis)
      std_translated = float(settings.convolve[0]) * nmu / murange # translate from eV to scipy
      gauss_window = signal.gaussian(nmu,std_translated)
      logger.info('Convoluting with: {} [eV] standard deviation.'.format(settings.convolve[0]))

    ylabel = self.saveData(command, *args)
    self.headerwritten = False

    response = self.owned[command][3]
    magnetic = self.owned[command][4]

    # Decide whether to show std for this command.
    # Raw Onsager and conductivity (c-): shown by default when available.
    # All other derived quantities: only with --fullstd, because the
    # extremal-envelope computation can produce very wide bands for
    # ill-conditioned quantities (thermal conductivity, Nernst, ...).
    _is_raw    = response and self.owned[command][0]
    _is_cond   = command.startswith('c-')
    _show_std  = (
      self.has_std
      and not getattr(settings, 'nostd', False)
      and (_is_raw or _is_cond or getattr(settings, 'fullstd', False))
    )
    # For derived quantities use lo/hi envelope; for raw use datastd directly.
    _use_envelope = _show_std and not (_is_raw or _is_cond)

    '''
    check arguments in more detail
    raise valuerrors if inconsistencies are detected
    '''
    if len(args) > 0:
      if magnetic:
        for icomb in args:
          if len(icomb) == 3 and str(icomb)[0] == '0': # gets added in the main file automatically
            raise ValueError("Invalid directional argument: incorrect argument length [use e.g. xyz, uxxz]")
          if len(icomb) < 3 or len(icomb) > 4:
            raise ValueError("Invalid directional argument: incorrect argument length [use e.g. xyz, uxxz]")
      else:
        for icomb in args:
          if len(icomb) < 2 or len(icomb) > 3:
            raise ValueError("Invalid directional argument: incorrect argument length [use e.g. xx, uxy]")

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

                if settings.diag: # skip non-diagonal elements
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

              outarray *= settings.scale
              if settings.convolve:
                outarray = signal.convolve(outarray, gauss_window, mode='same') / sum(gauss_window)

              # --- extract std / envelope slice for this spin/direction ---
              # std arrays are float64 (real); .real on float64 is a no-op.
              outarray_std = None   # for raw quantities: +/- half-width
              outarray_lo  = None   # for derived: lower envelope
              outarray_hi  = None   # for derived: upper envelope
              if _show_std:
                try:
                  if _use_envelope:
                    # derived: pull from lo/hi dicts
                    _lo_src = self.datastd_lo    if ispin >= 0 else self.datastd_lo_ss
                    _hi_src = self.datastd_hi    if ispin >= 0 else self.datastd_hi_ss
                    if command in _lo_src and _lo_src[command] is not None:
                      _lo = _lo_src[command]
                      _hi = _hi_src[command]
                      if ispin >= 0:
                        outarray_lo = _lo[:,ispin,idir1,idir2] if idir3 is None else _lo[:,ispin,idir1,idir2,idir3]
                        outarray_hi = _hi[:,ispin,idir1,idir2] if idir3 is None else _hi[:,ispin,idir1,idir2,idir3]
                      else:
                        outarray_lo = _lo[:,idir1,idir2] if idir3 is None else _lo[:,idir1,idir2,idir3]
                        outarray_hi = _hi[:,idir1,idir2] if idir3 is None else _hi[:,idir1,idir2,idir3]
                  else:
                    # raw / conductivity: pull from datastd
                    _std_src = self.datastd if ispin >= 0 else self.datastd_spinsum
                    if _std_src and command in _std_src and _std_src[command] is not None:
                      _s = _std_src[command]
                      if ispin >= 0:
                        outarray_std = _s[:,ispin,idir1,idir2].real if idir3 is None else _s[:,ispin,idir1,idir2,idir3].real
                      else:
                        outarray_std = _s[:,idir1,idir2].real if idir3 is None else _s[:,idir1,idir2,idir3].real
                except Exception:
                  outarray_std = outarray_lo = outarray_hi = None

              if settings.plot:
                line, = plt.plot(self.axis, outarray.real, label='{}.real [{}{}]'.format(command, icombdescr, ' - '+self.fname if settings.compare else ''))
                # uncertainty band
                if outarray_std is not None:
                  # raw/conductivity: mean +/- std (grey shading)
                  plt.fill_between(self.axis,
                                   outarray.real - outarray_std,
                                   outarray.real + outarray_std,
                                   alpha=0.25, color='gray',
                                   label=r'{} $\pm1\sigma$ [{}]'.format(command, icombdescr))
                elif outarray_lo is not None:
                  # derived: extremal envelope (grey shading)
                  plt.fill_between(self.axis, outarray_lo, outarray_hi,
                                   alpha=0.25, color='gray',
                                   label='{} envelope [{}]'.format(command, icombdescr))
                if settings.imag: plt.plot(self.axis, outarray.imag, label='{}.imag [{}{}]'.format(command, icombdescr, ' - '+self.fname if settings.compare else ''))
              else:
                if idir3 is None:
                  auxarray = np.zeros((self.nT,3), dtype=int)
                  auxarray[None,:] = np.array([ispin+1,idir1+1,idir2+1], dtype=int)

                  if outarray_std is not None:
                    # raw: append _std.real and _std.imag (zero) columns
                    _scol  = outarray_std[:,None]
                    _zeros = np.zeros_like(_scol)
                    if not self.headerwritten:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            _scol, _zeros, auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i',
                                 header='  {0}{1}, {2:>31}.real, {2:>24}.imag, {2:>24}_std.real, {2:>17}_std.imag,  is id1 id2'.format(
                                   self.axisname, self.axisunit, command))
                      self.headerwritten = True
                    else:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            _scol, _zeros, auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i',
                                 comments='', header='\n')
                  elif outarray_lo is not None:
                    # derived: append _lo and _hi columns
                    if not self.headerwritten:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            outarray_lo[:,None], outarray_hi[:,None], auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i',
                                 header='  {0}{1}, {2:>31}.real, {2:>24}.imag, {2:>24}_env.lo,  {2:>17}_env.hi,   is id1 id2'.format(
                                   self.axisname, self.axisunit, command))
                      self.headerwritten = True
                    else:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            outarray_lo[:,None], outarray_hi[:,None], auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i',
                                 comments='', header='\n')
                  else:
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
                  auxarray = np.zeros((self.nT,4), dtype=int)
                  auxarray[None,:] = np.array([ispin+1,idir1+1,idir2+1,idir3+1], dtype=int)

                  if outarray_std is not None:
                    _scol  = outarray_std[:,None]
                    _zeros = np.zeros_like(_scol)
                    if not self.headerwritten:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            _scol, _zeros, auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i %2i',
                                 header='  {0}{1}, {2:>31}.real, {2:>24}.imag, {2:>24}_std.real, {2:>17}_std.imag,  is id1 id2 id3'.format(
                                   self.axisname, self.axisunit, command))
                      self.headerwritten = True
                    else:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            _scol, _zeros, auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i %2i',
                                 comments='', header='\n')
                  elif outarray_lo is not None:
                    if not self.headerwritten:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            outarray_lo[:,None], outarray_hi[:,None], auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i %2i',
                                 header='  {0}{1}, {2:>31}.real, {2:>24}.imag, {2:>24}_env.lo,  {2:>17}_env.hi,   is id1 id2 id3'.format(
                                   self.axisname, self.axisunit, command))
                      self.headerwritten = True
                    else:
                      np.savetxt(self.textpipe,
                                 np.hstack((self.axis[:,None], outarray.real[:,None], outarray.imag[:,None],
                                            outarray_lo[:,None], outarray_hi[:,None], auxarray)),
                                 fmt='%25.15e %30.18e %30.18e %30.18e %30.18e %5i %2i %2i %2i',
                                 comments='', header='\n')
                  else:
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
      if settings.plot:
        plt.plot(self.axis, outarray, label='{}{}'.format(command, ' - '+self.fname if settings.compare else ''))
      else:
        np.savetxt(self.textpipe, np.hstack((self.axis[:,None], outarray[:,None])), header='{}{}, {}'.format \
        (self.axisname,self.axisunit, self.owned[command][1]))


    if ylabel is not None:
      return ylabel
    else:
      return ''


  def outputList(self, onsager=False):
    '''
    List the internally existing data sets.
    full=False does not list the raw responses (L0 ...)
    full=True  lists all datasets
    '''

    barlength = 80

    print('\n{:<18}  {}'.format('Key', 'Description'))
    print(barlength*u'\u2500')

    # quantities
    for (key, value) in self.owned.items():
      raw_dset, path, description, response, magnetic = value
      if not response:
        print('{:<18}  {}'.format(key, description))
    print(barlength*u'\u2500')

    if onsager:
      # raw responses
      for (key, value) in self.owned.items():
        raw_dset, path, description, response, magnetic = value
        if response and raw_dset and "Boltzmann" not in description:
          print('{:<18}  {}'.format(key, description))
      print(barlength*u'\u2500')
      if self.boltz:
        for (key, value) in self.owned.items():
          raw_dset, path, description, response, magnetic = value
          if response and raw_dset and "Boltzmann" in description:
            print('{:<18}  {}'.format(key, description))
        print(barlength*u'\u2500')
    else:
      # derived responses
      for (key, value) in self.owned.items():
        raw_dset, requirements, description, response, magnetic = value
        if response and not raw_dset and "Boltzmann" not in description:
          print('{:<18}  {}'.format(key, description))
      print(barlength*u'\u2500')
      if self.boltz:
        for (key, value) in self.owned.items():
          raw_dset, requirements, description, response, magnetic = value
          if response and not raw_dset and "Boltzmann" in description:
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
          print('{:<22}  {}'.format(i, self.config[i].decode("utf-8"))) # byte string to string
        except:
          print('{:<22}  {}'.format(i, self.config[i]))

      print(barlength*u'\u2500')


      if self.spins == 1:
        try:
          gapped = h5['.structure/bandgap/gapped'][()]
          gap    = h5['.structure/bandgap/gapsize'][()]
          print('gap [eV]:',gap)
        except:
          print('no gap')
      else:
        try:
          gappedup = h5['.structure/bandgap/up/gapped'][()]
          gapup    = h5['.structure/bandgap/up/gapsize'][()]
          print('up: gap [eV]:',gapup)
        except:
          print('up: no gap')
        try:
          gappeddn = h5['.structure/bandgap/up/gapped'][()]
          gapdn    = h5['.structure/bandgap/dn/gapsize'][()]
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
        print('temperature [K]: {}'.format(self.temp[0]))

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


  def _getResponseCombination(self, key, spinsum, std=False):
    '''
    Get the key numpy array from the file.

    Parameters
    ----------
    key     : str  -- HDF5 path to the sum dataset (e.g. L11/intra/sum)
    spinsum : bool -- if True sum over the spin axis (axis=1)
    std     : bool -- if True read sum_std instead of sum (float64 real array).
                      Returns None if the dataset does not exist.

    Notes on dtypes
    ---------------
    sum datasets are complex128 (Fortran complex(8) via HDF5 compound type).
    sum_std datasets are float64 (real), written by np.std(arr.real, ...) in
    the averaging script.  L11 is purely real; L12/L22 have small imaginary
    parts, so sum_std captures the real-part spread which is what lprint plots.
    '''

    actual_key = key.replace('/sum', '/sum_std') if std else key

    with h5py.File(self.fname,'r') as h5:
      if actual_key not in h5:
        return None   # std dataset absent -- caller must check has_std
      outputarray = h5['{}'.format(actual_key)][()]

    if spinsum:
      data = np.sum(outputarray, axis=1)
      return data

    else:
      # NOTE: sum_std is float64; the complex128 promotion below still works
      # because float64 assigns cleanly into complex128 with zero imaginary part.
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
    self.boltz = False
    for (key, value) in self.datasets.items():
      raw_dset, requirements, description, response, magnetic = value

      if raw_dset:
        continue
      for ireq in requirements:
        if ireq not in self.owned:
          break
      else:
        self.owned.update({key:value})
        if "Boltzmann" in description:
          self.boltz = True

    if self.boltz:
      logger.debug('Detected Boltzmann container.')

    if logger.isEnabledFor(logging.DEBUG):
      print('Owned datasets:')
      for key, value in self.owned.items():
        print('{0:<20}'.format(key), value)

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
        logger.debug('Trying to open: {}'.format(fi))
        with h5py.File(fi,'r') as hfi:
          if hfi.attrs['identifier'].decode("utf-8") == 'LRTCoutput':
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

    Additionally sets self.has_std (bool): True when the file contains
    sum_std datasets produced by the disorder averaging workflow.
    '''

    try:
      with h5py.File(self.fname,'r') as hfi:
        # we get a byte string here
        # hence we have to decode it to utf-8
        if hfi.attrs['identifier'].decode("utf-8") != 'LRTCoutput':
          raise IOError('{} is not an LRTC output file.'.format(self.fname))
        if hfi['.quantities'].attrs['mode'].decode("utf-8") == 'temp':
          self.mode = 'temp'
        elif hfi['.quantities'].attrs['mode'].decode("utf-8") == 'mu':
          self.mode = 'mu'

        self.spins = hfi['.structure/ispin'][()]
        self.ndim = hfi['.unitcell/ndim'][()]
        self.dims = hfi['.unitcell/dims'][()]
        self.dimmask2 = np.logical_and(self.dims[:,None], self.dims[None,:])
        self.dimmask3 = np.logical_and(np.logical_and(self.dims[:,None], self.dims[None,:])[:,:,None], self.dims[None,None,:])
        print('#   File: {} - Run mode: {} - {} crystal directions: {}'.format(self.fname, self.mode, self.ndim, np.array(["x","y","z"])[self.dims]))

        # Detect whether this file was produced by the disorder averaging
        # workflow (sum_std datasets present alongside the usual sum datasets).
        self.has_std = False
        for _probe in ['L11/intra/sum_std', 'L11/inter/sum_std',
                       'L11B/intra/sum_std', 'L11/intraBoltzmann/sum_std']:
          if _probe in hfi:
            self.has_std = True
            break
        if self.has_std:
          print('#   Standard-deviation data detected (disorder averaging workflow).')
    except:
      raise IOError('{} is not an LRTC output file.'.format(self.fname))

  def _get_axis(self, altaxis):
    '''
    Get the temperature and inverse temperature axis
    Also save the number of temperature steps
    '''
    try:
      with h5py.File(self.fname, 'r') as h5:
        self.temp    = h5['.quantities/tempAxis'][()]
        self.invtemp = h5['.quantities/betaAxis'][()]
        self.carrier = h5['.quantities/carrier'][()]
        self.mu      = h5['.quantities/mu'][()]
        self.mudft   = h5['.structure/mudft'][()]
        self.nT      = self.temp.shape[0]
    except:
      raise IOError('Incomplete LRTC output file')

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

  def outputDOS(self, plot, broadening=0.02, npoints=1001, emin=None, emax=None):
    '''
    Print DOS/NOS
    '''

    import matplotlib.pyplot as plt

    with h5py.File(self.fname,'r') as h5:
      mu  = h5['.structure/mudft'][()]
      spins = h5['.structure/ispin'][()]
      weights = h5['.structure/weights'][()]

      if spins==1:
        ene = h5['.structure/energies'][()]
        dosaxis, dos, nos = calcDOS(ene, weights, npoints=npoints, gamma=broadening, windowsize=1.1, emin=emin, emax=emax)
        del ene
      else:
        eneup = h5['.structure/energies/up'][()]
        dosaxisup, dosup, nosup = calcDOS(eneup, weights, npoints=npoints, gamma=broadening, windowsize=1.1, emin=emin, emax=emax)
        del eneup
        enedn = h5['.structure/energies/dn'][()]
        dosaxisdn, dosdn, nosdn = calcDOS(enedn, weights, npoints=npoints, gamma=broadening, windowsize=1.1, emin=emin, emax=emax)
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

    ''' sub select values for valid directions '''
    if self.ndim < 3 and self.ndim > 0: # 1D and 2D
      # select direction combinations
      if spins:
        outputmasked = data[:,:,self.dimmask2].copy()
        outputmasked = outputmasked.reshape((data.shape[0],data.shape[1],self.ndim,self.ndim))
      else:
        outputmasked = data[:,self.dimmask2].copy()
        outputmasked = outputmasked.reshape((data.shape[0],self.ndim,self.ndim))
    elif self.ndim == 3: # 3D
      outputmasked = data
    elif self.ndim == 0: # 0D
      outputmasked = data
      # there cannot be mixed directions in 0D per definition -> override 0
      outputmasked[...,[0,0,1,1,2,2],[1,2,0,2,0,1]] = 0.0

    ''' invert provided masked values '''
    if self.ndim > 0:
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
    else: # 0D : avoid inversion routines
      inverted = np.full_like(outputmasked, fill_value=0., dtype=np.complex128)
      if spins:
        for i in range(outputmasked.shape[0]):
          for j in range(outputmasked.shape[1]):
            inverted[i,j,[0,1,2],[0,1,2]] = 1. / outputmasked[i,j,[0,1,2],[0,1,2]]
      else:
        for i in range(outputmasked.shape[0]):
          inverted[i,[0,1,2],[0,1,2]] = 1. / outputmasked[i,[0,1,2],[0,1,2]]

      ''' hotfix for weird behavior in 0D '''
      inverted = np.nan_to_num(inverted, nan=0.0)

    ''' for 2D and 3D : expand inverted values back to 3 x 3 array '''
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
