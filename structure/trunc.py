#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np

from structure.dft   import DftCalculation

# rather bare-bone Exception classes
class TruncationError(IOError):
  def __init__(self, message):
    super(TruncationError, self).__init__()
    self.message = message
  def __str__(self):
    return self.message


def truncate(dftcalc, velcalc, energy1, energy2, absolute=False):
  '''
  Use the provided energy range (absolute or relative to the chemical potential,
  calculate the common bands, and truncate the energy window for the
  energies and optical elements. Also reduce the charge in the system accordingly.
  This method can only be applied to DFT calculations and not for models.
  '''

  if not isinstance(dftcalc, DftCalculation):
    raise IOError('Provided DFT calculation object for truncation is illegal.')

  # make relative energies absolute
  if not absolute:
    energy1 += dftcalc.mu
    energy2 += dftcalc.mu

  # to avoid user errors
  truncmin = min(energy1,energy2)
  truncmax = max(energy1,energy2)

  # raise an error if the provided truncation window
  # is non-sense
  if truncmax < dftcalc.mu \
     or truncmin > dftcalc.mu:
    raise TruncationError('Truncation window does not include Fermi level.')


  # find the band identifier
  # for the bands we want to retain for the provided energy window

  cutoffbottom = []
  for ispin in range(dftcalc.spins):
    for iband in range(dftcalc.energyBandMax):
      ene = dftcalc.energies[ispin][:,iband]
      bandmax = np.max(ene)

      if bandmax < truncmin:
        continue
      else:
        cutoffbottom.append(iband)
        break

  cutofftop    = []
  for ispin in range(dftcalc.spins):
    for iband in range(dftcalc.energyBandMax-1,-1,-1):
      ene = dftcalc.energies[ispin][:,iband]
      bandmin = np.min(ene)

      if bandmin > truncmax:
        continue
      else:
        cutofftop.append(iband+1)
        break

  # take the lowest and highest bands
  bmin = min(cutoffbottom)
  bmax = max(cutofftop)

  enemin = []
  enemax = []
  for ispin in range(dftcalc.spins):
    enemin.append(min(dftcalc.energies[ispin][:,bmin]))
    enemax.append(max(dftcalc.energies[ispin][:,bmax-1]))

  # number of cut-off bands
  cutlow  = bmin
  cuthigh = dftcalc.energyBandMax - bmax

  logger.debug('')
  logger.debug('spin band min max (# == included)')
  for ispin in range(dftcalc.spins):
    for iband in range(dftcalc.energyBandMax):
      ene = dftcalc.energies[ispin][:,iband]
      bandmin = np.min(ene)
      bandmax = np.max(ene)
      logger.debug('{:4} {:4} {:18.10f} {:18.10f} {}'.format(ispin,iband,bandmin,bandmax, '#' if iband in range(bmin,bmax) else '' ))
    logger.debug("")


  oldbandmax = dftcalc.energyBandMax
  oldcharge  = dftcalc.charge

  # truncate the energies
  for ispin in range(dftcalc.spins):
    dftcalc.energies[ispin] = dftcalc.energies[ispin][:,bmin:bmax]
    if dftcalc.gapped[ispin]:
      dftcalc.cb[ispin] -= bmin
      dftcalc.vb[ispin] -= bmin
    dftcalc.charge -= bmin*dftcalc.weightsum # remove charge equivalent to the cut-off bands

  dftcalc.energyBandMax = bmax - bmin

  # output truncation information (it rhymes)
  logger.info("Truncating window: {} - {} [eV]".format(truncmin, truncmax))
  logger.info("Truncating procedure resulted in:")
  logger.info("   window:  {} - {} [eV]".format(min(enemin),max(enemax)))
  logger.info("   [ bands that touch the limits get fully included ]")
  logger.info("   range:  1 - {:3<} ---> {:3>} - {:3<}".format(oldbandmax, bmin+1, bmax))
  logger.info("   bands:    {:3>}   --->   {:3<}".format(oldbandmax, dftcalc.energyBandMax))
  logger.info("   charge:  {:5.1f} ---> {:5.1f}".format(oldcharge, dftcalc.charge))


  # truncate possible derivatives
  # this is straight-forward: same treatment as the energies as we get the derivatives
  # and curvatures on the full energy range
  if velcalc is not None:
    for ispin in range(velcalc.spins):
      velcalc.velocities[ispin]  = velcalc.velocities[ispin][:,bmin:bmax]
      velcalc.curvatures[ispin]  = velcalc.curvatures[ispin][:,bmin:bmax]
      velcalc.opticalDiag[ispin] = velcalc.opticalDiag[ispin][:,bmin:bmax]
      velcalc.BopticalDiag[ispin] = velcalc.BopticalDiag[ispin][:,bmin:bmax]
      velcalc.opticalBandMin = 0
      velcalc.opticalBandMax = dftcalc.energyBandMax


  # general case of truncation procedure
  # energy truncation
  # 0                     :                      energyBandMax
  # |--------------------------------------------| energies
  # ->---->--->-
  #      bmin             :                 bmax
  #      |----------------------------------|      new energies
  #
  #
  # optical are on the original energy energy range but could look like
  #   |----------------------------------------|   outside
  #        |-----------------------------|         inside
  #   |----------------------------------|         one side out one in
  #        |-----------------------------------|   one side in one out
  # all 4 case have to be cut appropiately

  # nb: opticalbandmin and opticalbandmax are still in python notation from the input
  #  i.e. energies[:energyBandMax] -- optical[opticalBandMin:opticalBandMax]

  # truncate optical elements
  if dftcalc.opticdiag:

    opticalinterval = dftcalc.opticalBandMax - dftcalc.opticalBandMin # original optical interval

    # offset from optical band range to old energy band range
    if cutlow > dftcalc.opticalBandMin: # we cut into the lower optical interval
      opticalinterval -= (cutlow - dftcalc.opticalBandMin)
      optstart = (cutlow - dftcalc.opticalBandMin) # array access
      dftcalc.opticalBandMin = 0
    else:
      dftcalc.opticalBandMin -= cutlow
      optstart = 0

    if cuthigh > (oldbandmax - dftcalc.opticalBandMax): # we cut into the upper optical interval
      opticalinterval -= (cuthigh - (oldbandmax-dftcalc.opticalBandMax))

    optend = optstart + opticalinterval
    dftcalc.opticalBandMax = dftcalc.opticalBandMin + opticalinterval

    # perform the truncation of the arrays
    for ispin in range(dftcalc.spins):
      dftcalc.opticalDiag[ispin] = dftcalc.opticalDiag[ispin][:,optstart:optend,:]
      if dftcalc.opticfull:
        dftcalc.opticalMoments[ispin] = dftcalc.opticalMoments[ispin][:,optstart:optend,optstart:optend,:]
