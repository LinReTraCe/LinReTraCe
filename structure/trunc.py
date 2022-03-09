#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import logging
logger = logging.getLogger(__name__)

import numpy as np

from structure.dft   import DFTcalculation

# rather bare-bone Exception classes
class TruncationError(IOError):
  def __init__(self, message):
    super(TruncationError, self).__init__()
    self.message = message
  def __str__(self):
    return self.message


def truncate(dftcalc, btpinterp, energy1, energy2, absolute=False):
  '''
  Use the provided energy range (absolute or relative to the chemical potential,
  calculate the common bands, and truncate the energy window for the
  energies and optical elements. Also reduce the charge in the system accordingly.
  This method can only be applied to DFT calculations and not for models.
  '''

  if not isinstance(dftcalc, DFTcalculation):
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
        cutofftop.append(iband)
        break


  # take the lowest and highest bands
  # remember: these are python indices
  bmin = min(cutoffbottom)
  bmax = max(cutofftop)

  # number of cut-off bands
  cutlow  = bmin
  cuthigh = dftcalc.energyBandMax - bmax - 1

  logger.debug('')
  logger.debug('spin band min max (# == included)')
  for ispin in range(dftcalc.spins):
    for iband in range(dftcalc.energyBandMax):
      ene = dftcalc.energies[ispin][:,iband]
      bandmin = np.min(ene)
      bandmax = np.max(ene)
      logger.debug('{:4} {:4} {:18.10f} {:18.10f} {}'.format(ispin,iband,bandmin,bandmax, '#' if iband in range(bmin,bmax+1) else '' ))
    logger.debug("")


  oldbandmax = dftcalc.energyBandMax # for verbose output
  oldcharge  = dftcalc.charge        # for verbose output

  # truncate the energies
  for ispin in range(dftcalc.spins):
    dftcalc.energies[ispin] = dftcalc.energies[ispin][:,bmin:bmax+1] # +1 because of python ranges
    if dftcalc.gapped[ispin]:
      dftcalc.cb[ispin] -= bmin
      dftcalc.vb[ispin] -= bmin
    dftcalc.charge -= bmin*dftcalc.weightsum # remove charge equivalent to the cut-off bands

  dftcalc.energyBandMax = bmax - bmin + 1

  # output truncation information (it rhymes)
  logger.info("Truncating window: {} - {} [eV]".format(truncmin, truncmax))
  logger.info("Truncating procedure resulted in:")
  logger.info("   range:  1 - {:3<} ---> {:3>} - {:3<}".format(oldbandmax, bmin+1, bmax+1))
  logger.info("   bands:    {:3>}   --->   {:3<}".format(oldbandmax, dftcalc.energyBandMax))
  logger.info("   charge:  {:5.1f} ---> {:5.1f}".format(oldcharge, dftcalc.charge))


  # truncate possible derivatives
  # this is straight-forward: same treatment as the energies
  if btpinterp is not None:
    for ispin in range(btpinterp.spins):
      btpinterp.velocities[ispin]  = btpinterp.velocities[ispin][:,bmin:bmax+1]
      btpinterp.curvatures[ispin]  = btpinterp.curvatures[ispin][:,bmin:bmax+1]
      btpinterp.opticalDiag[ispin] = btpinterp.opticalDiag[ispin][:,bmin:bmax+1] # peierls
      btpinterp.opticalBandMin = 0
      btpinterp.opticalBandMax = dftcalc.energyBandMax


  # general case of truncation procedure
  # fuck this shit
  # 0                                            energyBandMax
  # |--------------------------------------------| energies
  #      bmin                               bmax+1
  #      |----------------------------------|      new energies
  # optical elements could be
  #   |----------------------------------------|   outside
  #        |-----------------------------|         inside
  #   |----------------------------------|         one side
  #        |-----------------------------------|   one side
  # all 4 case have to be cut appropiately

  # truncate optical elements
  if dftcalc.opticdiag:
    # offset from optical band range to old energy band range
    offsetlow  = dftcalc.opticalBandMin
    offsethigh = oldbandmax - dftcalc.opticalBandMax
    opticalinterval = dftcalc.opticalBandMax - dftcalc.opticalBandMin # original optical interval

    if cutlow > offsetlow:
      optstart = cutlow - offsetlow
    else:
      optstart = 0

    if cuthigh > offsethigh:
      optend = opticalinterval - (cuthigh - offsethigh)
    else:
      optend = opticalinterval

    opticalinterval = optend - optstart

    if optstart > 0: # lower truncation -> same starting index
      dftcalc.opticalBandMin = 0
    else:
      dftcalc.opticalBandMin -= cutlow

    dftcalc.opticalBandMax = dftcalc.opticalBandMin + opticalinterval

    # perform the truncation of the arrays
    for ispin in range(dftcalc.spins):
      dftcalc.opticalDiag[ispin] = dftcalc.opticalDiag[ispin][:,optstart:optend,:]
      if dftcalc.opticfull:
        dftcalc.opticalMoments[ispin] = dftcalc.opticalMoments[ispin][:,optstart:optend,optstart:optend,:]
