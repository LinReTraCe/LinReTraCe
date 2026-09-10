#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import math

rydberg2eV    = 13.605662285137
hartree2eV    = 2*rydberg2eV
# Exact inverse of hartree2eV.  Used when handing energies to BoltzTraP2,
# which works in Hartree, so that the eV -> Ha -> eV round trip is the
# identity rather than accurate to ~1e-7 (see KI-05).
eV2hartree    = 1.0/hartree2eV
bohr2angstrom = 0.529177
hbarJs        = 1.054571817e-34
hbareVs       = 6.582119569e-16
mekg          = 9.10938356e-31
echargeC      = 1.6021766208e-19

# Boltzmann constant in eV / K.
# Same value as the Fortran parameter kB in src_linretrace/params.f90 and as
# the local definitions in scattering/*.py and scripts/*.py.  All internal
# energies in LinReTraCe are in eV, hence beta = 1 / (kB_eV * T[K]) in eV^-1.
kB_eV         = 8.61733034e-5

w2kmom        = (1/bohr2angstrom * 1e20 * hbarJs * hbareVs / mekg) ** 2

if __name__ == '__main__':
  print(w2kmom)
