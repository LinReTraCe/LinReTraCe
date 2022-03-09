#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import math

rydberg2eV    = 13.605662285137
hartree2eV    = 2*rydberg2eV
bohr2angstrom = 0.529177
hbarJs        = 1.054571817e-34
hbareVs       = 6.582119569e-16
mekg          = 9.10938356e-31
echargeC      = 1.6021766208e-19

w2kmom        = (1/bohr2angstrom * 1e20 * hbarJs * hbareVs / mekg) ** 2

# print(math.pi * echargeC**2 * hbarJs * hbareVs**2 / mekg**2 *1e50 / bohr2angstrom**2)
# print(math.pi * echargeC / hbareVs * 1e10 * w2kmom)

# print(w2kmom)

# print(math.pi * echargeC / hbareVs * 1e10)
