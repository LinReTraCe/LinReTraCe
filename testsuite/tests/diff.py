#! /usr/bin/env python

import numpy as np
import h5py
import sys

def main():
  with h5py.File('testsuite/tests/Si_output.hdf5','r') as h5:
    temp = h5['.quantities/tempAxis'][()]
    L11 = h5['L11/intra/sum'][:,0,0,0]
    L12 = h5['L12/intra/sum'][:,0,0,0]
    cond_out = L11
    seeb_out = -L12/L11/temp

  with h5py.File('testsuite/tests/Si_output_comparison.hdf5','r') as h5:
    temp = h5['.quantities/tempAxis'][()]
    L11 = h5['L11/intra/sum'][:,0,0,0] # spin up xx
    L12 = h5['L12/intra/sum'][:,0,0,0] # spin up xx
    cond_comp = L11
    seeb_comp = -L12/L11/temp

  close_conductivity = np.allclose(cond_out, cond_comp, rtol=1e-4, atol=1e-6)
  close_seebeck      = np.allclose(seeb_out, seeb_comp, rtol=1e-4, atol=1e-6)

  print('Conductivity check: ', close_conductivity)
  print('Seebeeck     check: ', close_seebeck)

  if close_conductivity and close_seebeck:
    return 0
  else:
    return 1

if __name__ == '__main__':
  val = main()
  sys.exit(val)
