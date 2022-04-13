import numpy as np
'''
2D Monoclinic yz
'''

strsym = 'C2-yz'
nsym = 2
symop = np.zeros((nsym,3,3), dtype=np.float64)
invsymop = np.zeros_like(symop, dtype=np.float64)

symop[0,:,:] = np.array([[ 1, 0, 0],\
                         [ 0, 1, 0],\
                         [ 0, 0, 1]])
symop[1,:,:] = np.array([[ 1, 0, 0],\
                         [ 0,-1, 0],\
                         [ 0, 0,-1]])

for isym in range(nsym):
  invsymop[isym] = np.linalg.inv(symop[isym])

if __name__ == '__main__':
  symops = set([])
  for i in range(nsym):
    temp = tuple(['{0:1f}'.format(s) for s in symop[i].flatten()])
    symops.add(temp)
  print('number of ineqvivalent symop: ',len(symops))
