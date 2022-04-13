import numpy as np
'''
3D Tetragonal tight-binding lattices a=c!=b
a=b is important here, otherwise the symmetry operations below would be different
'''

strsym = 'D4H-xz'
nsym = 16
symop = np.zeros((nsym,3,3), dtype=np.float64)
invsymop = np.zeros_like(symop, dtype=np.float64)

symop[0,:,:]  = np.array([[-1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0,-1]])
symop[1,:,:]  = np.array([[-1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0, 1]])
symop[2,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0,-1]])
symop[3,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0, 1]])
symop[4,:,:]  = np.array([[ 0, 0,-1],\
                          [ 0,-1, 0],\
                          [-1, 0, 0]])
symop[5,:,:]  = np.array([[ 0, 0, 1],\
                          [ 0,-1, 0],\
                          [-1, 0, 0]])
symop[6,:,:]  = np.array([[ 0, 0,-1],\
                          [ 0, 1, 0],\
                          [-1, 0, 0]])
symop[7,:,:]  = np.array([[ 0, 0, 1],\
                          [ 0, 1, 0],\
                          [-1, 0, 0]])
symop[8,:,:]  = np.array([[ 0, 0,-1],\
                          [ 0,-1, 0],\
                          [ 1, 0, 0]])
symop[9,:,:]  = np.array([[ 0, 0, 1],\
                          [ 0,-1, 0],\
                          [ 1, 0, 0]])
symop[10,:,:] = np.array([[ 0, 0,-1],\
                          [ 0, 1, 0],\
                          [ 1, 0, 0]])
symop[11,:,:] = np.array([[ 0, 0, 1],\
                          [ 0, 1, 0],\
                          [ 1, 0, 0]])
symop[12,:,:] = np.array([[ 1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0,-1]])
symop[13,:,:] = np.array([[ 1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0, 1]])
symop[14,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0,-1]])
symop[15,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0, 1]])


for isym in range(nsym):
  invsymop[isym] = np.linalg.inv(symop[isym])


if __name__ == '__main__':
  symops = set([])
  for i in range(nsym):
    temp = tuple(['{0:1f}'.format(s) for s in symop[i].flatten()])
    symops.add(temp)
  print('number of ineqvivalent symop: ',len(symops))
