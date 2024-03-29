import numpy as np
'''
3D Cubic tight-binding lattices a=b=c
'''

strsym = 'OH'
nsym = 48
symop = np.zeros((nsym,3,3), dtype=np.float64)
invsymop = np.zeros_like(symop, dtype=np.float64)

symop[0,:,:]  = np.array([[-1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0,-1]])
symop[1,:,:]  = np.array([[-1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0, 1]])
symop[2,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 0,-1],\
                          [ 0,-1, 0]])
symop[3,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 0,-1],\
                          [ 0, 1, 0]])
symop[4,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 0, 1],\
                          [ 0,-1, 0]])
symop[5,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 0, 1],\
                          [ 0, 1, 0]])
symop[6,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0,-1]])
symop[7,:,:]  = np.array([[-1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0, 1]])
symop[8,:,:]  = np.array([[ 0,-1, 0],\
                          [-1, 0, 0],\
                          [ 0, 0,-1]])
symop[9,:,:]  = np.array([[ 0,-1, 0],\
                          [-1, 0, 0],\
                          [ 0, 0, 1]])
symop[10,:,:] = np.array([[ 0,-1, 0],\
                          [ 0, 0,-1],\
                          [-1, 0, 0]])
symop[11,:,:] = np.array([[ 0,-1, 0],\
                          [ 0, 0,-1],\
                          [ 1, 0, 0]])
symop[12,:,:] = np.array([[ 0,-1, 0],\
                          [ 0, 0, 1],\
                          [-1, 0, 0]])
symop[13,:,:] = np.array([[ 0,-1, 0],\
                          [ 0, 0, 1],\
                          [ 1, 0, 0]])
symop[14,:,:] = np.array([[ 0,-1, 0],\
                          [ 1, 0, 0],\
                          [ 0, 0,-1]])
symop[15,:,:] = np.array([[ 0,-1, 0],\
                          [ 1, 0, 0],\
                          [ 0, 0, 1]])
symop[16,:,:] = np.array([[ 0, 0,-1],\
                          [-1, 0, 0],\
                          [ 0,-1, 0]])
symop[17,:,:] = np.array([[ 0, 0,-1],\
                          [-1, 0, 0],\
                          [ 0, 1, 0]])
symop[18,:,:] = np.array([[ 0, 0,-1],\
                          [ 0,-1, 0],\
                          [-1, 0, 0]])
symop[19,:,:] = np.array([[ 0, 0,-1],\
                          [ 0,-1, 0],\
                          [ 1, 0, 0]])
symop[20,:,:] = np.array([[ 0, 0,-1],\
                          [ 0, 1, 0],\
                          [-1, 0, 0]])
symop[21,:,:] = np.array([[ 0, 0,-1],\
                          [ 0, 1, 0],\
                          [ 1, 0, 0]])
symop[22,:,:] = np.array([[ 0, 0,-1],\
                          [ 1, 0, 0],\
                          [ 0,-1, 0]])
symop[23,:,:] = np.array([[ 0, 0,-1],\
                          [ 1, 0, 0],\
                          [ 0, 1, 0]])
symop[24,:,:] = np.array([[ 0, 0, 1],\
                          [-1, 0, 0],\
                          [ 0,-1, 0]])
symop[25,:,:] = np.array([[ 0, 0, 1],\
                          [-1, 0, 0],\
                          [ 0, 1, 0]])
symop[26,:,:] = np.array([[ 0, 0, 1],\
                          [ 0,-1, 0],\
                          [-1, 0, 0]])
symop[27,:,:] = np.array([[ 0, 0, 1],\
                          [ 0,-1, 0],\
                          [ 1, 0, 0]])
symop[28,:,:] = np.array([[ 0, 0, 1],\
                          [ 0, 1, 0],\
                          [-1, 0, 0]])
symop[29,:,:] = np.array([[ 0, 0, 1],\
                          [ 0, 1, 0],\
                          [ 1, 0, 0]])
symop[30,:,:] = np.array([[ 0, 0, 1],\
                          [ 1, 0, 0],\
                          [ 0,-1, 0]])
symop[31,:,:] = np.array([[ 0, 0, 1],\
                          [ 1, 0, 0],\
                          [ 0, 1, 0]])
symop[32,:,:] = np.array([[ 0, 1, 0],\
                          [-1, 0, 0],\
                          [ 0, 0,-1]])
symop[33,:,:] = np.array([[ 0, 1, 0],\
                          [-1, 0, 0],\
                          [ 0, 0, 1]])
symop[34,:,:] = np.array([[ 0, 1, 0],\
                          [ 0, 0,-1],\
                          [-1, 0, 0]])
symop[35,:,:] = np.array([[ 0, 1, 0],\
                          [ 0, 0,-1],\
                          [ 1, 0, 0]])
symop[36,:,:] = np.array([[ 0, 1, 0],\
                          [ 0, 0, 1],\
                          [-1, 0, 0]])
symop[37,:,:] = np.array([[ 0, 1, 0],\
                          [ 0, 0, 1],\
                          [ 1, 0, 0]])
symop[38,:,:] = np.array([[ 0, 1, 0],\
                          [ 1, 0, 0],\
                          [ 0, 0,-1]])
symop[39,:,:] = np.array([[ 0, 1, 0],\
                          [ 1, 0, 0],\
                          [ 0, 0, 1]])
symop[40,:,:] = np.array([[ 1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0,-1]])
symop[41,:,:] = np.array([[ 1, 0, 0],\
                          [ 0,-1, 0],\
                          [ 0, 0, 1]])
symop[42,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 0,-1],\
                          [ 0,-1, 0]])
symop[43,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 0,-1],\
                          [ 0, 1, 0]])
symop[44,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 0, 1],\
                          [ 0,-1, 0]])
symop[45,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 0, 1],\
                          [ 0, 1, 0]])
symop[46,:,:] = np.array([[ 1, 0, 0],\
                          [ 0, 1, 0],\
                          [ 0, 0,-1]])
symop[47,:,:] = np.array([[ 1, 0, 0],\
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
