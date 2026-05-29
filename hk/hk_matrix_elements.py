import sympy
from sympy import *
from sympy.physics.matrices import *
import numpy as np
import json


#defines the symbols for Sympy 
symbol_names = ["kx", "ky", "kz"]
symbols_dict = {name: symbols(name) for name in symbol_names} #define the symbols 

kx=symbols_dict["kx"]
ky=symbols_dict["ky"]
kz=symbols_dict["kz"]

class HkInputError(Exception):
    pass

def Hk(d0,d1,d2,d3): #builds a 2x2 Hamiltonian using Pauli matrices 
    return d0*Matrix([[1,0],[0,1]])+d1*msigma(1)+d2*msigma(2)+d3*msigma(3)

def dHk(Hk): #takes the derivatives of a given k-space Hamiltonian in each direction
    dHkx=diff(Hk,kx)
    dHky=diff(Hk,ky)
    dHkz=diff(Hk,kz)
    return dHkx,dHky,dHkz

def gen_kmesh(b1,b2,b3,nkx,nky,nkz): #generates a reducible k-space mesh 

    k_points=[]

    dx=np.linspace(0,1,nkx,endpoint=False)
    dy=np.linspace(0,1,nky,endpoint=False)
    dz=np.linspace(0,1,nkz,endpoint=False)

    for i in dx:
        for j in dy:
            for k in dz:
                k_points.append(i*b1+j*b2+k*b3)

    return np.array(k_points)

def mapping_func(n,N):
    return floor(n/N), n%N #maps an input n to a given i,j coordinate of a matrix of 'size' N where N is the total number of components 

def write_ham_to_file(Hk,n_bands,kpoints,filename): #writes the Hamiltonian to a file which can be read by w2dynamics DMFT code
    
    nk=len(Hk)

    with open(filename,'w') as A:

        A.write(f"{nk:10d}{n_bands:6d}{n_bands:6d}\n") #write first line

        for ik in range(nk):
            
            Kx=kpoints[ik][0]
            Ky=kpoints[ik][1]
            Kz=kpoints[ik][2]

            A.write(f"{Kx:12.8f}{Ky:12.8f}{Kz:12.8f}\n") #write kpoints

            for band1 in range(n_bands):
                for band2 in range(n_bands):

                    Hk_real=Hk[ik,band1,band2].real
                    Hk_imag=Hk[ik,band1,band2].imag

                    A.write(f"{Hk_real:21.12E}") #write real part of Hk
                    A.write(f"{Hk_imag:21.12E}") #write imag part of Hk
            A.write("\n")

def Matrix_Elements_Gen(hk,b1,b2,b3,nkx,nky,nkz,to_file=False): #reads in matrix elements for a general size H(k)

    nk=nkx*nky*nkz #number of k-points

    matrix_size=len(hk) #finds the size of the Matrix object provided. This is the total number of elements inside the matrix

    n_bands=int(np.sqrt(matrix_size)) #number of bands of the system

    dhkx,dhky,dhkz=dHk(hk) #takes the symbolic derivative of the k-space Hamiltonian 

    #arrays to hold the lambdified functions
    hk_funcs=[]
    dhkx_funcs=[]
    dhky_funcs=[]
    dhkz_funcs=[]
    
    #generates the reducible k-space mesh on which H(k) and dH(k) will be evaulated 
    kpoints=gen_kmesh(b1,b2,b3,nkx,nky,nkz)

    #initialse the arrays to be loaded into later 
    Hk_array=np.zeros((nk,n_bands,n_bands),dtype=complex) # kpoint, index1, index2
    dHk_array=np.zeros((nk,n_bands,n_bands,3),dtype=complex) # kpoint, index1, index2, direction 
    
    energy_array=np.zeros((2,nk,n_bands)) #spin, kpoint, bands
    full_array=np.zeros((2,nk,n_bands,n_bands,9))#spin, kpoint, band1, band2, direction (note directions 6/7/8 are reserved for possible imaganiary inter-unit cell hoppings- so consistent with tb approach)
    diagonal_array=np.zeros((2,nk,n_bands,9))#spin, kpoint, bands, direction
    berry_array=np.zeros((2,nk,n_bands,n_bands,3)) #spin, kpoint, band1, band2, direction

    #take each component from the Hamiltonian matrix Hij and turn it into a function which can be accessed for evaluation and add to a list
    for i in range(matrix_size):
        hk_funcs.append(lambdify([kx,ky,kz],hk[i],'math'))
        dhkx_funcs.append(lambdify([kx,ky,kz],dhkx[i],'math'))
        dhky_funcs.append(lambdify([kx,ky,kz],dhky[i],'math'))
        dhkz_funcs.append(lambdify([kx,ky,kz],dhkz[i],'math'))

    for k in range(nk):
        for n in range(matrix_size):
            i,j=mapping_func(n,n_bands) #mapping from the index in a list 'n' to the indices of a matrix (i,j)

            #extract the nth component from the function list to allow for evaluation at current k-point 
            fhk=hk_funcs[n] 
            fhkx=dhkx_funcs[n]
            fhky=dhky_funcs[n]
            fhkz=dhkz_funcs[n]

            #evaluate the Hamiltonian and its derivatives at each k-point

            Hk_array[k,i,j]=fhk(kpoints[k][0],kpoints[k][1],kpoints[k][2])
            dHk_array[k,i,j,0]=fhkx(kpoints[k][0],kpoints[k][1],kpoints[k][2]) #kx derivative
            dHk_array[k,i,j,1]=fhky(kpoints[k][0],kpoints[k][1],kpoints[k][2]) #ky derivative
            dHk_array[k,i,j,2]=fhkz(kpoints[k][0],kpoints[k][1],kpoints[k][2]) #kz derivative 

    if to_file==True: #if write_to_file is true then the Hamiltonian is saved to txt file 
        write_ham_to_file(Hk_array,n_bands,kpoints,"hamk.txt")

    ek, U = np.linalg.eigh(Hk_array) #finds the eigenvalues 'ek' and the eigenvectors 'U' of the Hamiltonian

    for ik in range(nk): #reorders the eignevalues and eignevectors for consistency 
        ekk, Uk = ek[ik,:], U[ik,:,:]
        idx = ekk.argsort()
        ek[ik,:] = ekk[idx]
        U[ik,:,:] = Uk[:,idx]

    energy_array[0,...]=ek.real #load in the energy

    Uinv = np.linalg.inv(U) #find the inverse of the eigenvector U(k)

    vel = np.einsum('kab,kbci,kcd->kadi',Uinv,dHk_array,U) #rotates the dH(k) matrix into the band basis by U^(-1)(k)dH(k)U(k)

    #indices used for symmetric matrix elements- order of xx/yy/zz/xy/xz/yz
    index1=np.array([0,1,2,0,0,1]) 
    index2=np.array([0,1,2,1,2,2]) 

    #indices used for the anti-symmertic Berry curvature matrix elements- order of xy/xz/yz
    index3=np.array([0,0,1]) 
    index4=np.array([1,2,2])


    #load in the diagonal matrix elements first 
    for band in range(n_bands):
        for direc in range(3):
            diagonal_array[0,:,band,direc]=(vel[:,band,band,direc]*vel[:,band,band,direc]).real
            full_array[0,:,band,band,direc]=(vel[:,band,band,direc]*vel[:,band,band,direc]).real

    #load in the off-diago
    for band1 in range(n_bands):
        for band2 in range(n_bands):
            if band1==band2: #only want off-diagonal elements 
                continue
            for direc in range(6):
                if direc<3:
                    #full array stores 6 total directions, here xx/yy/zz
                    full_array[0,:,band1,band2,direc]=(vel[:,band1,band2,index1[direc]]*vel[:,band2,band1,index2[direc]]).real

                    # Berry-Curvature array only stores 3 directions xy/xz/yz
                    berry_array[0,:,band1,band2,direc]=-0.5*(vel[:,band1,band2,index3[direc]]*vel[:,band2,band1,index4[direc]]-vel[:,band2,band1,index3[direc]]*vel[:,band1,band2,index4[direc]]).imag

        
                else:
                    #load in directions xy/xz/yz
                    full_array[0,:,band1,band2,direc]=0.5*(vel[:,band1,band2,index1[direc]]*vel[:,band2,band1,index2[direc]]+vel[:,band2,band1,index1[direc]]*vel[:,band1,band2,index2[direc]]).real

    return energy_array,diagonal_array,full_array,berry_array

def load_hk(filename): #loads in the k-space hamiltonian from a given JSON file

    with open(filename,'r') as A:

      config = json.load(A)

      if config.get("d0")==None: #use direct Hamiltonian input format
            
            if config.get("Hk")==None:
                raise HkInputError("Hamiltonian is missing or incorrectly formatted") #raises an error if the wrong format is provided
            
            hk=config["Hk"]

            for i in range(len(hk)):
                  for j in range(len(hk)):
                      if type(hk[i][j])==str:#only need to parse the variable if it is a string
                            hk[i][j]=parse_expr(hk[i][j],local_dict=symbols_dict) #parses the strings into sympy variables

            hk=Matrix(hk) #turn it into the Sympy Matrix object- needed to be compatible with Matrix_Elements_Gen() 

      else: #use the Pauli matrix format
          
            ds=[config["d0"],config["d1"],config["d2"],config["d3"]]

            for i in range(4):
                if type(ds[i])==str: #only parse the variable if it is a string
                    ds[i]=parse_expr(ds[i], local_dict=symbols_dict)  #parses the Pauli matrix coeffient strings into sympy variables 

            d0=ds[0]
            d1=ds[1]
            d2=ds[2]
            d3=ds[3]
 
            hk=Hk(d0,d1,d2,d3) #builds the hamiltonian from 2x2 Pauli matrices 

      #reads in the reciprocal lattice vectors 

      if config.get("b1")==None or config.get("b2")==None or config.get("b3")==None:
           raise HkInputError("Reciprocal Lattice Vector is missing from input or incorrectly formatted") #raises an error if the wrong format provided 

      if type(config["b1"])==str: #if user specifies the recp lattice vectors as a complete string parse as one variable
          
          b1=np.array(parse_expr(config["b1"], local_dict=symbols_dict)) #parses the reciprocal lattice vector strings into sympy variables 
          b2=np.array(parse_expr(config["b2"], local_dict=symbols_dict))
          b3=np.array(parse_expr(config["b3"], local_dict=symbols_dict))

      else: #if user decides to specify individual components as strings then need to parse indivdual components
          
          bs=[config["b1"],config["b2"],config["b3"]]

          for i in range(3):
              for j in range(3):
                if type(bs[i][j])==str:
                    bs[i][j]=parse_expr(bs[i][j], local_dict=symbols_dict)

          b1=np.array(bs[0])
          b2=np.array(bs[1])
          b3=np.array(bs[2])

      if config.get("Name")==None:
          name=False
      else:
          name=config["Name"]

    return hk, b1, b2, b3, name
