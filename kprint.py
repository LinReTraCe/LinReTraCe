import numpy as np
import h5py 
import matplotlib.pyplot as plt
import argparse
from argparse import RawTextHelpFormatter
#from hk.hk_matrix_elements import *

def parse_args(args=None): #parses the arguments from command line
  parser = argparse.ArgumentParser(
    description='Argument parser for the Hk or TB pre-processing for LRTC', \
    formatter_class=RawTextHelpFormatter)
  
  parser.add_argument('file', type=str, help='LRTC file') 
  parser.add_argument('obj', type=str, help='Quantity to be plotted')
  parser.add_argument('direction', type=str, help='Direction to be plotted')
  parser.add_argument('path',nargs='?',default=None, type=str, help='path through BZ to be plotted')
  parser.add_argument('--BZ',required=False,action="store_true",help="Plot the BZ")
 
  return parser.parse_args(args)


def BZ_plot(b1,b2,spec_points,k_array): #plots the BZ using reciprocal lattice vectors and plots the path specified through the BZ 
    fig, ax = plt.subplots()
    ax.quiver(0,0,b1[0],b1[1],angles='xy',scale_units='xy',scale=1) #draws the b1 reciprocal lattice vector
    ax.quiver(0,0,b2[0],b2[1],angles='xy',scale_units='xy',scale=1) #draws the b2 reciprocal lattice vector
    ax.quiver(b1[0],b1[1],b2[0],b2[1],angles='xy',scale_units='xy',scale=1,color="gray") #draws the b1 edge of BZ
    ax.quiver(b2[0],b2[1],b1[0],b1[1],angles='xy',scale_units='xy',scale=1,color="gray") #draws the b2 edge of BZ

    for i in spec_points: #plots any special points in the BZ
        coord=special_point(i) #extract the coordiante of special points provided 
        ax.plot([coord[0]*b1[0]+coord[1]*b2[0]],[coord[0]*b1[1]+coord[1]*b2[1]],'o',color='red') 
        if i=="G" or i=="g":
            ax.text(-0.3*np.sign(b1[0]),0,r"$\Gamma$") #if Gamma point displays using the Greek symbol
        elif i=="P":
            ax.text((coord[0]*b1[0]+coord[1]*b2[0])+0.2*np.sign(b1[0]),coord[0]*b1[1]+coord[1]*b2[1],"K'")
        else:
            ax.text((coord[0]*b1[0]+coord[1]*b2[0])+0.2*np.sign(b1[0]),coord[0]*b1[1]+coord[1]*b2[1],i) #plots the text of the corresponding special point
    
    x=[]
    y=[]
    for i in range(len(k_array)): #append the x and y components of the provided list of k-points for plotting 
        x.append(k_array[i][0])
        y.append(k_array[i][1])

    ax.plot(x,y,'x',color="orange") #plot the provided path of k-points in BZ
             
    if b1[1]==0: #if the reciprocal lattice vectors are orthogonal then change the x and y-lims are handeled 
        ax.text(0.9*b1[0],0.9*b1[1],r"$\boldsymbol{b}_1$")
        ax.text(0.9*b2[0],0.9*b2[1],r"$\boldsymbol{b}_2$")
        ax.set_xlim(-0.1*np.sign(b1[0]),(2*np.pi+0.1)*np.sign(b1[0]))
        ax.set_ylim(-0.1*np.sign(b1[0]),(2*np.pi+0.1)*np.sign(b1[0]))
    else: #x and y lims for non-orthogonal reciprocal lattice vectors 
        ax.text(0.75*b1[0],b1[1],r"$\boldsymbol{b}_1$")
        ax.text(0.75*b2[0],b2[1],r"$\boldsymbol{b}_2$")
        ax.set_xlim(-0.1*np.sign(b1[0]),2.1*b1[0])
        ax.set_ylim(-1.1*b1[1],1.1*b1[1])
    ax.set_aspect("equal") #make sure plotting the true scale 
    ax.axis("off")
    plt.show()

def index_convert(ikx,jky,nkx,nky): #converts a coordinate in (kx,ky) space into the stored index value which can be used to access the correpsonding k-point from hdf5 file

    i=round(ikx*nkx)
    j=round(jky*nky)

    return i*nky+j

def line_eq_finder(X1,Y1,X2,Y2): #finds the equation of a line from two points in (kx,ky) space 

    M=np.array([[X1,1],[X2,1]])
    v=np.array([[Y1],[Y2]])

    #special cases for infinite or no gradient
    if X1==X2:
        return None, Y1
    if Y1==Y2:
        return 0, X1
    
    #find gradient and intercpet 
    else:
        Minv=np.linalg.inv(M)
        sol=np.dot(Minv,v)

        m=sol[0][0]
        c=sol[1][0]

    return m,c


def line_eq_output(m,c,xi): 

    return m*xi+c

def line_eq_output_y(m,c,yi):

    return (yi-c)/m

def line_gen(X1,Y1,X2,Y2,nkx,nky,endpoint=False): #generates a list of indices corresponding to a the line of points between (X1,Y1) and (X2,Y2)

    m,c=line_eq_finder(X1,Y1,X2,Y2)

    index_list=[] #stores the indices to be returned 
    
    if m==None: #special case for a line of constant kx

        Yend=max([Y1,Y2])
        Ystart=min([Y1,Y2])
        dy=1/nky #spacing between ky points (as fractions of b2)

        #generate a list of ky points (as fractions of b2)
        if endpoint==False:
            iy_list=np.arange(Ystart,Yend,dy)
        if endpoint==True: #if the endpoint is included need extra point
            n = int(round((Yend - Ystart) / dy)) + 1
            iy_list = np.linspace(Ystart, Yend, n)

        if Y1>Y2: #line always has to go from 1 to 2
            if endpoint==False:
                iy_list=np.arange(Ystart+dy,Yend+dy,dy)
            if endpoint==True: #if the endpoint is included need extra point
                n=int(round((Yend-Ystart)/dy))+1
                iy_list=np.linspace(Ystart,Yend,n)

            iy_list=np.flip(iy_list) #have to flip the array so the path goes from 1 to 2

        for iy in iy_list:
            index_list.append(index_convert(X1,iy,nkx,nky)) #use the ky points and constant kx point to generate a list of indices 
            
        return index_list
        
    if m==0: #special case for a line of constant ky 

        Xend=max([X1,X2])
        Xstart=min([X1,X2])
        dx=1/nkx #spacing between kx points (as fractions of b1)

        #generate a list of kx points (as fractions of b1)
        if endpoint==False:
            ix_list=np.arange(Xstart,Xend,dx)
        if endpoint==True: #if the endpoint is included need extra point
            n = int(round((Xend - Xstart) / dx)) + 1
            ix_list = np.linspace(Xstart, Xend, n)

        if X1>X2:
            if endpoint==False:
                ix_list=np.arange(Xstart+dx,Xend+dx,dx)
            if endpoint==True: #if the endpoint is included need extra point
                n=int(round((Xend-Xstart)/dx))+1
                ix_list=np.linspace(Xstart,Xend,n)

            ix_list=np.flip(ix_list) #have to flip the array so the path goes from 1 to 2

        for ix in ix_list:
            index_list.append(index_convert(ix,Y1,nkx,nky)) #use the kx points and constant ky point to generate a list of indices 

        return index_list

    if m>=1 or m<=-1: #for gradients of this ranges use a selection of kx-points to generate the indices 

        Xend=max([X1,X2])
        Xstart=min([X1,X2])
        dx=1/nkx #spacing between kx points (as fractions of b1)

        if endpoint==False:
            ix_list=np.arange(Xstart,Xend,dx)
            #n = int(np.floor((Xend - Xstart) / dx))
            #ix_list = np.linspace(Xstart, Xstart + (n - 1)*dx, n)
        if endpoint==True: #if the endpoint is included need extra point
            n = int(round((Xend - Xstart) / dx)) + 1
            ix_list = np.linspace(Xstart, Xend, n)

        if X1>X2:
            if endpoint==False:
                #n = int(np.floor((Xend - Xstart) / dx))
                #ix_list = np.linspace(Xstart + dx, Xstart + n*dx, n)            
                ix_list=np.arange(Xstart+dx,Xend+dx,dx)
            if endpoint==True: #if the endpoint is included need extra point
                n=int(round((Xend-Xstart)/dx))+1
                ix_list=np.linspace(Xstart,Xend,n)

            ix_list=np.flip(ix_list) #have to flip the array so the path goes from 1 to 2

        for ix in ix_list:
            index_list.append(index_convert(ix,line_eq_output(m,c,ix),nkx,nky)) #use the kx points and the equation of the line to generate a list of indices 
        return index_list

    if m<1 and m>-1: #for gradients of this ranges use a selection of ky-points to generate the indices 

        Yend=max([Y1,Y2])
        Ystart=min([Y1,Y2])
        dy=1/nky #spacing between ky points (as fractions of b2)

        if endpoint==False:
            iy_list=np.arange(Ystart,Yend,dy)
            #n = int(np.floor((Yend - Ystart) / dy))
            #iy_list = np.linspace(Ystart, Ystart + (n - 1)*dy, n)
        if endpoint==True: #if the endpoint is included need extra point
            n = int(round((Yend - Ystart) / dy)) + 1
            iy_list = np.linspace(Ystart, Yend, n)

        if Y1>Y2:
            if endpoint==False:
                iy_list=np.arange(Ystart+dy,Yend+dy,dy)
                #n = int(np.floor((Yend - Ystart) / dy))
                #iy_list = np.linspace(Ystart + dy, Ystart + n*dy, n) 
            if endpoint==True: #if the endpoint is included need extra point
                n=int(round((Yend-Ystart)/dy))+1
                iy_list=np.linspace(Ystart,Yend,n)

            iy_list=np.flip(iy_list) #have to flip the array so the path goes from 1 to 2

        for iy in iy_list:
            index_list.append(index_convert(line_eq_output_y(m,c,iy),iy,nkx,nky)) #use the ky points and the equation of the line to generate a list of indices 

        return index_list

   

def special_point(point): #converts the name of special points into coordinates in (kx,ky) space 

    if point=="G":
        return np.array([0,0])
    if point=="K":
        return np.array([2/3,1/3])
    if point=="P":
        return np.array([1/3,2/3])
    if point=="M":
        return np.array([1/2,1/2])
    if point=="g":
        return np.array([1,1])
    if point=="T":
        return np.array([0,1])
    if point=="B":
        return np.array([1,0])
    if point=="C":
        return np.array([0.5,1])
    if point=="D":
        return np.array([1,0.5])
    if point=="E":
        return np.array([0.5,0])
    if point=="F":
        return np.array([0,0.5])
    else:
        raise KeyError("Special Point not recognized")

def path_reader(path,nkx,nky):

    path_list=list(path) #seperates the string into individual components

    if path_list[-1]=="G": #if the path ends in the Gamma point, change this to the (1,1) Gamma point so path goes across BZ 
        path_list[-1]="g"

    coord_list=[special_point(i) for i in path_list] #converts string into coordinates in the (kx,ky) space for each special point
   
    indices=[] #stores the indices along the entire path
    path_lengths=[] #stores the length of each component of the path- used for locating special points 

    for i in range(len(path_list)-1): #a path made of n special points has n-1 components to it

        if i<len(path_list)-2:
            line_arr=line_gen(coord_list[i][0],coord_list[i][1],coord_list[i+1][0],coord_list[i+1][1],nkx,nky) #new component of the path- DOESN'T contain the endpoint
            path_lengths.append(len(line_gen(coord_list[i][0],coord_list[i][1],coord_list[i+1][0],coord_list[i+1][1],nkx,nky))) #lentgh of the current component of the path

        else: #for the final component of the path we need to include the endpoint to display correctly
            line_arr=line_gen(coord_list[i][0],coord_list[i][1],coord_list[i+1][0],coord_list[i+1][1],nkx,nky,True) #final component of the path INCLUDING the endpoint
            path_lengths.append(len(line_gen(coord_list[i][0],coord_list[i][1],coord_list[i+1][0],coord_list[i+1][1],nkx,nky,True))) #length of final path component

        for j in line_arr:
            indices.append(j) #add the new component of the path to the indices list 

    return indices,path_lengths

def obj_reader(obj): #reads an object as a string and outputs into numerical format

    obj_list=list(obj)

    if obj_list[0]=="E" or obj_list[0]=="e": #if energy is used as obj

        if obj_list[-1]=="y": #if no band choice specificed then plot all energy bands, Band1=None
            return ["Energy",None,None]
        else:
            band=int(obj_list[-1])-1 #if band choice is specified then Band1=band

            return ["Energy",band,None]
    
    if obj_list[0]=="B" or obj_list[0]=="b":

        if obj_list[-1]=="v" or obj_list[-1]=="e" or obj_list[-1]=="C" or obj_list[-1]=="c": #if no band choice specificed then plot all berry curv, Band1=None
            return ["BerryCurv",None,None]
        else:
            band1=int(obj_list[-2])-1
            band2=int(obj_list[-1])-1
            return ["BerryCurv",band1,band2] #if band choice is specified then specify choice numerically
    
    
    if obj_list[0]=="M" or obj_list[0]=="m":

        if obj_list[-1]=="M" or obj_list[-1]=="m": #if no band choice specificed then plot all matrix elements, Band1=None
            return ["M",None,None]
        else:

            band1=int(obj_list[-2])-1
            band2=int(obj_list[-1])-1

            return ["M",band1,band2]#if band choice is specified then specify choice numerically
    
    
def direction_reader(direc): #takes in a direction as a string and outputs as a number

    if direc=="xx":
        return [0,None] #first num is general direction
                        #second num used just for Berry-Curv which has different structure
    if direc=="yy":
        return [1,None]
    if direc=="zz":
        return [2,None]
    if direc=="xy":
        return [3,0]
    if direc=="xz":
        return [4,1]
    if direc=="yz":
        return [5,2]
 

def K_space_Plot(filename,obj,direction,path,BZ):


    with h5py.File(filename, 'r') as f:
        n_bands=int(f[".bands/energyBandMax"][()]) #number of bands of system
        nkp=int(f[".kmesh/nkp"][()]) #number of k-points sampled
        nkx=f[".kmesh/nkx"][()] #number of kx points
        nky=f[".kmesh/nky"][()] #number of ky points
        get_berry=f.get('kPoint/0000000001/BerryCurv') #checks if BerrCurv is present, if it is not then get_berry=None
        bs=f[".unitcell/kvec"][()] #loads in the reciprocal lattice vectors

        irreducible=f[".kmesh/irreducible"][()]
    
    if irreducible==True: #plotting code only works with REDUCIBLE k-meshes
        raise KeyError("K-mesh must be irreducible, please use lprint")
    
    indices,spec_point_loc=path_reader(path,nkx,nky) #indices=list of the index of the k-points sampled by the 'path'
                                                     #spec_point_loc= array containing size of each part of path [length path1, length path2...]  
    path_names=list(path)

    obj_name,Band1,Band2=obj_reader(obj) #reads in the object input, reads into the name="Energy"/"M"/"BerryCurv" and the bands indexs chosen to be plotted
                                         # if Band1==None then plot ALL of the chosen object
    k_array=[] #sampled k-point array
    
    with h5py.File(filename, 'r') as f: #append the sampled k-points to array
        for i in indices:
            if i>nkp-1: #if index is too big then reject it
                continue
            dx=f['.kmesh/points'][i][0]
            dy=f['.kmesh/points'][i][1]
            k_array.append(dx*bs[0]+dy*bs[1])
            
    delta_x=1/len(k_array) #distance between sampled points that we will use to plot on the 'x' axis
    
    direc_num=direction_reader(direction)[0] #converts a string input "xx"/"yy" etc to numerical input 0/1
    berry_direc=direction_reader(direction)[1] #as Berry curv only has 3 components need to convert differently
                                               #e.g "xy" = 0 for Berry

    xs=[] #stores the 'x' coordinate we will use to plot 
    xs_spec=[] #stores the 'x' coordinate of the special-points, this will be used for marking the path on the plot
    to_remove=[] #stores any indices which are too large and need to be removed from indices[]
    
    spec_indices=[0] #index corresponding to special points- the first index in the list will ALWAYS be the first special point
    current=0
    for i in range(len(spec_point_loc)):
        spec_indices.append(spec_point_loc[i]+current) #special points located at the start of each new compoenent of the path
        current+=spec_point_loc[i] 
    
    for ik in range(len(indices)):
        if indices[ik]>nkp-1:
            to_remove.append(indices[ik]) #save any incorrect indices to be removed
            continue
        if ik in spec_indices:
            xs_spec.append(ik*delta_x) #save the 'x' location of the special points
        xs.append(ik*delta_x) #save the 'x' location for all of the points in the path

    xs_spec.append(xs[-1]) #the firnal index will ALWAYS be the final special point in path
    
    for rem_index in to_remove:
        indices.remove(rem_index) #remove the incorrect indices 


    #here ALL of the corresponding object is plotted 
    if Band1==None: 

        #define the objects as arrays
        energy_hk=np.zeros((len(xs),n_bands))
        M_hk=np.zeros((len(xs),n_bands,n_bands))
        berry_hk=np.zeros((len(xs),n_bands,n_bands))


        #load in the objects from chosen file into the arrays        
        with h5py.File(filename, 'r') as f:
            for ik in range(len(indices)):
                for band1 in range(n_bands):
                    for band2 in range(n_bands):
                        if band1==band2:
                            energy_hk[ik,band1]=f['energies'][indices[ik]][band1]
                            M_hk[ik,band1,band1]=f['kPoint/{:010}/moments'.format(indices[ik]+1)][band1][band1][direc_num]
                        if band1!=band2:
                            M_hk[ik,band1,band2]=f['kPoint/{:010}/moments'.format(indices[ik]+1)][band1][band2][direc_num]

                            if berry_direc==None or get_berry==None: #if BerryCurv doesn't exist in the file OR if xx/yy/zz chosen then don't load the Berry curv
                                continue
                            else:
                                berry_hk[ik,band1,band2]=f['kPoint/{:010}/BerryCurv'.format(indices[ik]+1)][band1][band2][berry_direc]

        #plot all the energy bands
        if obj_name=="Energy":
            ymaxs=[] #stores max energy of all bands
            ymins=[] #stores min energy of all bands 
            for iband in range(n_bands): 
                plt.plot(xs,energy_hk[:,iband],label=f"Band{iband+1}") #plot each band in turn
                ymaxs.append(max(energy_hk[:,iband])) #add max energy
                ymins.append(min(energy_hk[:,iband])) #add min energy
            ymin=min(ymins) #pick out the smallest possible energy
            ymax=max(ymaxs) #pick out the biggest possible energy

            if ymin==0 and ymax==0: #if provided object is totally zero then set some limits for sake of plotting
                ymin=-0.5
                ymax=0.5

            #plot the special points as vertical lines
            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")     
                elif path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])  #labels of the special points
            plt.ylabel(r"$\epsilon$")

            ax = plt.gca() #aquires the most recent axis
            ax.get_xaxis().set_visible(False) #hides the x-axis
            plt.legend()
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array) #plots the BZ with path displayed
            plt.show()

        #plot the intra and inter band matrix elements
        if obj_name=="M":
            ymaxs=[] #stores max energy of all bands
            ymins=[] #stores min energy of all bands 
            for iband1 in range(n_bands):
                for iband2 in range(n_bands):
                    plt.plot(xs,M_hk[:,iband1,iband2],label=f"M{iband1+1}{iband2+1}{direction}") #plot all matrix elements in turn 
                ymaxs.append(max(M_hk[:,iband1,iband2])) 
                ymins.append(min(M_hk[:,iband1,iband2]))
            ymin=min(ymins)
            ymax=max(ymaxs)

            if ymin==0 and ymax==0:
                ymin=-0.5
                ymax=0.5

            #plot the special points as vertical lines
            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")    
                elif path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])  
            plt.ylabel(r"M")

            ax = plt.gca()#aquires the most recent axis
            ax.get_xaxis().set_visible(False)#hides the x-axis
            plt.legend()
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array)
            plt.show()

        #plot the 'berry curvature' anti-sym matrix elements 
        if obj_name=="BerryCurv":
            ymaxs=[]
            ymins=[]
            for iband1 in range(n_bands):
                for iband2 in range(n_bands):
                    if iband1==iband2: #only want to plot the 'off-diagonal' elements as if iband1=iband2 the Omega=0 
                        continue
                    else:
                        plt.plot(xs,berry_hk[:,iband1,iband2],label=rf"$\Omega${iband1+1}{iband2+1}{direction}")
                        ymaxs.append(max(berry_hk[:,iband1,iband2]))
                        ymins.append(min(berry_hk[:,iband1,iband2]))
            ymin=min(ymins)
            ymax=max(ymaxs)

            if ymin==0 and ymax==0:
                ymin=-0.5
                ymax=0.5

            #plot the special points as vertical lines
            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")
                    continue     
                if path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])  
            plt.ylabel(r"$\Omega$")

            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
            plt.legend()
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array)
            plt.show()
        

    else: #if Band1!=None then we plot a specific choice of object 

        if obj_name=="Energy":
            energy_hk=[]
            with h5py.File(filename, 'r') as f:
                for ik in range(len(indices)):
                    energy_hk.append(f['energies'][indices[ik]][Band1]) #plots energy of band Band1, e.g Energy1
            ymin=min(energy_hk)
            ymax=max(energy_hk)

            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")     
                if path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])       
            plt.plot(xs,energy_hk)

            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
            plt.ylabel(f"E{Band1+1}")
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array)
            plt.show()

        if obj_name=="M":
            M_hk=[]
            with h5py.File(filename, 'r') as f:
                for ik in range(len(indices)):
                    M_hk.append(f['kPoint/{:010}/moments'.format(indices[ik]+1)][Band1][Band2][direc_num]) #plots matrix element of band Band1 Band2, e.g M12
            ymin=min(M_hk)
            ymax=max(M_hk)

            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")     
                if path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])  
            plt.plot(xs,M_hk)

            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
            plt.ylabel(f"M{Band1+1}{Band2+1}{direction}")
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array)
            plt.show()

        if obj_name=="BerryCurv":
            berry_hk=[]
            with h5py.File(filename, 'r') as f:
                for i in range(len(indices)):
                    berry_hk.append(f['kPoint/{:010}/BerryCurv'.format(indices[i]+1)][Band1][Band2][berry_direc]) #plots berry curve element of band Band1 Band2, e.g Omega12
            ymin=min(berry_hk)
            ymax=max(berry_hk)

            for i in range(len(xs_spec)):
                plt.vlines(xs_spec[i],ymin,ymax,color="gray",linestyle="dashed")
                if path_names[i]=="G" or path_names[i]=="g":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,r"$\Gamma$")     
                if path_names[i]=="P":
                    plt.text(xs_spec[i]+0.01,ymin*1.08,"K'")
                else:
                    plt.text(xs_spec[i]+0.01,ymin*1.08,path_names[i])  
            plt.plot(xs,berry_hk)

            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
            plt.ylabel(rf"$\Omega${Band1+1}{Band2+1}{direction}")
            if BZ==True:
                BZ_plot(bs[0],bs[1],path_names,k_array)
            plt.show()



args = parse_args()

if args.path==None: #allows for energy to be plotted without the need for a direction to be specified
    args.path=args.direction #if no direction is asked then the path will actually be stored in the direction arg so we swap this
    args.direction="xx" #set a generic direction, for energy this doesn't effect the plotting
    name=list(args.obj)
    if name[0]== "E" or name[0]== "e": #check that it was energy that was called
        pass
    else: #if a different obj was called the direction is REQUIRED so we raise an error
        raise KeyError("No direction provided")
    
K_space_Plot(args.file,args.obj,args.direction,args.path,args.BZ)
