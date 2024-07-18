import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from pNbody import*
from TessFunction import*


def Tetrahedron6DVolume(vectorList): # je ne suis pas sur dutout que la formule utiliser soit juste
    vectorList= np.array(vectorList)
    M=np.zeros((6,6))
    v0 = vectorList[0]
    for jj in range(6):
        M[jj,:]= vectorList[jj+1]-v0
    #print(f'DEBUG :: M=\n{M}')
    return np.abs(1/5040*np.linalg.det(M))

def Simplex4DVolume(vectorList): # je ne suis pas sur dutout que la formule utiliser soit juste
    vectorList= np.array(vectorList)
    M=np.zeros((4,4))
    v0 = vectorList[0]
    for jj in range(4):
        M[jj,:]= vectorList[jj+1]-v0
    #print(f'DEBUG :: M=\n{M}')
    return np.abs(1/24*np.linalg.det(M))

def Tetrahedron_vol(vectorList):# je ne suis pas sur dutout que la formule utiliser soit juste
    vectorList= np.array(vectorList)
    M=np.zeros((3,3))
    v0 = vectorList[0]
    for jj in range(3):
        M[jj,:]= vectorList[jj+1]-v0
    #print(f'DEBUG :: M=\n{M}')
    return np.abs(1/6*np.linalg.det(M))


def Delaunay_density_3d(tri):
    points = tri.points
    dens = []
    NpointFraq = points_per_triangle(tri)
    for ii in range(len(tri.simplices)):
        x1 = points[tri.simplices[ii]][0]
        x2 = points[tri.simplices[ii]][1]
        x3 = points[tri.simplices[ii]][2]
        x4 = points[tri.simplices[ii]][3]
        dens.append(NpointFraq[ii]*1/Tetrahedron_vol([x1,x2,x3,x4]))
    return dens

def voronoi_volumes(vor):
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
        indices = vor.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            #vol[i] = ConvexHull(vor.vertices[indices]).volume
            #print(f'DEBUG :: vor.vertices[indices]={vor.vertices[indices]}')
            #print(np.shape(vor.vertices[indices]))
            vol[i] = volume6D(vor.vertices[indices])
    return vol

def volume6D(points):
    tri = Delaunay(points)
    return np.sum(Delaunay_volume_6d(tri))

def Delaunay_volume_6d(tri):
    points=tri.points
    vol = []
    for ii in range(len(tri.simplices)):
        x1 = points[tri.simplices[ii]][0]
        x2 = points[tri.simplices[ii]][1]
        x3 = points[tri.simplices[ii]][2]
        x4 = points[tri.simplices[ii]][3]
        x5 = points[tri.simplices[ii]][4]
        x6 = points[tri.simplices[ii]][5]
        x7 = points[tri.simplices[ii]][6]
        vol.append(Tetrahedron6DVolume([x1,x2,x3,x4,x5,x6,x7]))
    return vol



##############################################################################################
# Dernière fonction de SimAnalyse
##############################################################################################


def find_next_ID(jj,nbNow,nbNext):
    xweight = 10
    r = []
    for ii in range(len(nbNow.pos)):
        r.append(np.sqrt(xweight*(nbNow.pos[jj][0]-nbNext.pos[ii][0])**2+xweight*(nbNow.pos[jj][1]-nbNext.pos[ii][1])**2+xweight*(nbNow.pos[jj][2]-nbNext.pos[ii][2])**2+(nbNow.vel[jj][0]-nbNext.vel[ii][0])**2+(nbNow.vel[jj][1]-nbNext.vel[ii][1])**2+(nbNow.vel[jj][2]-nbNext.vel[ii][2])**2))
    return np.argmin(r)

def track_ID(ii,nb_sim):
    ID = [ii]
    for tt in range(len(nb_sim)-1):
        ID.append(find_next_ID(ID[-1],nb_sim[tt],nb_sim[tt+1]))
    return ID



def create_time_list(Nstep,dt):
    return np.linspace(0,Nstep,Nstep)*dt

def create_particle_list(particleIndexList,nb_sim,norm_pos=1,norm_vel =1):
    NtimeStep = len(nb_sim)
    Nparticle = len(particleIndexList)
    pointsList = []
    
    #IDList_List = [track_ID(particleIndexList[0],nb_sim),track_ID(particleIndexList[1],nb_sim),track_ID(particleIndexList[2],nb_sim),track_ID(particleIndexList[3],nb_sim),track_ID(particleIndexList[4],nb_sim),track_ID(particleIndexList[5],nb_sim),track_ID(particleIndexList[6],nb_sim)]
    IDList_List = []
    for ii in range(len(particleIndexList)):
        IDList_List.append(track_ID(particleIndexList[ii],nb_sim))


    for ii in range(NtimeStep):
        points = []
        for particleIndex in IDList_List:

            addPoint = [nb_sim[ii].pos[particleIndex[ii]][0]/norm_pos,nb_sim[ii].pos[particleIndex[ii]][1]/norm_pos,nb_sim[ii].pos[particleIndex[ii]][2]/norm_pos,nb_sim[ii].vel[particleIndex[ii]][0]/norm_vel,nb_sim[ii].vel[particleIndex[ii]][1]/norm_vel,nb_sim[ii].vel[particleIndex[ii]][2]/norm_vel]
            points.append(addPoint)
        pointsList.append(points)
    return pointsList

def create_tri_list(pointsList):
    triList= []
    for ii in range(len(pointsList)):
        triList.append(Delaunay(pointsList[ii]))
    return triList

def Delaunay_volume_6d(tri):
    vol = []
    for ii in range(len(tri.simplices)):
        vol.append(Tetrahedron6DVolume(tri.points))
    return vol

def volume_list_simulation(triList):
    volListSim = []
    for tt in range(len(triList)):
        volListSim.append(Delaunay_volume_6d(triList[tt]))
    return volListSim

def volume_list_one_simplices_simulation(triList):
    volListSim = []
    for tt in range(len(triList)):
        volListSim.append(Tetrahedron6DVolume(triList[tt].points))
    return volListSim



def Axs_Analyse_Evolve_Volume(nb_sim,simplexIndex,tri,axs):
    
    x_scale = np.std(nb_sim[0].pos[:,0])
    v_scale = np.std(nb_sim[0].vel[:,0])

    time = create_time_list(len(nb_sim),1)

    particleIndexList = tri.simplices[simplexIndex]

    pointsList = create_particle_list(particleIndexList,nb_sim,norm_pos=x_scale,norm_vel=v_scale)
    triList = create_tri_list(pointsList)
    volList = volume_list_one_simplices_simulation(triList)

    axs[0].plot(time, volList)
    #axs.set_ylim([0,max(volList)])
    axs[0].set_xlabel('time [# step]')
    axs[0].set_ylabel('simplex volume')
    axs[1].plot(time, volList-np.mean(volList))
    axs[1].set_xlabel('time [# step]')
    axs[1].set_ylabel(r'V - <V>$_t$')
    return time,pointsList,volList


def draw1FrameHalo(tt,nb_sim,particleIndex):

    x = [nb_sim[tt].pos[particleIndex][0]]
    y = [nb_sim[tt].pos[particleIndex][1]]
    vx = [nb_sim[tt].vel[particleIndex][0]]
    vy = [nb_sim[tt].vel[particleIndex][1]]
    """x = [pos_sim[tt][particleIndex][0]]
    y = [pos_sim[tt][particleIndex][1]]
    vx = [nb_sim[tt].vel[particleIndex][0]]
    vy = [nb_sim[tt].vel[particleIndex][1]]"""


    xlim = 1.20
    vlim = 2

    fig, axs = plt.subplots(1,2,figsize=(7.5, 2.5))
    axs[0].plot(nb_sim[tt].pos[:,0], nb_sim[tt].pos[:,1], '.',markersize='1')
    #axs[0].plot(pos_sim[tt][:,0], pos_sim[tt][:,1], '.',markersize='1')
    axs[0].set_xlim([-xlim,xlim])
    axs[0].set_ylim([-xlim,xlim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    axs[0].set_title('position space')
    axs[0].plot(x,y,'or')
    axs[1].plot(nb_sim[tt].vel[:,0], nb_sim[tt].vel[:,1], '.',markersize='1')
    axs[1].set_xlim([-vlim,vlim])
    axs[1].set_ylim([-vlim,vlim])
    axs[1].plot(vx,vy,'or')
    axs[1].set_xlabel('x velocity')
    axs[1].set_ylabel('y velocity')
    axs[1].set_title('velocity space')
    
    img = fig2img(fig)
    plt.close()
    return  img

def draw1FrameHalo5Particle(tt,nb_sim,particleIndexList):

    xlim = 1.20
    vlim = 2

    fig, axs = plt.subplots(1,2,figsize=(7.5, 2.5))
    axs[0].plot(nb_sim[tt].pos[:,0], nb_sim[tt].pos[:,1], '.',markersize='1')
    #axs[0].plot(pos_sim[tt][:,0], pos_sim[tt][:,1], '.',markersize='1')
    axs[0].set_xlim([-xlim,xlim])
    axs[0].set_ylim([-xlim,xlim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    axs[0].set_title('position space')
    
    axs[1].plot(nb_sim[tt].vel[:,0], nb_sim[tt].vel[:,1], '.',markersize='1')
    axs[1].set_xlim([-vlim,vlim])
    axs[1].set_ylim([-vlim,vlim])
    
    axs[1].set_xlabel('x velocity')
    axs[1].set_ylabel('y velocity')
    axs[1].set_title('velocity space')

    for ParticleIndex in particleIndexList:
        x = [nb_sim[tt].pos[ParticleIndex[tt]][0]]
        y = [nb_sim[tt].pos[ParticleIndex[tt]][1]]
        vx = [nb_sim[tt].vel[ParticleIndex[tt]][0]]
        vy = [nb_sim[tt].vel[ParticleIndex[tt]][1]]
        axs[1].plot(vx,vy,'o')
        axs[0].plot(x,y,'o')
    """x = [pos_sim[tt][particleIndex][0]]
    y = [pos_sim[tt][particleIndex][1]]
    vx = [nb_sim[tt].vel[particleIndex][0]]
    vy = [nb_sim[tt].vel[particleIndex][1]]"""

    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGifHallo():
    return 0

def makeGifHallo5Particle():
    return 0


#========================================================================================================================================
#                                   Nouvelle façon de calculer les volumes, 2 dimension
#========================================================================================================================================
def is_inside_2d(points_list,point,dx,dy):
    for ii in range(len(points_list)):
        if np.abs(point[0]-points_list[ii][0])<dx and np.abs(point[1]-points_list[ii][1])<dy:
            #print(f'DEBUG :: np.abs(point[0]-points_list[ii][0])={np.abs(point[0]-points_list[ii][0])}')
            #print(f'DEBUG :: np.abs(point[1]-points_list[ii][1]={np.abs(point[1]-points_list[ii][1])}')
            return True
    return False

def volumeDiscretSpace_2D(points,N):
    x = []
    y = []
    for ii in range(len(points)):
        x.append(points[ii][0])
        y.append(points[ii][1])
    
    xMax = max(x)
    xMin = min(x)
    yMax = max(y)
    yMin = min(y)
    xMat = np.linspace(xMin-0.25*np.abs(xMin),xMax+0.25*np.abs(xMax),N)
    yMat = np.linspace(yMin-0.25*np.abs(yMin),yMax+0.25*np.abs(yMax),N)
    dx = xMat[1]-xMat[0]
    dy = yMat[1]-yMat[0]
    #VolMat = np.zeros((N,N))
    vol_tot = 0
    for ii in range(N):
        for jj in range(N):
            xx,yy = xMat[ii],yMat[jj]
            if is_inside_2d(points,[xx,yy],dx,dy):
                #VolMat[ii,jj]=1
                vol_tot = vol_tot+dx*dy
            #else:
                #VolMat[ii,jj]=0

    return vol_tot

def volumeDiscretSpace_2D_test(points,N):
    
    plt.figure(figsize=(8, 8))
    plt.axis('equal')
    x = []
    y = []
    for ii in range(len(points)):
        x.append(points[ii][0])
        y.append(points[ii][1])
    
    xMax = max(x)
    xMin = min(x)
    yMax = max(y)
    yMin = min(y)
    xMat = np.linspace(xMin-0.25*np.abs(xMin),xMax+0.25*np.abs(xMax),N)
    yMat = np.linspace(yMin-0.25*np.abs(yMin),yMax+0.25*np.abs(yMax),N)
    dx = xMat[1]-xMat[0]
    dy = yMat[1]-yMat[0]
    VolMat = np.zeros((N,N))
    for ii in range(N):
        for jj in range(N):
            xx,yy = xMat[ii],yMat[jj]
            if is_inside_2d(points,[xx,yy],dx,dy):
                VolMat[ii,jj]=1
                plt.fill([xx,xx,xx+dx,xx+dx],[yy,yy+dy,yy+dy,yy],'r')
            else:
                VolMat[ii,jj]=0
                plt.fill([xx,xx,xx+dx,xx+dx],[yy,yy+dy,yy+dy,yy],'b')
    plt.scatter(x,y,c='k')
    plt.show()
    print(f'DEBUG :: dx={dx}')
    print(f'DEBUG :: dx={dy}')
    return VolMat


def volume_evolve_time_2d(points_tList,N):
    Nstep = len(points_tList)
    vol_tList = np.zeros(Nstep)
    for tt in range(Nstep):
        print(f'frame={tt} of {Nstep}', end='\r', flush=True)
        vol_tList[tt] = volumeDiscretSpace_2D(points_tList[tt],N)
    return vol_tList

def convert_points_list(xParticle,vParticle):
    points = []
    for tt in range(len(xParticle[0])):
        tmp = []
        for ii in range(len(xParticle)):
            tmp.append([xParticle[ii][tt],vParticle[ii][tt]])
        points.append(tmp)
    return points

def draw1FrameManyParticleHO1D_volume(tt,points_tList,xLim,vLim):
    fig, axs = plt.subplots(figsize=(5, 5))
    x=[]
    v=[]
    for ii in range(len(points_tList[tt])):
        x.append(points_tList[tt][ii][0])
        v.append(points_tList[tt][ii][1])


    points = points_tList[tt]
    N=20
    y = v
    
    xMax = max(x)
    xMin = min(x)
    yMax = max(y)
    yMin = min(y)
    xMat = np.linspace(xMin-0.25*np.abs(xMin),xMax+0.25*np.abs(xMax),N)
    yMat = np.linspace(yMin-0.25*np.abs(yMin),yMax+0.25*np.abs(yMax),N)
    dx = xMat[1]-xMat[0]
    dy = yMat[1]-yMat[0]
    VolMat = np.zeros((N,N))
    for ii in range(N):
        for jj in range(N):
            xx,yy = xMat[ii],yMat[jj]
            if is_inside_2d(points,[xx,yy],dx,dy):
                VolMat[ii,jj]=1
                axs.fill([xx,xx,xx+dx,xx+dx],[yy,yy+dy,yy+dy,yy],'r')
            else:
                VolMat[ii,jj]=0
                axs.fill([xx,xx,xx+dx,xx+dx],[yy,yy+dy,yy+dy,yy],'b')



    axs.scatter(x, v,s=10)
    axs.set_xlim([-xLim,xLim])
    axs.set_ylim([-vLim,vLim])
    axs.set_xlabel('x position')
    axs.set_ylabel('x velocity')
    axs.set_title('Harmonic oscillator phase space trajectory')

    img = fig2img(fig)
    plt.close()
    return  img

def create_gif_HO1D_Nparticle_volume(gif_Name,points_tList,Nframe = 100, xLim=1,vLim=1):

    print(f'DEBUG :: create gif')
    img0 =draw1FrameManyParticleHO1D_volume(0,points_tList,xLim,vLim)
    print(f'DEBUG :: image 0')
    imgList = [img0]

    freq = len(points_tList)//Nframe

    if freq<1:
        freq=1

    #count = 0
    for ii in range(len(points_tList)):
        if ii%freq==0:
            print(f'frame={ii} of {len(points_tList)}', end='\r', flush=True)
            imgList.append(draw1FrameManyParticleHO1D_volume(ii,points_tList,xLim,vLim))

    img0.save(gif_Name, save_all=True, append_images=imgList)


#========================================================================================================================================
#                                   Nouvelle façon de calculer les volumes, 4 dimension
#========================================================================================================================================


#========================================================================================================================================
#                                   Autre
#========================================================================================================================================
def change_list_to_nparray(SimList):
    Xvec = np.zeros(np.shape(SimList))

    for ii in range(len(SimList)):
        for jj in range(3):
            Xvec[ii,jj] = SimList[ii][jj]
    return Xvec