import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*


from concaveDelaunayRefinement import*

from TessFunction import*
from functionInWork import*


import matplotlib
matplotlib.rcParams.update({'font.size': 14})

def compute_std(nb):
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    for ii in range(len(nb.num)):
        x.append(nb.pos[ii][0])
        y.append(nb.pos[ii][1])
        z.append(nb.pos[ii][2])
        vx.append(nb.vel[ii][0])
        vy.append(nb.vel[ii][1])
        vz.append(nb.vel[ii][2])
    std_x = np.std(x)
    std_y = np.std(y)
    std_z = np.std(z)
    std_vx = np.std(vx)
    std_vy = np.std(vy)
    std_vz = np.std(vz)
    return std_x,std_y,std_z,std_vx,std_vy,std_vz

def PosVel_to_w(nb):
    Nparticle = len(nb.pos)
    w = []

    for ii in range(Nparticle):
        w.append([nb.pos[ii][0],nb.pos[ii][1],nb.pos[ii][2],nb.vel[ii][0],nb.vel[ii][1],nb.vel[ii][2]])
    return w

def create_w_tList_cell(nb_tList,id_particle,tri0):
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][id_particle]:tri0.vertex_neighbor_vertices[0][id_particle+1]]

    nb_id = []
    for ii in id_particle_in_cell:
        nb_id.append(nb_tList[0].num[ii])

    w_tList = []
    for tt in range(len(nb_tList)):
        w_tList.append(PosVel_to_w(nb_tList[tt].selectp(nb_id)))
    return w_tList

def scale_nb(nb,std_List):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = std_List[0],std_List[1],std_List[2],std_List[3],std_List[4],std_List[5]
    nb.pos[:,0] = nb.pos[:,0]/std_x
    nb.pos[:,1] = nb.pos[:,1]/std_y
    nb.pos[:,2] = nb.pos[:,2]/std_z
    nb.vel[:,0] = nb.vel[:,0]/std_vx
    nb.vel[:,1] = nb.vel[:,1]/std_vy
    nb.vel[:,2] = nb.vel[:,2]/std_vz
    return nb

def scale_nb_tList(nb_tList):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = compute_std(nb_tList[0])
    for tt in range(len(nb_tList)):
        scale_nb(nb_tList[tt],[std_x,std_y,std_z,std_vx,std_vy,std_vz])
    return nb_tList



def create_w_tList_cell_with_scaling(nb_tList,id_particle,tri0):
    
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][id_particle]:tri0.vertex_neighbor_vertices[0][id_particle+1]]

    nb_id = []
    for ii in id_particle_in_cell:
        nb_id.append(nb_tList[0].num[ii])

    
    w_tList = []
    for tt in range(len(nb_tList)):
        w_tList.append(PosVel_to_w(nb_tList[tt].selectp(nb_id)))
    return w_tList

def find_center_3D(pointsList):
    return 1/4*(pointsList[0]+pointsList[1]+pointsList[2]+pointsList[3])

def find_center_6D(pointsList):
    return 1/7*(pointsList[0]+pointsList[1]+pointsList[2]+pointsList[3]+pointsList[4]+pointsList[5]+pointsList[6])

def find_center_points(index_concave_delaunay,tri):
    X_center = []
    for ii in index_concave_delaunay:
        center = find_center_6D(tri.points[tri.simplices[ii]])
        X_center.append([center[0],center[1],center[2],center[3],center[4],center[5]])
    return X_center

def find_center_points_tList(index_concave_delaunay_tList,tri_tList):
    Center_position_tList = []
    for tt in range(len(tri_tList)):
        Center_position_tList.append(find_center_points(index_concave_delaunay_tList[tt],tri_tList[tt]))
    return Center_position_tList


def draw1Frame_Configuration_Space(tt,points_tList,xLim):
    x = []
    y = []
    z = []
    for ii in range(len(points_tList[tt])):
        x.append(points_tList[tt][ii][0])
        y.append(points_tList[tt][ii][1])
        z.append(points_tList[tt][ii][2])

    fig, axs = plt.subplots(1,2,figsize=(15.0, 5.0))
    axs[0].scatter(x, y, marker='.',s=1)
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    #axs[0].plot(x,y,'or')
    axs[1].scatter(x,z, marker='.',s=1)
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    #axs[1].plot(vx,vy,'or')
    axs[1].set_xlabel('x position')
    axs[1].set_ylabel('z position')
    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGif_Configuration_Space(points_tList,name,Nframe=100,xLim=5):
    img0 = draw1Frame_Configuration_Space(0,points_tList,xLim)
    
    Nstep = len(points_tList)


    imgList = [img0]
    freq = Nstep//Nframe

    if freq<1:
        freq=1
    for ii in range(Nstep):
        if ii%freq==0:
            imgList.append(draw1Frame_Configuration_Space(ii,points_tList,xLim))
    img0.save(name, save_all=True, append_images=imgList)

def draw1Frame_Full_Configuration_Space(nb,xLim):

    fig, axs = plt.subplots(1,2,figsize=(15.0, 5.0))
    axs[0].scatter(nb.pos[:,0], nb.pos[:,1], marker='.',s=1)
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    #axs[0].plot(x,y,'or')
    axs[1].scatter(nb.pos[:,0],nb.pos[:,2], marker='.',s=1)
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    #axs[1].plot(vx,vy,'or')
    axs[1].set_xlabel('x position')
    axs[1].set_ylabel('z position')
    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGif_Full_Configuration_Space(nb_tList,name,Nframe=100,xLim=5):
    img0 = draw1Frame_Full_Configuration_Space(nb_tList[0],xLim)
    
    Nstep = len(nb_tList)


    imgList = [img0]
    freq = Nstep//Nframe

    if freq<1:
        freq=1
    for tt in range(Nstep):
        if tt%freq==0:
            imgList.append(draw1Frame_Full_Configuration_Space(nb_tList[tt],xLim))
    img0.save(name, save_all=True, append_images=imgList)

def axs_particle_traj(axs,nb_tList,part_id,xLim):
    nb = nb_tList[0]
    x = []
    y = []
    z = []
    for tt in range(len(nb_tList)):
        nb_tmp = nb_tList[tt].selectp([part_id])
        x.append(nb_tmp.pos[0][0])
        y.append(nb_tmp.pos[0][1])
        z.append(nb_tmp.pos[0][2])

    axs[0].scatter(nb.pos[:,0], nb.pos[:,1], marker='.',s=1)
    axs[0].plot(x,y,ls='solid',lw=2,color='red')
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')

    axs[1].scatter(nb.pos[:,0],nb.pos[:,2], marker='.',s=1)
    axs[1].plot(x,z,ls='solid',lw=2,color='red')
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    axs[1].set_xlabel('x position')
    axs[1].set_ylabel('z position')

def draw1Frame_Cell_Configuration_Space(nb,points,xLim):
    w = np.zeros((len(points),3))
    for ii in range(len(points)):
        w[ii,:] = np.array([points[ii][0],points[ii][1],points[ii][2]])

    fig, axs = plt.subplots(1,2,figsize=(15, 5))
    axs[0].scatter(nb.pos[:,0], nb.pos[:,1], marker='.',s=1)
    axs[0].scatter(w[:,0], w[:,1], marker='.',s=1)
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    #axs[0].plot(x,y,'or')
    axs[1].scatter(nb.pos[:,0],nb.pos[:,2], marker='.',s=1)
    axs[1].scatter(w[:,0],w[:,2], marker='.',s=1)
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    #axs[1].plot(vx,vy,'or')
    axs[1].set_xlabel('x position')
    axs[1].set_ylabel('z position')
    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGif_Cell_Configuration_Space(nb_tList,points_in_cell_tList,name,Nframe=100,xLim=5):
    img0 = draw1Frame_Cell_Configuration_Space(nb_tList[0],points_in_cell_tList[0],xLim)
    
    Nstep = len(nb_tList)

    imgList = [img0]
    freq = Nstep//Nframe

    if freq<1:
        freq=1
    for tt in range(Nstep):
        if tt%freq==0:
            imgList.append(draw1Frame_Cell_Configuration_Space(nb_tList[tt],points_in_cell_tList[tt],xLim))
    img0.save(name, save_all=True, append_images=imgList)

def draw1Frame_Cell_and_center_Configuration_Space(nb,points,points_cell_center,xLim):
    w = np.zeros((len(points),3))
    w_center = np.zeros((len(points_cell_center),3))
    for ii in range(len(points)):
        w[ii,:] = np.array([points[ii][0],points[ii][1],points[ii][2]])
    for jj in range(len(points_cell_center)):
        w_center[jj,:] = np.array([points_cell_center[jj][0],points_cell_center[jj][1],points_cell_center[jj][2]])

    fig, axs = plt.subplots(1,2,figsize=(15, 5))
    #axs[0].scatter(nb.pos[:,0], nb.pos[:,1], marker='.',s=1)
    axs[0].scatter(w[:,0], w[:,1], marker='.',s=1)
    axs[0].scatter(w_center[:,0], w_center[:,1], marker='.',s=1)
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x position')
    axs[0].set_ylabel('y position')
    #axs[1].scatter(nb.pos[:,0],nb.pos[:,2], marker='.',s=1)
    axs[1].scatter(w[:,0],w[:,2], marker='.',s=1)
    axs[1].scatter(w_center[:,0], w_center[:,2], marker='.',s=1)
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    axs[1].set_xlabel('x position')
    axs[1].set_ylabel('z position')
    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGif_Cell_and_center_Configuration_Space(nb_tList,points_in_cell_tList,points_cell_center_tList,name,Nframe=100,xLim=5):
    img0 = draw1Frame_Cell_and_center_Configuration_Space(nb_tList[0],points_in_cell_tList[0],points_cell_center_tList[0],xLim)
    
    Nstep = len(nb_tList)

    imgList = [img0]
    freq = Nstep//Nframe

    if freq<1:
        freq=1
    for tt in range(Nstep):
        if tt%freq==0:
            imgList.append(draw1Frame_Cell_and_center_Configuration_Space(nb_tList[tt],points_in_cell_tList[tt],points_cell_center_tList[tt],xLim))
    img0.save(name, save_all=True, append_images=imgList)


def draw1Frame_Cell_and_center_velocity_Space(nb,points,points_cell_center,xLim):
    w = np.zeros((len(points),3))
    w_center = np.zeros((len(points_cell_center),3))
    for ii in range(len(points)):
        w[ii,:] = np.array([points[ii][3],points[ii][4],points[ii][5]])
    for jj in range(len(points_cell_center)):
        w_center[jj,:] = np.array([points_cell_center[jj][3],points_cell_center[jj][4],points_cell_center[jj][5]])

    fig, axs = plt.subplots(1,2,figsize=(15, 5))
    #axs[0].scatter(nb.pos[:,0], nb.pos[:,1], marker='.',s=1)
    axs[0].scatter(w[:,0], w[:,1], marker='.',s=1)
    axs[0].scatter(w_center[:,0], w_center[:,1], marker='.',s=1)
    axs[0].set_xlim([-xLim,xLim])
    axs[0].set_ylim([-xLim,xLim])
    axs[0].set_xlabel('x velocity')
    axs[0].set_ylabel('y velocity')
    #axs[1].scatter(nb.pos[:,0],nb.pos[:,2], marker='.',s=1)
    axs[1].scatter(w[:,0],w[:,2], marker='.',s=1)
    axs[1].scatter(w_center[:,0], w_center[:,2], marker='.',s=1)
    axs[1].set_xlim([-xLim,xLim])
    axs[1].set_ylim([-xLim,xLim])
    axs[1].set_xlabel('x velocity')
    axs[1].set_ylabel('z velocity')
    
    img = fig2img(fig)
    plt.close()
    return  img

def makeGif_Cell_and_center_velocity_Space(nb_tList,points_in_cell_tList,points_cell_center_tList,name,Nframe=100,xLim=5):
    img0 = draw1Frame_Cell_and_center_velocity_Space(nb_tList[0],points_in_cell_tList[0],points_cell_center_tList[0],xLim)
    
    Nstep = len(nb_tList)

    imgList = [img0]
    freq = Nstep//Nframe

    if freq<1:
        freq=1
    for tt in range(Nstep):
        if tt%freq==0:
            imgList.append(draw1Frame_Cell_and_center_velocity_Space(nb_tList[tt],points_in_cell_tList[tt],points_cell_center_tList[tt],xLim))
    img0.save(name, save_all=True, append_images=imgList)