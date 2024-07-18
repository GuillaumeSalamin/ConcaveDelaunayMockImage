import numpy as np
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

def download_simulation(name,trans,rot):


    nb = Nbody(name,ftype='swift')

    nb.translate(trans)
    nb.rotate(angle=rot)
    return nb

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
    return [std_x,std_y,std_z,std_vx,std_vy,std_vz]

def scale_nb(nb,std_List):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = std_List[0],std_List[1],std_List[2],std_List[3],std_List[4],std_List[5]
    nb.pos[:,0] = nb.pos[:,0]/std_x
    nb.pos[:,1] = nb.pos[:,1]/std_y
    nb.pos[:,2] = nb.pos[:,2]/std_z
    nb.vel[:,0] = nb.vel[:,0]/std_vx
    nb.vel[:,1] = nb.vel[:,1]/std_vy
    nb.vel[:,2] = nb.vel[:,2]/std_vz
    return nb

def PosVel_to_w(nb):
    Nparticle = len(nb.pos)
    w = []

    for ii in range(Nparticle):
        w.append([nb.pos[ii][0],nb.pos[ii][1],nb.pos[ii][2],nb.vel[ii][0],nb.vel[ii][1],nb.vel[ii][2]])
    return w

def find_neighbors(tri):
    ans = []
    for ii in range(len(tri.points)):
        tmp = tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][ii]:tri.vertex_neighbor_vertices[0][ii+1]]
        ans.append(tmp)
    return ans

#=============================================================================================================
#               fonction pour HaloPartialImage
#=============================================================================================================
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*

from concaveDelaunayRefinement import*
from toolMockImage import*

import pickle

import matplotlib
matplotlib.rcParams.update({'font.size': 14})

def triangleArea(x1,x2,x3):
    T =1/2*np.abs((x1[0]-x3[0])*(x2[1]-x1[1])-(x1[0]-x2[0])*(x3[1]-x1[1]))
    return T

def delaunay_volume(tri):
    points=tri.points
    vol = []
    for ii in range(len(tri.simplices)):
        x1 = points[tri.simplices[ii]][0]
        x2 = points[tri.simplices[ii]][1]
        x3 = points[tri.simplices[ii]][2]
        vol.append(triangleArea(x1,x2,x3))
    return vol

def distance6D(x1,x2):
    return np.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2+(x1[3]-x2[3])**2+(x1[4]-x2[4])**2+(x1[5]-x2[5])**2)


def simplex_projection2d(points,ax1,ax2):
    pro1 = [points[0][ax1],points[0][ax2]]
    pro2 = [points[1][ax1],points[1][ax2]]
    pro3 = [points[2][ax1],points[2][ax2]]
    pro4 = [points[3][ax1],points[3][ax2]]
    pro5 = [points[4][ax1],points[4][ax2]]
    pro6 = [points[5][ax1],points[5][ax2]]
    pro7 = [points[6][ax1],points[6][ax2]]
    return [pro1,pro2,pro3,pro4,pro5,pro6,pro7]

def create_w_tList_cell(nb_tList,neighbors,id_particle,Rmax):
    
    nb0 = nb_tList[0]
    id_particle_in_cell = neighbors[id_particle]
    w0 =   [nb0.pos[id_particle][0],nb0.pos[id_particle][1],nb0.pos[id_particle][2],nb0.vel[id_particle][0],nb0.vel[id_particle][1],nb0.vel[id_particle][2]]
    nb_id = []
    for ii in id_particle_in_cell:
        w =   [nb0.pos[ii][0],nb0.pos[ii][1],nb0.pos[ii][2],nb0.vel[ii][0],nb0.vel[ii][1],nb0.vel[ii][2]]
        dist_tmp = distance6D(w0,w)
        if dist_tmp<Rmax:
            nb_id.append(nb_tList[0].num[ii])

    
    w_tList = []
    for tt in range(len(nb_tList)):
        w_tList.append(PosVel_to_w(nb_tList[tt].selectp(nb_id)))
    return w_tList

def Scale_triangle(vertex):
    for ii in range(len(vertex)):
        vertex[ii][0] = vertex[ii][0]/100
        tmp = vertex[ii][1]/100
        vertex[ii][1] = vertex[ii][2]/100
        vertex[ii][2] = tmp
        vertex[ii][3] = vertex[ii][3]/20
        vertex[ii][4] = vertex[ii][4]/20
        vertex[ii][5] = vertex[ii][5]/20
    return vertex

def Scale_triangle_2d(vertex):
    for ii in range(len(vertex)):
        vertex[ii][0] = vertex[ii][0]/100
        vertex[ii][1] = vertex[ii][1]/100
    return vertex

def projection_6d_to_2d(tri,id_list,ax1,ax2):
    tri2d_list = []
    alpha_list = []
    #V6d_list = []
    #V2d_list = []
    jj=0
    n=len(id_list)
    for ii in id_list:
        #print(f'simplex {jj} of {n}', end='\r')
        jj=jj+1
        #if ii in id_list:
        V = Tetrahedron6DVolume(tri.points[tri.simplices[ii]])
        w2d = simplex_projection2d(tri.points[tri.simplices[ii]],ax1,ax2)
        tri2d = Delaunay(w2d)
        V2d_ListSimplices = delaunay_volume(tri2d)
        V2d = np.sum(V2d_ListSimplices)
        for ii in range(len(V2d_ListSimplices)):
            #V2d_ListSimplices[ii] = V2d_ListSimplices[ii]*V/V2d
            V2d_ListSimplices[ii] = V/V2d
        alpha = V2d_ListSimplices
        tri2d_list.append(tri2d)
        alpha_list.append(alpha)
        #V6d_list.append(V)
        #V2d_list.append(V2d)
    return tri2d_list,alpha_list#V6d_list,V2d_list

def VertexList2d(tri,id_list,ax1,ax2):
    triangle_List = []
    alpha = []
    tri2d_list,alpha_list = projection_6d_to_2d(tri,id_list,ax1,ax2)
    for ii in range(len(tri2d_list)):
        for jj in range(len(tri2d_list[ii].simplices)):
            vertices = tri2d_list[ii].points[tri2d_list[ii].simplices[jj]]
            triangle_List = [*triangle_List,*Scale_triangle_2d(vertices)]
            alpha.append(alpha_list[ii][jj])
    return triangle_List,alpha