import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*

import time

import sys
sys.path.append('/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis')

from concaveDelaunayRefinement import*
from toolMockImage import*
from TessFunction import*
from functionInWork import*

import pickle

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


def create_w_tList_cell_v8(nb_tList,neighbors,id_particle):
    
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = neighbors[id_particle]

    nb_id = []
    for ii in id_particle_in_cell:
        nb_id.append(nb_tList[0].num[ii])

    
    w_tList = []
    for tt in range(len(nb_tList)):
        w_tList.append(PosVel_to_w(nb_tList[tt].selectp(nb_id)))
    return w_tList

def create_w_tList_cell_v9(nb_tList,neighbors,id_particle,Rmax):
    
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

def create_w_tList_cell_v7(nb_tList,id_particle,tri0,Rmax,Nmax):
    nb0 = nb_tList[0]
    nb_id = []
    dist = []
    w0 = [nb0.pos[id_particle][0],nb0.pos[id_particle][1],nb0.pos[id_particle][2],nb0.vel[id_particle][0],nb0.vel[id_particle][1],nb0.vel[id_particle][2]]
    for ii in range(len(nb0.num)):
        w =   [nb0.pos[ii][0],nb0.pos[ii][1],nb0.pos[ii][2],nb0.vel[ii][0],nb0.vel[ii][1],nb0.vel[ii][2]]
        dist_tmp = distance6D(w0,w)
        if dist_tmp<Rmax:
            nb_id.append(nb_tList[0].num[ii])
            dist.append(dist_tmp)
    if len(dist)>Nmax+1:
        dist_lim = (sorted(dist))[Nmax]
        nb_id = np.array(nb_id)
        dist = np.array(dist)
        nb_id=nb_id[dist<dist_lim]
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


print('scaling data')
nb_scale_tList = scale_nb_tList(download_simulation_v2())

file_path = '/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonSaveVariable/neighbors_v2.pickle'
with open(file_path, 'rb') as file:
    # Deserialize and retrieve the variable from the file
    loaded_data = pickle.load(file)

print("The variable 'data' has been loaded successfully.")

neighbors = loaded_data

file_path = '/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/DelaunayTri/tri0.pickle'
with open(file_path, 'rb') as file:
    # Deserialize and retrieve the variable from the file
    loaded_data = pickle.load(file)

print("The variable 'data' has been loaded successfully.")

tri0 = loaded_data

filePara = 0

parameter_a_list = []
for ii in range(10):
    parameter_a_list.append(filePara*10+ii)

for parameter_a in parameter_a_list:
    cell_id_list = []
    triangle_List = []
    opacity_list = []
    data = []
    for ii in range(len(nb_scale_tList[0].num)):
        if ii%100==parameter_a:
            cell_id_list.append(ii)

    n = len(cell_id_list)
    Nmin = 20
    jj=0

    
    #triangle_List2 = []
    
    for p_id in cell_id_list:
        jj=jj+1
        print(f'DEBUG :: iteration {jj} of {len(cell_id_list)}')
        print(f'DEBBUG :: cell creation')
        #w_scale_tList = create_w_tList_cell_v8(nb_scale_tList,neighbors,p_id)
        w_scale_tList = create_w_tList_cell_v9(nb_scale_tList,neighbors,p_id,Rmax=50)
        if len(w_scale_tList[0])>Nmin:
            print('DEBUG : Concave Delaunay Algorithm')
            index_conserved_simplex_tList,tri_tList,flag_tList = Volume_evolution_Liouville_condition_concave_delaunay_refinements_test(w_scale_tList,error_max=5,verbos=False)
            flag = flag_tList[-1]

            index_List = index_conserved_simplex_tList[-1]
            tri = tri_tList[-1]
            V_cons = delaunay_volume_6D(tri_tList[0])

            if flag==0:
                tmp,opa = VertexList2d(tri,index_List,ax1=0,ax2=2)
                opa = opa/V_cons
                triangle_List = [*triangle_List,*tmp]
                opacity_list = [*opacity_list,*opa]
                #for ii in range(len(tri.simplices)):
                    #if ii in index_List:
                        #vertices = tri.points[tri.simplices[ii]]
                        #triangle_List2 = [*triangle_List2,*Scale_triangle(vertices)]

    print(f'number of triangle :: {len(triangle_List)}')
    name = f'Full_Image_{parameter_a}'

    
    for ii in range(len(triangle_List)):
        data.append([triangle_List[ii][0],triangle_List[ii][1],opacity_list[ii//3]])


    file_path_data = f'/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/totalHalo/triangle3/{name}_data.pickle'
    with open(file_path_data, 'wb') as file:
        # Serialize and write the variable to the file
        pickle.dump(data, file)

    print(f'The variable "data" has been saved successfully. \n File name :: {file_path_data}')