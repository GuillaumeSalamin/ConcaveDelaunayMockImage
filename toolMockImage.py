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

def scale_nb(nb,std_List):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = std_List[0],std_List[1],std_List[2],std_List[3],std_List[4],std_List[5]
    nb.pos[:,0] = nb.pos[:,0]/std_x
    nb.pos[:,1] = nb.pos[:,1]/std_y
    nb.pos[:,2] = nb.pos[:,2]/std_z
    nb.vel[:,0] = nb.vel[:,0]/std_vx
    nb.vel[:,1] = nb.vel[:,1]/std_vy
    nb.vel[:,2] = nb.vel[:,2]/std_vz
    return nb,[std_x,std_y,std_z,std_vx,std_vy,std_vz]

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