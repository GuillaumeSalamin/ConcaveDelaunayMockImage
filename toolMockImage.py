import numpy as np
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

"""
Ce fichier contient tout les fonction nécessaire à la création du mock Image avec la méthode du concace Delaunay refinement.
Les méthode propre au Concave Delaunay refinement sont dans un autre fichier (concaveDelaunayRefinement.py). 
Ce fichier contient tout les fonctions nécessaire à la création de Mock Image uniquement. Se sont principalement des function 'utilitaire'.
"""

def download_simulation(name,trans,rot):
    """
    def download_simulation(name,trans,rot):
    Download un fichier hdf5 contenant une snapshot de Nbody simulation.
    La function utilise pNbody et la simulation doit provenir du code SWIFT.
    - name : string du path au ficheir hdf5
    - trans : vecteur [x,y,z] de la translation de la simulation avec pNbody, nb = nb.translate(trans) (le code a besoin que les particules soit centrer en (0,0,0))
    - rota : angle de rotation de pNbody, nb = nb.rotate(angle=rot)
    """


    nb = Nbody(name,ftype='swift')

    nb = nb.translate(trans)
    nb = nb.rotate(angle=rot)
    return nb

def compute_std(nb):
    """
    def compute_std(nb):
    fonctino calculant la standar deviation de x,y,z,vx,vy,vz. Les positions et vitesses
    sont argument et un Nbody object de pNbody
    retourne un vecteur [std_x,std_y,std_z,std_vx,std_vy,std_vz]
    """
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
    """
    def scale_nb(nb,std_List):
    Prend un Nbody object et fait un sorte que la dispersion dans les 6 dimension de l'espace de phase soit 1.
    Scaling choisie pour construire des volumes dans l'espace de phase.
    - nb : un Nbody object de pNbody
    - std_List : une liste [std_x,std_y,std_z,std_vx,std_vy,std_vz]. La fonction compute_std(nb) permet de l'obtenir.
    """
    std_x,std_y,std_z,std_vx,std_vy,std_vz = std_List[0],std_List[1],std_List[2],std_List[3],std_List[4],std_List[5]
    nb.pos[:,0] = nb.pos[:,0]/std_x
    nb.pos[:,1] = nb.pos[:,1]/std_y
    nb.pos[:,2] = nb.pos[:,2]/std_z
    nb.vel[:,0] = nb.vel[:,0]/std_vx
    nb.vel[:,1] = nb.vel[:,1]/std_vy
    nb.vel[:,2] = nb.vel[:,2]/std_vz
    return nb

def PosVel_to_w(nb):
    """
    def PosVel_to_w(nb):
    prend un Nbody object de pNbody en entrée.
    retourne un array de array6 [w0,w1,w2,...] avec w_i un array de taille 6 [x,y,z,vx,vy,vz]
    """
    Nparticle = len(nb.pos)
    w = []

    for ii in range(Nparticle):
        w.append([nb.pos[ii][0],nb.pos[ii][1],nb.pos[ii][2],nb.vel[ii][0],nb.vel[ii][1],nb.vel[ii][2]])
    return w

def find_neighbors(tri):
    """
    def find_neighbors(tri):
    Prend entrée un triangulation de Delaunay (scipy.spatial) et revoie un array telle que :
    ans[i] est la liste des points voisine de i dans la triangulation
    i étant le points au coordonné tri.points[i]
    """
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
    """
    def triangleArea(x1,x2,x3):
    
    Calcule l'aire d'un triangle (2d) donc les sommets sont aux positions, x1,x2,x3 (x_i un array/list 2d)
    """
    T =1/2*np.abs((x1[0]-x3[0])*(x2[1]-x1[1])-(x1[0]-x2[0])*(x3[1]-x1[1]))
    return T

def delaunay_volume(tri):
    """
    def delaunay_volume(tri):
    
    Compute the list of volume of each triangle of a Delaunay triangulation. This function is for 2d triangulation
    Input :
        - tri : a Delaunay triangulation (scipy)
    return :
        - vol : a list of float
    vol[i] is the volume of simplices tri.simplices[i]
    """
    points=tri.points
    vol = []
    for ii in range(len(tri.simplices)):
        x1 = points[tri.simplices[ii]][0]
        x2 = points[tri.simplices[ii]][1]
        x3 = points[tri.simplices[ii]][2]
        vol.append(triangleArea(x1,x2,x3))
    return vol

def distance6D(x1,x2):
    """
    Compute the distance between x1 and x2, two points in 6 dimension
    """
    return np.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2+(x1[3]-x2[3])**2+(x1[4]-x2[4])**2+(x1[5]-x2[5])**2)


def simplex_projection2d(points,ax1,ax2):
    """
    def simplex_projection2d(points,ax1,ax2):
    For a simplex in 6D (with 7 vertex) computes the projection of the point in a 2d plane of axe1,axe2
    ax=0 -> x
    ax=1 -> y
    ax=2 -> z
    ax=3 -> vx
    ax=4 -> vy
    ax=5 -> vz
    return a vector of 7 two-dimensional coordinates. ([[x1,y1],[x2,y2],..,[x7,y7]])
    """
    pro1 = [points[0][ax1],points[0][ax2]]
    pro2 = [points[1][ax1],points[1][ax2]]
    pro3 = [points[2][ax1],points[2][ax2]]
    pro4 = [points[3][ax1],points[3][ax2]]
    pro5 = [points[4][ax1],points[4][ax2]]
    pro6 = [points[5][ax1],points[5][ax2]]
    pro7 = [points[6][ax1],points[6][ax2]]
    return [pro1,pro2,pro3,pro4,pro5,pro6,pro7]

def create_w_tList_cell(nb_tList,neighbors,id_particle,Rmax):
    """
    def create_w_tList_cell(nb_tList,neighbors,id_particle,Rmax):
    
    Input :
        - nb_tList : A list of Nbody object. for exemple a List of frame ordered with time of the simulation
        - neighbors : List of neighbors for the triangulation. neighbors[i] is neighbors of vertex i. Computed with script CreateTri.py and stored in a pickle file
        - id_particle : particle central of the cell
        - Rmax : maximum distance from the central particle to one of his neighbors for the neighbors to be in the cell. Ideally Rmax = infty
    Return :
        - w_tList : List of List of position of each particle in each frame
     w_tList = [[List of position in frame 0],[List of position in frame 1],...]
     List of position in frame i = [...,[6 dimension array of phase space position of jth particle in the cell in frame i],..]

     Compute the list of position of the particle in cell  id_particle for each frame given.
     The frame given should be at least (and most of the time) two frame. One initial frame, another frame whom we want to make an image.
    """
    
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

def Scale_triangle_2d(vertex,std_List,ax1,ax2,halfBoxSize):
    """
    def Scale_triangle_2d(vertex,std_List,ax1,ax2,halfBoxSize):
    
    this function scale the triangle's vertex to be projected. First std_list is use again to scale each coordinate to his true value. Then halfBoxSize is used such as all particle between [-halfBoxSize,halfBoxSize] is between [-1,1] to be prjected by pyOpenGl.
    The image will include particle with coordinate in [-halfBoxSize,halfBoxSize].
    input:
        - vertex : list of vertex of the projected triangle. List of 2 dimensional array
        - std_List : List [std_x,std_y,std_z,std_vx,std_vy,std_vz] compute by compute_std(nb)
        - ax1 : axe of x coordinate of the image
        - ax2 : axe of y coordinate of the image
        - halfBoxSize : the half of the boxSize of the image (in Unit of the simulation)
    return:
        - List of vertex scale by halfBoxSize from the original value in the simulation
    """
    for ii in range(len(vertex)):
        vertex[ii][0] = vertex[ii][0]/halfBoxSize*std_List[ax1]
        vertex[ii][1] = vertex[ii][1]/halfBoxSize*std_List[ax2]
    return vertex

def projection_6d_to_2d(tri,id_list,ax1,ax2):
    """
    def projection_6d_to_2d(tri,id_list,ax1,ax2):
    
        fonction to compute the projection

    input: 
        - tri : Delaunay triangulation of the image frame
        - id_list : id_list of the simplices in the concave refinement obtain with concaveDelaunayRefinement.py
        - ax1 : axe of x coordinate of the image
        - ax2 : axe of y coordinate of the image
    return:
        - tri2d_list: List of Delaunay triangulation in two dimension. The triangle of the triangulation are the triangle to display
        - for each triangle, a parameter alpha (matrix of float) representing it's intensity (abitrary unit). (it's a matrix alpha[trianglulation][triangle of the triangulation])
    """
    tri2d_list = []
    alpha_list = []
    jj=0
    n=len(id_list)
    for ii in id_list:
        jj=jj+1
        V = Tetrahedron6DVolume(tri.points[tri.simplices[ii]])
        w2d = simplex_projection2d(tri.points[tri.simplices[ii]],ax1,ax2)
        tri2d = Delaunay(w2d)
        V2d_ListSimplices = delaunay_volume(tri2d)
        V2d = np.sum(V2d_ListSimplices)
        for ii in range(len(V2d_ListSimplices)):
            V2d_ListSimplices[ii] = V/V2d
        alpha = V2d_ListSimplices
        tri2d_list.append(tri2d)
        alpha_list.append(alpha)
    return tri2d_list,alpha_list

def VertexList2d(tri,std_List,id_list,ax1,ax2,halfBoxSize):
    """
    input: 
        - tri : Delaunay triangulation of the image frame
        - id_list : id_list of the simplices in the concave refinement obtain with concaveDelaunayRefinement.py
        - ax1 : axe of x coordinate of the image
        - ax2 : axe of y coordinate of the image
        - halfBoxSize: the image will contain particle in [-halfBoxSize,halfBoxSize]
    return:
        - triangle_List: List of triangle (three two-dimensonal array) to display
        - for each triangle, a parameter alpha (float) representing it's intensity (abitrary unit)
    """
    triangle_List = []
    alpha = []
    tri2d_list,alpha_list = projection_6d_to_2d(tri,id_list,ax1,ax2)
    for ii in range(len(tri2d_list)):
        for jj in range(len(tri2d_list[ii].simplices)):
            vertices = tri2d_list[ii].points[tri2d_list[ii].simplices[jj]]
            triangle_List = [*triangle_List,*Scale_triangle_2d(vertices,std_List,ax1,ax2,halfBoxSize)]
            alpha.append(alpha_list[ii][jj])
    return triangle_List,alpha