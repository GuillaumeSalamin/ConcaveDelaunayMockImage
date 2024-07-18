# ce fichier ne fait que contenir des fonctions copier d'ailleurs, la version à jour et travailler ce trouve allieurs

##############################################################################################
# Concave Delaunay
##############################################################################################


import numpy as np
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

##############################################################################################
# Fonction trouver sur github pour calculer les cercles circoncrit
##############################################################################################

def hypercircumsphere(points, return_R=False):
    # First, ensure that we are working with a numpy array
    if not isinstance(points, np.ndarray):
        points = np.array(points)
    # Also ensure that we are working with a floating point type
    if not isinstance(points[0][0], np.floating):
        points = np.array(points, dtype=float)
        
    # Check that for N dimensions we have N+1 points
    shape = points.shape
    if len(shape) != 2 or shape[0] != shape[1] + 1:
        raise ValueError("Exactly N+1 points needed for N dimensions")

    # Compute the  directions of the vectors between the points.
    # Implicitly these are the normals of the planes that are normal to them
    directions = points[1:] - points[0]

    # If the points all lie in some lower-dimensional subspace then there is
    # no solution to this problem
    if np.linalg.det(directions) == 0:
        raise ValueError(f"Points all lie in a subspace smaller than {shape[1]} dimensions")
        
    # Compute the mid-points along these connections
    midpoints = (points[1:] + points[0])/2

    # A plane can be defined by its normal and a constant value for the dot
    # product of the points and that nomrmal. We want the planes that go through
    # the midpoint between each pair of points, so the constant can bbe found
    # by just evaluating at the midpoint
    plane_k = (directions * midpoints).sum(axis=1)
    
    # The center of the circumsphere is where the planes intersect, so we solve
    # for the simultanious set of plane equations
    answer = np.linalg.solve(directions, plane_k)

    if return_R:
        return (answer, ((points[0] - answer) ** 2).sum() ** 0.5)
    else:
        return answer    

##############################################################################################
# 6D Algorithm
##############################################################################################

def Tetrahedron6DVolume(vectorList): 
    """
    Compute the 6 dimensional volume formed by the simplex of the 7 points given in input 
    """
    vectorList= np.array(vectorList)
    M=np.zeros((6,6))
    v0 = vectorList[0]
    for jj in range(6):
        M[jj,:]= vectorList[jj+1]-v0
    #print(f'DEBUG :: M=\n{M}')
    return np.abs(1/5040*np.linalg.det(M))

def delaunay_volume_6D(tri):
    """"
    Comput the volume as the sum of all simplex in the triangulation. This is equivalent to the volume of the convex hull

    Input :
    - Delaunay tesselation (scipy) : The delaunay tesselation of the set of point 
    return
    - float : the volume
    """
    vol = 0
    for ii in range(len(tri.simplices)):
        vol += Tetrahedron6DVolume(tri.points[tri.simplices[ii]])
    return vol

def find_radius_list(tri):
    """
    Take in entry a delaunay tesselation (scipy) and return a list as the simplex tri.simplices[ii] as a circumcircle radius of r_list[ii]. r_list the outuput of this function.
    """
    r_list = []
    for ii in range(len(tri.simplices)):
        _, Radius = hypercircumsphere(tri.points[tri.simplices[ii]], return_R=True)
        r_list.append(Radius)
    return r_list

def find_volume_list(tri):
    """
    Take in entry a delaunay tesselation (scipy) and return a list of volume of each simplex.
    """
    vol_list = []
    for ii in range(len(tri.simplices)):
        vol_list.append(Tetrahedron6DVolume(tri.points[tri.simplices[ii]]))
    return vol_list

def concave_delaunay_volume_6D_alpha_volume_known(alpha,volume_list,radius_list):
    """"
    def concave_delaunay_volume_6D_alpha_volume_known(alpha,volume_list,radius_list):
    Compute the volume of the shape formed by a set of point and a concave parameter alpha.
    Take the delaunay tesselation and compute the volume of the circumcircle with radius less than 1/alpha. (Concave refinement)
    Note : if alpha = 0, this is just the volume of the convex hull or the volume of the normal delaunay tesselation.
    Here to speed up calculation, the radius of circumcircle and volume of each simplices are input, so already compute.
    Algo
    vol=0
    for each simplex in the delaunay triangulation
        R<- simplex circumcircle
        if alpha<1/R
            vol+=simplex volume
    return vol

    Input :
    - Delaunay tesselation (scipy) : The delaunay tesselation of the set of point 
    - float : alpha, the concave parameter
    - list of float : Radius_list, the list of radius for each simplices (to speed up computation)
    - list of float : volume_list, the list of volume for each simplices (to speed up computation)
    return
    - float : the volume of the shape of concave parameter alpha
    """
    vol = 0
    for ii in range(len(radius_list)):
        Radius = radius_list[ii]
        if alpha<1/Radius:
            vol += volume_list[ii]
    return vol

def find_alpha_concave_delaunay_simplices_index(tri,alpha,circumcircle_radius_list):
    """
    Compute the index of the simplex inside a concave delaunay tesselation of parameter concave alpha
    """
    index_list = []
    for ii in range(len(tri.simplices)):
        Radius = circumcircle_radius_list[ii]
        if alpha<1/Radius:
            index_list.append(ii)
    return index_list






def Find_simplex_list_target_volume_concave_delaunay_refinement(target_volume,tri,error_max):
    """
    def Find_simplex_list_target_volume_concave_delaunay_refinement(target_volume,tri,error_max):
    Find the shape with the targeted volume by using concave refinement of the delaunay tesselation take in input.
    Input :
        - float : A volume (the target)
        - Delaunay tesselation (scipy) : the delaunay tesselation who will be concavely refined
        - float : error_max, pourcentage (in %) of maximum difference between volum output and target volume
    output :
        - List of int : a list of simplex index forming the final shape with volume=target_volume
    """
    flag = 0

    circumcircle_radius_list = find_radius_list(tri)
    vol_list = find_volume_list(tri)
    alpha = 1/np.median(circumcircle_radius_list)
    alpha0 = np.copy(alpha)
    vol = 0
    maxIter = 200
    NIter = 0

    dalpha=0.5 # step size of alpha in fraction of it's original value
    sign_change = 0 # to avoid loop
    while (np.abs(target_volume-vol)/(target_volume)*100>error_max) and NIter<maxIter:
        NIter +=1
        vol = concave_delaunay_volume_6D_alpha_volume_known(alpha,vol_list,circumcircle_radius_list)
        if vol<target_volume:
            alpha = alpha-dalpha*alpha0
            if sign_change==1:
                dalpha=dalpha/2
            sign_change = -1
        if vol>target_volume:
            alpha = alpha+dalpha*alpha0
            if sign_change==-1:
                dalpha=dalpha/2
            sign_change = 1
        if NIter%50==0 and NIter!=0:
            dalpha=np.random.random_sample()*alpha0
    if NIter >= maxIter or (np.abs(target_volume-vol)/(target_volume)*100>error_max):
        flag = 1
    index_List = find_alpha_concave_delaunay_simplices_index(tri,alpha,circumcircle_radius_list)
    #print(f'DV/V = {np.abs(target_volume-vol)/(target_volume)*100}%')
    return index_List, flag

def Volume_evolution_Liouville_condition_concave_delaunay_refinements_test(w_tList,error_max,verbos=True):
    """
    def Volume_evolution_Liouville_condition_concave_delaunay_refinements(w_tList):
    take an list of set of points and use the Liouville theorem to comput the volume at each different time frame.
    The volume is compute for the first frime as the convexe hull of the set of points, in the other frame the algorithm will search the shape with same volume with delaunay concave refinement.
    Delaunay concave refinement start from the convexe hull and make the shape more and more concave by removing simplex with biggest circumcircle radius.
    input :
        w_tList : a list of time frame. Each time frame is a list of particle position in the form of 6 dimensional array
        float : error_max, pourcentage (in %) of maximum difference between volum output and target volume
    return :
        List of simplex in each time frame
        List of scipy Delaunay tesselation. One for each frame
        List of flag : 0 if convergence bewlow error_max, 1 otherwise
    """

    Volume = 0 # initial volume who will be the target value for all time frame volume
    tri_tList = [] # Delaunay tesselation of all time frame
    Nstep = len(w_tList)
    if verbos:
        print(f'DEBUG :: Volume_evolution_Liouville_condition_concave_delaunay_refinements_test :: Creation des triangulation de Delaunay')
    for tt in range(Nstep):
        tri_tList.append(Delaunay(w_tList[tt]))
    if verbos:
        print(f'DEBUG :: Volume_evolution_Liouville_condition_concave_delaunay_refinements_test :: Calcule du phase-space volume')
    tri0 = tri_tList[0]
    Volume = delaunay_volume_6D(tri0)
    if verbos:
        print(f'DEBUG :: Volume_evolution_Liouville_condition_concave_delaunay_refinements_test :: phase-space volume V = {Volume}')
    index_concave_delaunay_tList = []
    flag_tList = []

    for tt in range(Nstep):
        if tt ==0:
            index_concave_delaunay_tList.append(range(len(tri0.simplices)))
        else:
            if verbos:
                print(f'Calcul du concave delaunay frame {tt}')
            index_List,flag = Find_simplex_list_target_volume_concave_delaunay_refinement(Volume,tri_tList[tt],error_max)
            index_concave_delaunay_tList.append(index_List)
            flag_tList.append(flag)


    return index_concave_delaunay_tList,tri_tList,flag_tList


##############################################################################################
# Projection
##############################################################################################


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

def projection2d(vec,ax1,ax2):
    return [vec[ax1],vec[ax2]]

def simplex_projection2d(points,ax1,ax2):
    pro1 = [points[0][ax1],points[0][ax2]]
    pro2 = [points[1][ax1],points[1][ax2]]
    pro3 = [points[2][ax1],points[2][ax2]]
    pro4 = [points[3][ax1],points[3][ax2]]
    pro5 = [points[4][ax1],points[4][ax2]]
    pro6 = [points[5][ax1],points[5][ax2]]
    pro7 = [points[6][ax1],points[6][ax2]]
    return [pro1,pro2,pro3,pro4,pro5,pro6,pro7]

def projection_6d_to_2d(tri,id_list,ax1,ax2):
    tri2d_list = []
    V6d_list = []
    V2d_list = []
    jj=0
    n=len(id_list)
    for ii in id_list:
        #print(f'simplex {jj} of {n}', end='\r')
        jj=jj+1
        #if ii in id_list:
        V = Tetrahedron6DVolume(tri.points[tri.simplices[ii]])
        w2d = simplex_projection2d(tri.points[tri.simplices[ii]],ax1,ax2)
        tri2d = Delaunay(w2d)
        V2d = np.sum(delaunay_volume(tri2d))
        tri2d_list.append(tri2d)
        V6d_list.append(V)
        V2d_list.append(V2d)
    return tri2d_list,V6d_list,V2d_list

def see_projection(tri,id_list,ax1=0,ax2=1,resolution=256,boxeSize=1,trans=[0,0]):
    #print('computing projection')
    tri2d_list, V6d_list,V2d_list = projection_6d_to_2d(tri,id_list,ax1,ax2)
    #print(tri2d_list)
    data = np.zeros( (resolution,resolution))
    xList = np.linspace(-boxeSize,boxeSize,resolution)+trans[0]
    yList = np.linspace(-boxeSize,boxeSize,resolution)+trans[1]
    #xList = np.linspace(-75,50,resolution)
    #yList = np.linspace(-225,-100,resolution)
    #print('')
    q=0
    for ii in range(resolution):
        for jj in range(resolution):
            #print(f'pixel ({ii},{jj})', end='\r')
            for kk in range(len(tri2d_list)):
                tri = tri2d_list[kk]
                id_simplex = tri.find_simplex([xList[ii],yList[jj]])
                #print(id_simplex)
                if id_simplex>=0:
                    #print(tri.points[tri.simplices[id_simplex]])
                    V2d = triangleArea(*tri.points[tri.simplices[id_simplex]])
                    q = (V6d_list[kk])*V2d/V2d_list[kk]
                    #print(q)
                    data[resolution-1-jj,ii]+=q
    return data


##############################################################################################
# Notebook
##############################################################################################



def download_simulation_v2():
    nb_tList = []
    Nsim = 67
    for ii in range(Nsim):
        if ii in [0,66]:
            txt_insert = ''
            if ii<10:
                txt_insert = f'000{ii}'
            elif ii<100:
                txt_insert = f'00{ii}'
            elif ii<1000:
                txt_insert = f'0{ii}'

            nb_tList.append(Nbody(f'/srv/astro/projects/clastro/revaz/ARRAKIHS/SWIFT_UFDS_TidalStripping/Sims/DW103/snap/snapshot_{txt_insert}_stars.hdf5',ftype='swift'))

    for tt in range(len(nb_tList)):
        nb_tList[tt].translate([-500,-500,-500])
    return nb_tList

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

def unscale_w(points_tList,nb_tList):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = compute_std(nb_tList[0])
    #ans = []
    for tt in range(len(points_tList)):
        #tmp1 = []
        for ii in range(len(points_tList[tt])):
            #tmp2 = [0,0,0,0,0,0]
            points_tList[tt][ii][0] = points_tList[tt][ii][0]*std_x
            points_tList[tt][ii][1] = points_tList[tt][ii][1]*std_y
            points_tList[tt][ii][2] = points_tList[tt][ii][2]*std_z
            points_tList[tt][ii][3] = points_tList[tt][ii][3]*std_vx
            points_tList[tt][ii][4] = points_tList[tt][ii][4]*std_vy
            points_tList[tt][ii][5] = points_tList[tt][ii][5]*std_vz
            #tmp2[0] = points_tList[tt][ii][0]*std_x
            #tmp2[1] = points_tList[tt][ii][1]*std_y
            #tmp2[2] = points_tList[tt][ii][2]*std_z
            #tmp2[3] = points_tList[tt][ii][3]*std_vx
            #tmp2[4] = points_tList[tt][ii][4]*std_vy
            #tmp2[5] = points_tList[tt][ii][5]*std_vz
            #tmp1.append(tmp2)
        #ans.append(tmp1)
    return points_tList

def unscale_tri(tri_tList,nb_tList):
    std_x,std_y,std_z,std_vx,std_vy,std_vz = compute_std(nb_tList[0])
    for tt in range(len(tri_tList)):
        for ii in range(len(tri_tList[tt].points)):
            tri_tList[tt].points[ii][0] = tri_tList[tt].points[ii][0]*std_x
            tri_tList[tt].points[ii][1] = tri_tList[tt].points[ii][1]*std_y
            tri_tList[tt].points[ii][2] = tri_tList[tt].points[ii][2]*std_z
            tri_tList[tt].points[ii][3] = tri_tList[tt].points[ii][3]*std_vx
            tri_tList[tt].points[ii][4] = tri_tList[tt].points[ii][4]*std_vy
            tri_tList[tt].points[ii][5] = tri_tList[tt].points[ii][5]*std_vz
    return tri_tList

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

def distance6D(x1,x2):
    return np.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2+(x1[3]-x2[3])**2+(x1[4]-x2[4])**2+(x1[5]-x2[5])**2)

def create_tri0(id_trunc):
    print('scaling data')
    nb_scale_tList = scale_nb_tList(download_simulation_v2())
    print('make full tessellation')
    tri0 = Delaunay(PosVel_to_w(nb_scale_tList[0].selectp(id_trunc)))
    return tri0

def create_w_tList_cell_v2(nb_tList,id_trunc,id_particle,tri0):
    
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][id_particle]:tri0.vertex_neighbor_vertices[0][id_particle+1]]

    nb_id = []
    for ii in id_particle_in_cell:
        nb = nb_tList[0].selectp(id_trunc)
        nb_id.append(nb.num[ii])

    
    w_tList = []
    for tt in range(len(nb_tList)):
        nb = nb_tList[tt].selectp(id_trunc)
        w_tList.append(PosVel_to_w(nb.selectp(nb_id)))
    return w_tList

def create_w_tList_cell_v3(nb_tList,id_trunc,id_particle,tri0,Rmax):
    
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][id_particle]:tri0.vertex_neighbor_vertices[0][id_particle+1]]

    nb_id = []
    for ii in id_particle_in_cell:
        if distance6D(tri0.points[id_particle],tri0.points[ii])<Rmax:
            nb = nb_tList[0].selectp(id_trunc)
            nb_id.append(nb.num[ii])
    
    w_tList = []
    for tt in range(len(nb_tList)):
        nb = nb_tList[tt].selectp(id_trunc)
        w_tList.append(PosVel_to_w(nb.selectp(nb_id)))
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

def add_mid_point(index_concave_delaunay,tri):
    mid_point = []
    for ii in index_concave_delaunay:
        center = find_center_6D(tri.points[tri.simplices[ii]])
        vertex = tri.points[tri.simplices[ii]]
        for ii in range(7):
            mid_point.append(1/2*(vertex[ii]+center))
        mid_point.append(center)
    for ii in range(len(tri.points)):
        mid_point.append(tri.points[ii])
    return mid_point

def find_cell_points_tList(index_concave_delaunay_tList,tri_tList):
    cell_position_tList = []
    for tt in range(len(tri_tList)):
        cell_position_tList.append(add_mid_point(index_concave_delaunay_tList[tt],tri_tList[tt]))
    return cell_position_tList









##############################################################################################
# Dernier ajout
##############################################################################################

def create_w_tList_cell_v4(nb_tList,id_particle,tri0,Rmax):
    
    #tri0 = Delaunay(PosVel_to_w(nb_tList[0]))
    id_particle_in_cell = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][id_particle]:tri0.vertex_neighbor_vertices[0][id_particle+1]]

    nb_id = []
    for ii in id_particle_in_cell:
        if distance6D(tri0.points[id_particle],tri0.points[ii])<Rmax:
            nb_id.append(nb_tList[-1].num[ii])
    
    w_tList = []
    for tt in range(len(nb_tList)):
        w_tList.append(PosVel_to_w(nb_tList[tt].selectp(nb_id)))
    return w_tList


##############################################################################################
# Fonction en parrallele
##############################################################################################
def makeImageMatrixOneparticule(p_id,nb_scale_tList,tri0_scale,Rmax):
    w_scale_tList = create_w_tList_cell_v4(nb_scale_tList,p_id,tri0_scale,Rmax)

    #print('DEBUG : Concave Delaunay Algorithm')
    index_conserved_simplex_tList,tri_tList,flag_tList = Volume_evolution_Liouville_condition_concave_delaunay_refinements_test(w_scale_tList,error_max=5,verbos=False)
    flag = flag_tList[-1]

    resolution = 20

    if flag==0:
        data = see_projection(tri_tList[-1],index_conserved_simplex_tList[-1],ax1=0,ax2=2,resolution=resolution,boxeSize=90,trans=[0,0])
        sum = np.sum(data)
        if not sum==0:
            data = data/sum
        return data
    return np.zeros((resolution,resolution))