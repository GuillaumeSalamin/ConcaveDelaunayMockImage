import numpy as np
from scipy.spatial import Delaunay

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
    print(f'DV/V = {np.abs(target_volume-vol)/(target_volume)*100}%')
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