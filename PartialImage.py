import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*

from DrawFunction import*

from concaveDelaunayRefinement import*
from toolMockImage import*

import pickle
import yaml
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

import sys


if __name__ == "__main__":
    
    with open("paraImage.yaml") as stream:
        try:
            para = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    name = para["SimulationName"]
    trans = para["translation"]
    rot = para["rotation"]
    triName = para["triSaveName"]
    NeighborsName = para["NeighborsSaveName"]
    InitialFrame = para["InitialFrameName"]
    MatrixFolder = para["MatrixFolder"]
    axe1 = para["axe1"]
    axe2 = para["axe2"]
    Nmin = para["Nmin"]
    Rmax = para["Rmax"]
    halfBoxSize = para["halfBoxSize"]


    print('scaling data')
    nb0 = download_simulation(name,trans,rot)
    std_list = compute_std(nb0)
    nb0_scaled = scale_nb(nb0,std_list)

    nb = download_simulation(name,trans,rot)
    nb_scaled = scale_nb(nb,std_list)
   
    nb_scaled_tList = [nb0_scaled,nb_scaled]

    file_path = NeighborsName
    with open(file_path, 'rb') as file:
        # Deserialize and retrieve the variable from the file
        loaded_data = pickle.load(file)

    print("The variable 'Neighbors' has been loaded successfully.")

    neighbors = loaded_data

    file_path = triName
    with open(file_path, 'rb') as file:
        # Deserialize and retrieve the variable from the file
        loaded_data = pickle.load(file)

    print("The variable 'Tri' has been loaded successfully.")

    tri0 = loaded_data

    filePara = int(sys.argv[1])

    

    parameter_a_list = []
    for ii in range(10):
        parameter_a_list.append(filePara*10+ii)

    for parameter_a in parameter_a_list:
        cell_id_list = []
        ImageMatrix = np.zeros((600,800))
        for ii in range(len(nb_scaled_tList[0].num)):
            if ii%100==parameter_a:
                cell_id_list.append(ii)

        n = len(cell_id_list)
        jj=0

        
        for p_id in cell_id_list:
            jj=jj+1
            print(f'DEBUG :: iteration {jj} of {len(cell_id_list)}')
            print(f'DEBBUG :: cell creation')
            w_scale_tList = create_w_tList_cell(nb_scaled_tList,neighbors,p_id,Rmax=Rmax)
            flag = 0
            if len(w_scale_tList[0])<Nmin:
                flag=1
            else:
                print('DEBUG : Concave Delaunay Algorithm')
                index_conserved_simplex_tList,tri_tList,flag_tList = Volume_evolution_Liouville_condition_concave_delaunay_refinements(w_scale_tList,error_max=5,verbos=False)
                flag = flag_tList[-1]
                

                index_List = index_conserved_simplex_tList[-1]
                tri = tri_tList[-1]
                V_cons = delaunay_volume_6D(tri_tList[0])

                if flag==0:
                    tmp,opa = VertexList2d(tri,std_list,index_List,axe1,axe2,halfBoxSize)
                    opa = opa/V_cons
                    for ii in range(len(tmp)):
                        mat = triangleDraw(tmp[ii],opa[ii])
                        ImageMatrix +=mat
                    
            if flag==1:
                print('flag=1')

        name = f'Partial_Image_{parameter_a}'

        file_path_data = f'{MatrixFolder}{name}_data.pickle'
        with open(file_path_data, 'wb') as file:
            # Serialize and write the variable to the file
            pickle.dump(ImageMatrix, file)

        print(f'The variable "data" has been saved successfully. \n File name :: {file_path_data}')