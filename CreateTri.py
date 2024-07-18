from toolMockImage import*
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import pickle 
import yaml

#===============================================================================================
# Computation of the Delaunay triangulation
#===============================================================================================

print('starting')
print('Download parameter')

with open("paraImage.yaml") as stream:
    try:
        data = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

name = data["SimulationName"]
trans = data["translation"]
rot = data["rotation"]

print('Download Data')

#Download data
nb = download_simulation(name,trans,rot)
nb_scaled,std = scale_nb(nb)

print('Make triangulation')
tri0 = Delaunay(PosVel_to_w(nb_scaled))

print('Delaunay triangulation done')
file_path = data["triSaveName"]

print('save triangulation')
with open(file_path, 'wb') as file:
    # Serialize and write the variable to the file
    pickle.dump(tri0, file)

print("The Delaunay triangulation has been saved successfully.")

#===============================================================================================
#   Computation of Neighbors of each node
#===============================================================================================

print('Computation of Neighbors list')

neighbors = find_neighbors(tri0)

file_path_neighbors = data["NeighborsSaveName"]

with open(file_path_neighbors, 'wb') as file2:
    pickle.dump(neighbors, file2)

print("The Neighbors list has been saved successfully.")