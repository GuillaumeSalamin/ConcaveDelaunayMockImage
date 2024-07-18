import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import pickle

file_path_tri = '/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/DelaunayTri/tri0.pickle'
with open(file_path_tri, 'rb') as file:
    # Deserialize and retrieve the variable from the file
    loaded_data = pickle.load(file)

print("The variable 'data' has been loaded successfully.")

tri0 = loaded_data

from collections import defaultdict

def find_neighbors(tess):
    neighbors = defaultdict(set)

    for simplex in tess.simplices:
        for idx in simplex:
            other = set(simplex)
            other.remove(idx)
            neighbors[idx] = neighbors[idx].union(other)
    return neighbors

def find_neighbors_v2(tri):
    ans = []
    for ii in range(len(tri.points)):
        tmp = tri0.vertex_neighbor_vertices[1][tri0.vertex_neighbor_vertices[0][ii]:tri0.vertex_neighbor_vertices[0][ii+1]]
        ans.append(tmp)
    return ans

neighbors = find_neighbors_v2(tri0)

file_path = '/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonSaveVariable/neighbors_v2.pickle'

with open(file_path, 'wb') as file2:
    # Serialize and write the variable to the file
    pickle.dump(neighbors, file2)

print("The variable 'data' has been saved successfully.")