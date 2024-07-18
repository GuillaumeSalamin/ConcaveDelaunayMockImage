import sys
sys.path.append('/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript')
from toolMockImage import*
import numpy as np
import pickle 

print('starting')
print('Download Data')

#Download data
nb_scale_tList = scale_nb_tList(download_simulation_v2())

print('Make triangulation')
tri0 = Delaunay(PosVel_to_w(nb_scale_tList[0]))

print('Delaunay triangulation done')
file_path = '/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/DelaunayTri/tri0.pickle'

print('save triangulation')
with open(file_path, 'wb') as file:
    # Serialize and write the variable to the file
    pickle.dump(tri0, file)

print("The variable 'data' has been saved successfully.")