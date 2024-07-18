import numpy as np
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from pNbody import*
from toolMockImage import*
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

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



