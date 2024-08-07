# ConcaveDelaunayMockImage

Code pour utilise crée des mock image depuis des Nbody simulation. 
Tout les paramètre doivent être rentrée dans paraImage.yaml

SimulationName: La snapshot utiliser pour faire l'image
InitialFrameName: La situation inital de la simulation (snapshot a t=0)
translation: translation du système pour que la distribution de particules soit centrée sur l'origine
triSaveName: Path (location + nom) du fichier contenant la triangulation de Delaunay
NeighborsSaveName: Path (location + nom) du fichier contenant la liste de voisins de chaque particules
triangleFolder: Folder (path du fichier) ou les fichier triangles (des fichier contenant les donner utiliser pour dessiner l'image final) seront sauvegarder

axe1: premier axe de l'image
axe2: second axe de l'image (1=x, 2=y, 3=z)
halfBoxSize: taille de l'image (en unité de la simualtion)

Rmax: paramètre déterminant la taille maximum des cellules crée dans l'algorithm (le plus grand possible est la meilleur option mais cela peut causer des problèmes numérique) 
Nmin: nombre minimum de particules dans une cellules pour la considéré comme valide

ImagePath: Path (location+nom) de l'image qui sera crée


Le script CreateTri.py doit être lancé en premier.
D'une fois se fichier ayant terminer sont exécution, PartialImage.py doit être lancer avec les paramètre de 0 à 9
python CreateTri.py
d'une fois terminer 
python PartialImage.py 0
python PartialImage.py 1
python PartialImage.py 2
python PartialImage.py 3
....
python PartialImage.py 9

cela créera tout les fichier triangles nécessaire à la création de l'image
d'une fois tout les fichier triangle crée, exécuter BuildImage.py
python BuildImage.py