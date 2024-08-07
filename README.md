# ConcaveDelaunayMockImage

Code pour utilise crée des mock image depuis des Nbody simulation. 
Tout les paramètre doivent être rentrée dans paraImage.yaml

### paraImage.yaml <br />
SimulationName: La snapshot utiliser pour faire l'image <br />
InitialFrameName: La situation inital de la simulation (snapshot a t=0)<br />
translation: translation du système pour que la distribution de particules soit centrée sur l'origine<br />
triSaveName: Path (location + nom) du fichier contenant la triangulation de Delaunay<br />
NeighborsSaveName: Path (location + nom) du fichier contenant la liste de voisins de chaque particules<br />
triangleFolder: Folder (path du fichier) ou les fichier triangles (des fichier contenant les donner utiliser pour dessiner l'image final) seront sauvegarder<br />
axe1: premier axe de l'image<br />
axe2: second axe de l'image (1=x, 2=y, 3=z)<br />
halfBoxSize: taille de l'image (en unité de la simualtion)<br />

Rmax: paramètre déterminant la taille maximum des cellules crée dans l'algorithm (le plus grand possible est la meilleur option mais cela peut causer des problèmes numérique) 
Nmin: nombre minimum de particules dans une cellules pour la considéré comme valide

ImagePath: Path (location+nom) de l'image qui sera crée

### Comment exécuter le code<br />
Le script CreateTri.py doit être lancé en premier.<br />
D'une fois se fichier ayant terminer sont exécution, PartialImage.py doit être lancer avec les paramètre de 0 à 9<br />
python CreateTri.py<br />
d'une fois terminer <br />
python PartialImage.py 0<br />
python PartialImage.py 1<br />
python PartialImage.py 2<br />
python PartialImage.py 3<br />
....<br />
python PartialImage.py 9<br />
<br />
cela créera tout les fichier triangles nécessaire à la création de l'image
d'une fois tout les fichier triangle crée, exécuter BuildImage.py
python BuildImage.py
