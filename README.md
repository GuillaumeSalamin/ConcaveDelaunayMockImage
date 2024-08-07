# ConcaveDelaunayMockImage

Code pour utilise crée des mock images depuis des Nbody simulation. 
Tous les paramètres doivent être rentrés dans paraImage.yaml

### paraImage.yaml <br />
SimulationName: La snapshot utiliser pour faire l'image <br />
InitialFrameName: La situation initiale de la simulation (snapshot a t=0)<br />
translation: translation du système pour que la distribution de particules soit centrée sur l'origine<br />
triSaveName: Path (location + nom) du fichier contenant la triangulation de Delaunay<br />
NeighborsSaveName: Path (location + nom) du fichier contenant la liste de voisins de chaque particule<br />
triangleFolder: Folder (path du fichier) ou les fichiers triangles (des fichiers contenant les données utilisées pour dessiner l'image finale) seront sauvegardés<br />
axe1: premier axe de l'image<br />
axe2: second axe de l'image (1=x, 2=y, 3=z)<br />
halfBoxSize: taille de l'image (en unité de la simulation)<br />

Rmax: paramètre déterminant la taille maximum des cellules créée dans l'algorithme (le plus grand possible est la meilleure option, mais cela peut causer des problèmes numériques) 
Nmin: nombre minimum de particules dans une cellule pour la considérer comme valide

ImagePath: Path (location+nom) de l'image qui sera créée

### Comment exécuter le code<br />
Le script CreateTri.py doit être lancé en premier.<br />
D'une fois se fichier ayant terminé sont exécution, PartialImage.py doit être lancé avec les paramètres de 0 à 9<br />
python CreateTri.py<br />
d'une fois terminer <br />
python PartialImage.py 0<br />
python PartialImage.py 1<br />
python PartialImage.py 2<br />
python PartialImage.py 3<br />
....<br />
python PartialImage.py 9<br />
<br />
cela créera tous les fichiers triangles nécessaires à la création de l'image
d'une fois tous les fichiers triangle crée, exécuter BuildImage.py<br />
python BuildImage.py

## Prérequis<br />
Python3<br />
pNbody<br />
numpy<br />
matplolib<br />
scipy<br />
pickle<br />
yaml<br />
glfw<br />
pyOpenGl<br />
