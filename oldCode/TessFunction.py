import numpy as np
from scipy.spatial import Delaunay

def dist6d(px1,py1,pz1,vx1,vy1,vz1,px2,py2,pz2,vx2,vy2,vz2,nb):

   stdpx = nb.x().std()
   stdpy = nb.y().std()
   stdpz = nb.z().std()
   stdvx = nb.vx().std()
   stdvy = nb.vy().std()
   stdvz = nb.vz().std()

   print(stdpx,stdpy,stdpz,stdvx,stdvy,stdvz)

   dpx2 = ((px1-px2)/stdpx)**2
   dpy2 = ((py1-py2)/stdpy)**2
   dpz2 = ((pz1-pz2)/stdpz)**2
   dvx2 = ((vx1-vx2)/stdvx)**2
   dvy2 = ((vy1-vy2)/stdvy)**2
   dvz2 = ((vz1-vz2)/stdvz)**2

   d = np.sqrt(dpx2 + dpy2 + dpz2 + dvx2 + dvy2 + dvz2)

   return d

def NRandom(N,d):
    w = np.zeros((N,d))
    for ii in range(d):
        w[:,ii] = np.random.random(N)
    return w

from astropy import units as u

def NRandom2D(N,xscale=1,vscale=1,mscale=1):
    """Create a set of N particle with random position, velocity and mass
    take as argument
    -N number of particle created
    -x scale for the random generation of position (pc)
    -v scale for the random generation of velocity (km/s)
    -m scale for the random generation of mass"""

    vscale = vscale*u.km.to(u.pc)/u.s.to(u.yr)

    x = (np.random.random_sample(N)-0.5)*xscale
    y = (np.random.random_sample(N)-0.5)*xscale
    vx = (np.random.random_sample(N)-0.5)*vscale
    vy = (np.random.random_sample(N)-0.5)*vscale
    m = np.random.random_sample(N)*mscale
    return x,y,vx,vy,m

def NRandom3D(N,xscale=1,vscale=1,mscale=1):
    """Create a set of N particle with random position, velocity and mass
    take as argument
    -N number of particle created
    -x scale for the random generation of position
    -v scale for the random generation of velocity
    -m scale for the random generation of mass"""
    x = [] 
    v = [] 
    m = []
    for ii in range(N):
        x.append(np.random.random_sample(3)*xscale)
        v.append((np.random.random_sample(3)-0.5)*vscale)
        m.append(np.random.random_sample(1)[0]*mscale)
    return x,v,m

def createW(x,y):
    """Combine two vector in the good format to use scipy.delauny or scipy.voronoi"""
    N=len(x)
    w = np.zeros((N,2))
    for ii in range(N):
        w[ii,0]=x[ii]
        w[ii,1]=y[ii]
    return w


##############################################################################################
# Fonction to change the format of the different object containing data
##############################################################################################
def MakeShapeForPlot(Tab):
    Ndim = len(Tab[0])
    x=[]
    y=[]
    z=[]
    for step in Tab:
        xLine=[]
        yLine=[]
        zLine=[]
        for object in step:
            if Ndim>0:
                xLine.append(object[0])
                if Ndim>1:
                    yLine.append(object[1])
                    if Ndim>2:
                        zLine.append(object[2])
        x.append(xLine)
        y.append(yLine)
        z.append(zLine) 
    return x,y,z
    
##############################################################################################
# Fonction to make gif and display Voronoi/Delaunay tesselation
##############################################################################################

from PIL import Image
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

def fig2img(fig):
    """Convert a Matplotlib figure to a PIL Image and return it"""
    import io
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    return img

def draw1FrameFullInfo(ii,xList,yList,vxList,vyList,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    """def draw1FrameFullInfo(ii,xList,yList,vxList,vyList,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    Draw an image with 4 graphe, two position space Delaunay/Voronoi and two velocity space Delaunay/Vornoi, return only one pillow image
    take in argument, 
    - the index of the frame you want to draw
    - X vector position x of each particle
    - Y vector position y of each particle
    - VX vector velocity x of each particle
    - VY vector velocity y of each particle
    - the for next argument are the limit value of the plot (define the window draw)
    """
    x,y,vx,vy = xList[ii],yList[ii],vxList[ii],vyList[ii]
    wX = createW(x,y)
    wV = createW(vx,vy)

    triX = Delaunay(wX)
    triV= Delaunay(wV)
    vorX = Voronoi(wX)
    vorV = Voronoi(wV)

    fig, axs = plt.subplots(2, 2, figsize=(15, 10))

    axs[0,0].triplot(wX[:,0], wX[:,1], triX.simplices)
    axs[0,0].plot(wX[:,0], wX[:,1], 'o')
    axs[0,0].set_title('Delauny position_x-position_y tesselation')
    axs[0,0].set_xlabel('position x')
    axs[0,0].set_ylabel('position y')
    axs[0,0].set_xlim(xlimM,xlimP)
    axs[0,0].set_ylim(xlimM,xlimP)

    voronoi_plot_2d(vorX, ax=axs[0,1], show_vertices=False, line_colors='blue',line_width=2, line_alpha=0.6, point_size=2)
    axs[0,1].scatter(wX[:, 0], wX[:, 1], c='green', label='Points', zorder=5,s=5)
    axs[0,1].set_title('Voronoi position_x-position_y tesselation')
    axs[0,1].set_xlabel('position x')
    axs[0,1].set_ylabel('position y')
    axs[0,1].set_xlim(xlimM,xlimP)
    axs[0,1].set_ylim(xlimM,xlimP)

    axs[1,0].triplot(wV[:,0], wV[:,1], triV.simplices)
    axs[1,0].plot(wV[:,0], wV[:,1], 'o')
    axs[1,0].set_title('Delauny velocity_x-velocity_y tesselation')
    axs[1,0].set_xlabel('velocity x')
    axs[1,0].set_ylabel('velocity y')
    axs[1,0].set_xlim(vlimM,vlimP)
    axs[1,0].set_ylim(vlimM,vlimP)

    voronoi_plot_2d(vorV, ax=axs[1,1], show_vertices=False, line_colors='blue',line_width=2, line_alpha=0.6, point_size=2)
    axs[1,1].scatter(wV[:, 0], wV[:, 1], c='green', label='Points', zorder=5,s=5)
    axs[1,1].set_title('Voronoi velocity_x-velocity_y tesselation')
    axs[1,1].set_xlabel('velocity x')
    axs[1,1].set_ylabel('velocity y')
    axs[1,1].set_xlim(vlimM,vlimP)
    axs[1,1].set_ylim(vlimM,vlimP)

    img = fig2img(fig)
    plt.close()
    return  img
    #fig.savefig(f'/home/guillaume/TPIV/image/gif/frame{ii}.png')

def makeGifFullInfo(xList,yList,vxList,vyList,name,Nframe=100,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    """makeGifFullInfo(xList,yList,vxList,vyList,name,Nframe=100,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7)
    Take the result of a simulation (a table (or a np_array of np_array) of position and velocity) and make a gif from it
    argument are
    - Table of x position
    - Table of y position
    - Table of x velocity
    - Table of y velocity
    - name of the gif file create (don't forget to put .gif at the end)
    - number of frame draw (decrease to draw faster, increase to have a smooth gif)
    - 4 last argument are the window draw in position and velocity space (limit value for the axis)"""
    img0 = draw1FrameFullInfo(0,xList,yList,vxList,vyList)
    
    imgList = [img0]

    #freq = ceil(len(xList)/Nframe)
    freq = len(xList)//Nframe

    if freq<1:
        freq=1

    #count = 0
    for ii in range(len(xList)):
        if ii%freq==0:
            #print(f'DEBUG :: ii={ii}')
            #count=count+1
            imgList.append(draw1FrameFullInfo(ii,xList,yList,vxList,vyList,xlimM,xlimP,vlimM,vlimP))
    #imgList.append(draw1FrameFullInfo(-1,xList,yList,vxList,vyList,xlimM,xlimP,vlimM,vlimP))
    #print(f'DEBUG :: total frame used={count}')
    img0.save(name, save_all=True, append_images=imgList)

def draw1FrameRealSpace(ii,xList,yList,axeLim):
    x,y = xList[ii],yList[ii]

    fig, axs = plt.subplots(figsize=(3, 3))


    axs.plot(x, y, 'o')
    axs.set_title('Real space particle')
    axs.set_xlabel('position x [pc]')
    axs.set_ylabel('position y [pc]')
    axs.set_xlim(-axeLim,axeLim)
    axs.set_ylim(-axeLim,axeLim)


    img = fig2img(fig)
    plt.close()
    return  img

def makeGifRealSpace(xList,yList,name,Nframe=100,axeLim=1):
    img0 = draw1FrameRealSpace(0,xList,yList,axeLim)
    
    imgList = [img0]

    #freq = ceil(len(xList)/Nframe)
    freq = len(xList)//Nframe

    if freq<1:
        freq=1

    #count = 0
    for ii in range(len(xList)):
        if ii%freq==0:
            #print(f'DEBUG :: ii={ii}')
            #count=count+1
            imgList.append(draw1FrameRealSpace(ii,xList,yList,axeLim))
    img0.save(name, save_all=True, append_images=imgList)

def createImageRealSpace(ii,xList,yList):
    x,y = xList[ii],yList[ii]

    fig, axs = plt.subplots(figsize=(4, 4))

    for ii in range(len(x)):
        axs.plot(x[ii], y[ii], 'o')
    axs.set_xlabel('position x [pc]')
    axs.set_ylabel('position y [pc]')


    img = fig2img(fig)
    plt.close()
    return  img

def createListRealTrajectorySingleParticleXY(ParticleIndex,XTab):
    Nstep = len(XTab)
    x=[]
    y=[]
    for ii in range(Nstep):
        x.append(XTab[ii][ParticleIndex][0])
        y.append(XTab[ii][ParticleIndex][1])

    return x,y


##############################################################################################
# Voronoi and Delaunay density
##############################################################################################
from scipy.spatial import ConvexHull
import matplotlib as mpl
import matplotlib.cm as cm

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

from scipy.spatial import ConvexHull

"""La fonct tion Voronoi:volume se trouvait ici"""
def voronoi_volumes(vor):
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
        indices = vor.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(vor.vertices[indices]).volume
            #print(f'DEBUG :: vor.vertices[indices]={vor.vertices[indices]}')
            #print(np.shape(vor.vertices[indices]))
            #vol[i] = volume6D(vor.vertices[indices])
    return vol


from PIL import Image

def PlotVoronoiDensity(vor,log_scale=True,cmap=cm.Blues_r):
    volVoronoi = voronoi_volumes(vor)
    if log_scale:
        for ii in range(len(volVoronoi)):
            volVoronoi[ii] = np.log(volVoronoi[ii])

    # find min/max values for normalization
    minima = min(volVoronoi)
    maxima = max(volVoronoi[volVoronoi<1e10])#attention cette ligne de code peut poser problème dans certain cas extrème

    # normalize chosen colormap
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig, axs = plt.subplots()

    # plot Voronoi diagram, and fill finite regions with color mapped from vol value
    voronoi_plot_2d(vor, ax=axs, show_points=True, show_vertices=False, line_colors='blue',line_width=2, line_alpha=0.6, point_size=2)
    plt.close()
    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            axs.fill(*zip(*polygon), color=mapper.to_rgba(volVoronoi[r]))
    return fig

def DrawVoronoiDensity(vor,log_scale=True):
    img = fig2img(PlotVoronoiDensity(vor,log_scale=log_scale))
    plt.close()
    return img

def PlotDelaunayDensity(tri,points,log_scale=True,cmap=cm.Blues_r):
    volDelaunay = delaunay_volume(tri)

    if log_scale:
        for ii in range(len(volDelaunay)):
            volDelaunay[ii] = np.log(volDelaunay[ii])

    # find min/max values for normalization
    minima = min(volDelaunay)
    maxima = max(volDelaunay)

    # normalize chosen colormap
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    fig, axs = plt.subplots()

    # plot Voronoi diagram, and fill finite regions with color mapped from vol value
    axs.triplot(points[:,0], points[:,1], tri.simplices)
    axs.plot(points[:,0], points[:,1], 'o',markersize=1)
    for rr in range(len(tri.simplices)):
        polygon = plt.Polygon(points[tri.simplices[rr]],color=mapper.to_rgba(volDelaunay[rr]))
        axs.add_patch(polygon)
    return fig

def DrawDelaunayDensity(vor,log_scale=True):
    img = fig2img(PlotDelaunayDensity(vor,log_scale=log_scale))
    plt.close()
    return img


##############################################################################################
# Voronoi and Delaunay density -> draw fram
##############################################################################################
def draw1FrameFullInfoDensity(ii,xList,yList,vxList,vyList,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    """def draw1FrameFullInfo(ii,xList,yList,vxList,vyList,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    Draw an image with 4 graphe, two position space Delaunay/Voronoi and two velocity space Delaunay/Vornoi, return only one pillow image
    take in argument, 
    - the index of the frame you want to draw
    - X vector position x of each particle
    - Y vector position y of each particle
    - VX vector velocity x of each particle
    - VY vector velocity y of each particle
    - the for next argument are the limit value of the plot (define the window draw)
    """
    x,y,vx,vy = xList[ii],yList[ii],vxList[ii],vyList[ii]
    wX = createW(x,y)
    wV = createW(vx,vy)

    triX = Delaunay(wX)
    triV= Delaunay(wV)
    vorX = Voronoi(wX)
    vorV = Voronoi(wV)

    fig, axs = plt.subplots(2, 2, figsize=(15, 10))

    AxsVoronoiDensity(axs[0,1],vorX,cmap='inferno')
    AxsDelaunayDensity(axs[0,0],triX,wX,cmap='inferno')
    AxsVoronoiDensity(axs[1,1],vorV,cmap='inferno')
    AxsDelaunayDensity(axs[1,0],triV,wV,cmap='inferno')

    axs[0,0].set_title('Delauny position_x-position_y tesselation')
    axs[0,0].set_xlabel('position x')
    axs[0,0].set_ylabel('position y')
    axs[0,0].set_xlim(xlimM,xlimP)
    axs[0,0].set_ylim(xlimM,xlimP)


    axs[0,1].set_title('Voronoi position_x-position_y tesselation')
    axs[0,1].set_xlabel('position x')
    axs[0,1].set_ylabel('position y')
    axs[0,1].set_xlim(xlimM,xlimP)
    axs[0,1].set_ylim(xlimM,xlimP)


    axs[1,0].set_title('Delauny velocity_x-velocity_y tesselation')
    axs[1,0].set_xlabel('velocity x')
    axs[1,0].set_ylabel('velocity y')
    axs[1,0].set_xlim(vlimM,vlimP)
    axs[1,0].set_ylim(vlimM,vlimP)

 
    axs[1,1].set_title('Voronoi velocity_x-velocity_y tesselation')
    axs[1,1].set_xlabel('velocity x')
    axs[1,1].set_ylabel('velocity y')
    axs[1,1].set_xlim(vlimM,vlimP)
    axs[1,1].set_ylim(vlimM,vlimP)

    img = fig2img(fig)
    plt.close()
    return  img


def makeGifFullInfoDensity(xList,yList,vxList,vyList,name,Nframe=100,xlimM=-2,xlimP=4,vlimM=-7,vlimP=7):
    img0 = draw1FrameFullInfoDensity(0,xList,yList,vxList,vyList,xlimM=xlimM,xlimP=xlimP,vlimM=vlimM,vlimP=vlimP)
    
    imgList = [img0]
    freq = len(xList)//Nframe

    if freq<1:
        freq=1
    for ii in range(len(xList)):
        if ii%freq==0:

            imgList.append(draw1FrameFullInfoDensity(ii,xList,yList,vxList,vyList,xlimM,xlimP,vlimM,vlimP))
    img0.save(name, save_all=True, append_images=imgList)

##############################################################################################
# Voronoi and Delaunay density -> fonction experimental
##############################################################################################

def AxsDelaunayDensity(axs,tri,points,log_scale=True,cmap=cm.Blues_r,display_point=True):
    volDelaunay = delaunay_volume(points,tri)

    if log_scale:
        for ii in range(len(volDelaunay)):
            volDelaunay[ii] = -np.log(volDelaunay[ii])

    # find min/max values for normalization
    minima = min(volDelaunay)
    maxima = max(volDelaunay)

    # normalize chosen colormap
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    # plot Voronoi diagram, and fill finite regions with color mapped from vol value
    axs.triplot(points[:,0], points[:,1], tri.simplices)
    if display_point:
        axs.plot(points[:,0], points[:,1], 'o',markersize=1)
    for rr in range(len(tri.simplices)):
        polygon = plt.Polygon(points[tri.simplices[rr]],color=mapper.to_rgba(volDelaunay[rr]))
        axs.add_patch(polygon)

def AxsVoronoiDensity(axs,vor,log_scale=True,cmap=cm.Blues_r,display_point=True,display_vertice=True):
    volVoronoi = voronoi_volumes(vor)
    if log_scale:
        for ii in range(len(volVoronoi)):
            volVoronoi[ii] = np.log(volVoronoi[ii])

    # find min/max values for normalization
    minima = min(volVoronoi)
    maxima = max(volVoronoi[volVoronoi<1e10])#attention cette ligne de code peut poser problème dans certain cas extrème

    # normalize chosen colormap
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    # plot Voronoi diagram, and fill finite regions with color mapped from vol value
    lineWidth = 2.0
    if not display_vertice:
        lineWidth=0.0
    voronoi_plot_2d(vor, ax=axs, show_points=display_point, show_vertices=display_vertice, line_colors='blue',line_width=lineWidth, line_alpha=0.6, point_size=2)
    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            axs.fill(*zip(*polygon), color=mapper.to_rgba(-volVoronoi[r]))

##############################################################################################
# Voronoi and Delaunay density -> fonction prenant en compt le nombre de triangle autour de chaque points et le nombre de point par triangle
##############################################################################################
"""def number_of_triangle_per_points(tri):
    Npoints = len(tri.points)
    Ntriangle = len(tri.simplices)
    trianglePerPoints = np.zeros(Npoints)
    for ii in range(Npoints):
        PointID=ii
        ans = 0
        for jj in range(Ntriangle):
            for kk in range(len(tri.simplices[jj])):
                if PointID==tri.simplices[jj][kk]:
                    ans = ans+1
        trianglePerPoints[ii]=ans 
    return trianglePerPoints


def points_per_triangle(tri):
    NumberPoints = number_of_triangle_per_points(tri)
    Ntriangle = len(tri.simplices)
    NpointsFraq = np.zeros(Ntriangle)
    for ii in range(len(tri.simplices)):
        points_fraq = 0
        for kk in range(len(tri.simplices[ii])):
            points_fraq = points_fraq+1/NumberPoints[tri.simplices[ii][kk]]
        NpointsFraq[ii] = points_fraq
    return NpointsFraq"""

#========================================================================================
# nouvelle version des fonction qui était en N^2
#========================================================================================

def find_neighbors(pindex, triang):
    return triang.vertex_neighbor_vertices[1][triang.vertex_neighbor_vertices[0][pindex]:triang.vertex_neighbor_vertices[0][pindex+1]]

def find_all_neighbors(tri):
    neig = []
    for ii in range(len(tri.points)):
        neig.append(find_neighbors(ii,tri))
    return neig

def find_number_neighbors(tri):
    Nneighbours = []
    for ii in range(len(tri.points)):
        Nneighbours.append(len(find_neighbors(ii,tri)))
    return Nneighbours

def find_fraq_particle(tri):
    fraq = []
    triList = tri.simplices
    NneighboursList = find_number_neighbors(tri)
    for ii in range(len(tri.simplices)):
        ans = 0
        for jj in range(len(triList[ii])):
            ans = ans + 1/NneighboursList[triList[ii][jj]]
        fraq.append(ans)
    return fraq
        
def points_per_triangle(tri):
    return find_fraq_particle(tri)
#========================================================================================
# nouvelle version des fonction qui était en N^2
#========================================================================================

def center_triangle_2d(x1,x2,x3):
    Xtab = [x1,x2,x3]
    xCenter = 1/3*(Xtab[0][0]+Xtab[1][0]+Xtab[2][0])
    yCenter = 1/3*(Xtab[0][1]+Xtab[1][1]+Xtab[2][1])
    return xCenter,yCenter

def vec_float_to_str(vec):
    strVec = []
    for ii in range(len(vec)):
        strVec.append(f'{vec[ii]:.2f}')
    return strVec

def delaunay_density2d(tri):
    points = tri.points
    NpointsFraq = points_per_triangle(tri)
    density = []
    for ii in range(len(tri.simplices)):
        x1 = points[tri.simplices[ii]][0]
        x2 = points[tri.simplices[ii]][1]
        x3 = points[tri.simplices[ii]][2]
        density.append(NpointsFraq[ii]*1/triangleArea(x1,x2,x3))
    return density

def Axs_Delaunay_Voronoi_density(axs,tri,vor,log_scale=True,cmap=cm.Blues_r,display_point=True,display_vertice=True,twodCond=False):
    if twodCond:
        densityDelaunay = delaunay_density2d(tri)
    else:
        densityDelaunay = delaunay_density2d(tri) # ne pas oublier de change sa

    densityVoronoi = voronoi_volumes(vor)
    for ii in range(len(densityVoronoi)):
        densityVoronoi[ii] = 1/densityVoronoi[ii]

    if log_scale:
        for ii in range(len(densityDelaunay)):
            densityDelaunay[ii] = np.log(densityDelaunay[ii])
        for jj in range(len(densityVoronoi)):
            densityVoronoi[jj] = np.log(densityVoronoi[jj])


    minima = min(densityDelaunay)
    maxima = max(densityDelaunay)

    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    points = tri.points
    # draw Delaunay
    axs[0].triplot(points[:,0], points[:,1], tri.simplices)
    if display_point:
        axs[0].plot(points[:,0], points[:,1], 'o',markersize=1)
    for rr in range(len(tri.simplices)):
        polygon = plt.Polygon(points[tri.simplices[rr]],color=mapper.to_rgba(densityDelaunay[rr]))
        axs[0].add_patch(polygon)
    
    #draw Voronoi
    lineWidth = 2.0
    if not display_vertice:
        lineWidth=0.0
    voronoi_plot_2d(vor, ax=axs[1], show_points=display_point, show_vertices=display_vertice, line_colors='blue',line_width=lineWidth, line_alpha=0.6, point_size=2)
    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            axs[1].fill(*zip(*polygon), color=mapper.to_rgba(densityVoronoi[r]))


def SavePoints(fileName,Points):
    with open(fileName, 'wb') as f:
        np.save(f,Points)

def LoadPoints(fileName):
    with open(fileName, 'rb') as f:
        Points= np.load(f)
    return Points

##############################################################################################
# Fonction crée pour géré ce qui se passe en 6 dimension et les changement de format
##############################################################################################

def change_format_points_2d(nbPoints,axs1=0,axs2=1):
    pointsX = []
    pointsV = []
    for ii in range(len(nbPoints.pos[:,0])):
        pointsX.append([nbPoints.pos[ii][axs1],nbPoints.pos[ii][axs2]])
        pointsV.append([nbPoints.vel[ii][axs1],nbPoints.vel[ii][axs2]])
    return pointsX,pointsV

def change_format_points_6d(nb,scale_pos=1.0,scale_vel=1.0):
    points = []
    for ii in range(len(nb.pos[:,0])):
        points.append([1/scale_pos*nb.pos[ii][0],1/scale_pos*nb.pos[ii][1],1/scale_pos*nb.pos[ii][2],1/scale_vel*nb.vel[ii][0],1/scale_vel*nb.vel[ii][1],1/scale_vel*nb.vel[ii][2]])
    return points


