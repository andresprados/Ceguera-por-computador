import numpy as np
import misc
from numpy import linalg as la
from scipy import ndimage
from matplotlib import pyplot as ppl
import maxflow
from maxflow.fastmin import aexpansion_grid_step

def plothline_TODO(line, axes = None):
    """Plot a line given its homogeneous coordinates.
    
    Parameters
    ----------
    line : array_like
        Homogeneous coordinates of the line.
    axes : AxesSubplot
        Axes where the line should be plotted. If not given,
        line will be plotted in the active axis.
    """
    if axes == None:
        axes = ppl.gca()
    
    [x0, x1, y0, y1] = axes.axis()
    #     (x0, y0) ._____________________. (x1, y0)
    #              |                     |
    #              |                     |
    #              |                     |
    #              |                     |
    #              |                     |
    #              |                     |
    #     (x0, y1) .---------------------. (x1, y1)
    
    # TODO: Compute the intersection of the line with the image
    # borders.
    
    # TODO: Plot the line with axes.plot.
    #axes.plot(...)
    
    axes.axis([x0, x1, y0, y1])

def plot_epipolar_lines_TODO(image1, image2, F):
    """Ask for points in one image and draw the epipolar lines for those points.
    
    Parameters
    ----------
    image1 : array_like
        First image.
    image2 : array_like
        Second image.
    F : array_like
        3x3 fundamental matrix from image1 to image2.
    """
    # Prepare the two images.
    fig = ppl.gcf()
    fig.clf()
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(image1)
    ax1.axis('image')
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(image2)
    ax2.axis('image')
    ppl.draw()
    
    ax1.set_xlabel("Choose points in left image (or right click to end)")
    point = ppl.ginput(1, timeout=-1, show_clicks=False, mouse_pop=2, mouse_stop=3)
    while len(point) != 0:
        # point has the coordinates of the selected point in the first image.
        point = np.hstack([np.array(point[0]), 1])
        ax1.plot(point[0], point[1], '.r')
        
        # TODO: Determine the epipolar line.
        
        # TODO: Plot the epipolar line with plothline (the parameter 'axes' should be ax2).
        #plothline(..., axes=ax2)
        
        ppl.draw()
        # Ask for a new point.
        point = ppl.ginput(1, timeout=-1, show_clicks=False, mouse_pop=2, mouse_stop=3)
    
    ax1.set_xlabel('')
    ppl.draw()

def plot_correspondences_TODO(image1, image2, S, H1, H2):
    """
    Ask for points in the first image and plot their correspondences in
    the second image.
    
    Parameters
    ----------
    image1, image2 : array_like
        The images (before rectification)
    S : array_like
        The matrix of disparities.
    H1, H2 : array_like
        The homographies which rectify both images.
    """
    # Prepare the two images.
    fig = ppl.gcf()
    fig.clf()
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(image1)
    ax1.axis('image')
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(image2)
    ax2.axis('image')
    ppl.draw()
    
    ax1.set_xlabel("Choose points in left image (or right click to end)")
    point = ppl.ginput(1, timeout=-1, show_clicks=False, mouse_pop=2, mouse_stop=3)
    while len(point) != 0:
        # point has the coordinates of the selected point in the first image.
        point = np.c_[np.array(point), 1].T
        ax1.plot(point[0,:], point[1,:], '.r')
        
        # TODO: Determine the correspondence of 'point' in the second image.
        
        # TODO: Plot the correspondence with ax2.plot.
        #ax2.plot(...)
        
        ppl.draw()
        # Ask for a new point.
        point = ppl.ginput(1, timeout=-1, show_clicks=False, mouse_pop=2, mouse_stop=3)
    
    ax1.set_xlabel('')
    ppl.draw()

def main():
    pass

"""Cargado de las matrices de proyección"""
cameras=np.load("cameras.npz")
P1=cameras["left"]
P2=cameras["right"]

"""Función de Reconstrucción de puntos"""
"""Los puntos ingresados deben estar en coordenadas homogéneneas"""
"""Y las dimensiones de las matrices de puntos deben ser 3xN"""

def reconstruct(points1,points2,P1,P2):
    X_array=np.zeros(4)
    for i in range(len(points1.transpose())):
        M=np.zeros((4,4))
        for r in range(2):
            for s in range(4):
                ars=P1[r][s]-(P1[2][s])*(points1[r][i]/points1[2][i])
                M[r][s]=ars
        for r in range(2):
            for s in range(4):
                brs=P2[r][s]-(P2[2][s])*(points2[r][i]/points2[2][i])
                M[r+2][s]=brs
        U,S,V=la.svd(M, full_matrices=True)
        Vt=V.transpose()
        X=Vt[:,3]
        X=X/X[3]
        X_array=np.vstack((X_array,X))
    X_array=np.delete(X_array,0,0)
    X_array=X_array.transpose()
    return X_array

"""Función para el cálculo de Matriz Fundamental"""
def projmat2f(P1,P2):
    second_terma=np.zeros((3,3))
    A=P1[0:3,0:3]
    b=np.vstack(P1[:,3])
    B=P2[0:3,0:3]
    d=np.vstack(P2[:,3])
    first_term=(la.inv(B)).transpose()
    second_termb=np.dot(la.inv(B),d)-np.dot(la.inv(A),b)
    second_terma[1,0]=second_termb[2]
    second_terma[0,1]=-second_termb[2]
    second_terma[2,0]=-second_termb[1]
    second_terma[0,2]=second_termb[1]
    second_terma[2,1]=second_termb[0]
    second_terma[1,2]=-second_termb[0]
    final_term=la.inv(A)
    F=np.dot(first_term,np.dot(second_terma,final_term))
    return F

def f2projmat(F):
    #Cálculo de e' (epipolo de la segunda imagen)
    U,S,V=la.svd(F.transpose(), full_matrices=True)
    Vt=V.transpose()
    print Vt
    e_prima=Vt[:,2]
    print e_prima
    e_prima=e_prima/e_prima[2]
    print e_prima
    #Cálculo de la matriz skew de e'
    e_skew=np.zeros((3,3))
    e_skew[1,0]=e_prima[2]
    e_skew[0,1]=-e_prima[2]
    e_skew[2,0]=-e_prima[1]
    e_skew[0,2]=e_prima[1]
    e_skew[2,1]=e_prima[0]
    e_skew[1,2]=-e_prima[0]    
    #Cálculo de la primera matriz de proyección
    P1=np.concatenate((np.identity(3),np.zeros((3,1))),axis=1)
    z=np.dot(np.vstack(e_prima),(np.array([1,2,3])).reshape(1,3))
    print z
    #Cálculo de la segunda matriz de proyección
    P2=np.concatenate((np.dot(e_skew,F)+z,1000*np.vstack(e_prima)),axis=1)
    return P1,P2

"""Reconstrucción de la imágen, visualización 3D"""
img1=ndimage.imread("minoru_can_left.jpg")
img2=ndimage.imread("minoru_can_right.jpg")
points1,points2=misc.askpoints(img1,img2)
X=reconstruct(points1,points2,P1,P2)
ppl.figure(1)
ppl.subplot(311)
misc.plot3D(X[0,:],X[1,:],X[2,:])

"""Reproyección de puntos obtenidos en imagenes originales"""
new_points1=np.zeros((3,1))
new_points2=np.zeros((3,1))
for i in range(len(X.transpose())):
    coord1=np.vstack(np.dot(P1,X[:,i]))
    coord2=np.vstack(np.dot(P2,X[:,i]))
    new_points1=np.concatenate((new_points1,coord1),axis=1)
    new_points2=np.concatenate((new_points2,coord2),axis=1)
new_points1=np.delete(new_points1,0,1)
new_points2=np.delete(new_points2,0,1)
ppl.figure(2)
ppl.subplot(312)
misc.plothom(new_points1)
ppl.imshow(img1)
ppl.show()
ppl.figure(3)
ppl.subplot(313)
misc.plothom(new_points2)
ppl.imshow(img2)
ppl.show()

F=projmat2f(P1,P2)



    
    
    