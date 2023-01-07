"""
@author: Karan Jakhar
Collection of functions specific to the project
"""
import numpy as np
import h5py
import math
from eqsdiscovery.derivativeLibrary import derivativeLibrary

def loadmat (DATA_DIR, key, Tarr=None, transposeData=False):
    """
    Loads a variable from .mat file as numpy array
    """ 
    # DATA_DIR (string): Directory of .mat file '/application/.../test.mat'
    # key (string): key / name of variable to be loaded
    # Tarr (n x 1 vector): index of elements to be loaded  
    # transposeData (boolean) : Take transpose of data
    
    data = h5py.File(DATA_DIR, 'r')
    if (transposeData == False):
        A = np.array(data[key])
    else:
        A = np.transpose(data[key])
    if (Tarr == None):
        return A
    else:
        B = A[:,:,Tarr]
        return B
    
def initializeGrid(N):
    # Generate for 2D grid size of NxN
    pi = math.pi
    Lx = 2*pi
    Ly = 2*pi

    dx = float(Lx)/N
    dy = float(Ly)/N

    x = np.arange(1-int(N/2),1+int(N/2))*dx
    y = np.arange(1-int(N/2),1+int(N/2))*dy

    kx = np.zeros(N)
    ky = np.zeros(N)

    kx[0:int(N/2)+1]=2*pi*np.arange(0,N/2+1)/float(Lx)
    kx[int(N/2)+1:N]=-kx[int(N/2)-1:0:-1]

    ky[0:int(N/2)+1]=2*pi*np.arange(0,N/2+1)/float(Ly)
    ky[int(N/2)+1:N]=-ky[int(N/2)-1:0:-1]

    N = len(ky)
    
    return x, y, kx, ky
    
def corr2(a,b):
    # Correlation coefficient of N x N x T array
    a = a - np.mean(a)
    b = b - np.mean(b)

    r = (a*b).sum() / np.sqrt((a*a).sum() * (b*b).sum());
    
    return r

def RMSE(a,b):
    # Root mean square error (absolute error)
        
    a_2norm = np.linalg.norm(a.flatten(), ord=2)
    ab_2norm = np.linalg.norm(a.flatten()-b.flatten(), ord=2)
    
    RMSE = ab_2norm/a_2norm
    
    return RMSE

def energyTransferS(theta, S1, S2, S3):
    # Energy transfer (P_S) using stress terms
    
    P = -np.multiply((S1-S3),theta['Ux'])-np.multiply(S2,theta['Uy']+theta['Vx'])    
        
    return P


def interScaleTransferPIvorS(Vor, S1, S2, S3):
    # Interscale transfer (T) using stress terms
    
    N = S1.shape[0]
    T = S1.shape[2]
    x, y, kx, ky = initializeGrid(N)
    
    thetaS1, thetaS1Name = derivativeLibrary(S1, kx, ky, T, orderMax=2, degree=1, fname='S1', returnf=False)
    thetaS2, thetaS2Name = derivativeLibrary(S2, kx, ky, T, orderMax=2, degree=1, fname='S2', returnf=False)
    thetaS3, thetaS3Name = derivativeLibrary(S3, kx, ky, T, orderMax=2, degree=1, fname='S3', returnf=False)
    thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, T, orderMax=2, degree=1, fname='Vor', returnf=False)
    
    PIvor = -thetaS1[:,:,:,5] + thetaS3[:,:,:,5] + thetaS2[:,:,:,3] - thetaS2[:,:,:,4]
    
    T = np.sign(thetaVor[:,:,:,3]+thetaVor[:,:,:,4]) * PIvor
    
    return T


def interScaleTransferPIuS(theta, S1, S2, S3):
    # Interrscale transfer (Tu) using stress terms
    
    N = S1.shape[0]
    T = S1.shape[2]
    x, y, kx, ky = initializeGrid(N)
    
    thetaS1, thetaS1Name = derivativeLibrary(S1, kx, ky, T, orderMax=2, degree=1, fname='S1', returnf=False)
    thetaS2, thetaS2Name = derivativeLibrary(S2, kx, ky, T, orderMax=2, degree=1, fname='S2', returnf=False)
    
    PIu = thetaS1[:,:,:,1] + thetaS2[:,:,:,2]
    T = np.sign(theta['Uxx']+theta['Uyy']) * PIu
        
    return T

def interScaleTransferPIvS(theta, S1, S2, S3):
    # Interscale transfer (Tv) using stress terms
    
    N = S1.shape[0]
    T = S1.shape[2]
    x, y, kx, ky = initializeGrid(N)
    
    thetaS2, thetaS2Name = derivativeLibrary(S2, kx, ky, T, orderMax=2, degree=1, fname='S2', returnf=False)
    thetaS3, thetaS3Name = derivativeLibrary(S3, kx, ky, T, orderMax=2, degree=1, fname='S3', returnf=False)
    
    PIv = thetaS2[:,:,:,1] + thetaS3[:,:,:,2]
    T = np.sign(theta['Vxx']+theta['Vyy']) * PIv
        
    return T

def interScaleTransferPIvorSvor(Vor, Svor1, Svor2):
    # Enstrophy transfer (Z) using stress terms (vorticity)
    
    N = Svor1.shape[0]
    T = Svor1.shape[2]
    x, y, kx, ky = initializeGrid(N)
    
    thetaSvor1, thetaSvor1Name = derivativeLibrary(Svor1, kx, ky, T, orderMax=1, degree=1, fname='Svor1', returnf=False)
    thetaSvor2, thetaSvor2Name = derivativeLibrary(Svor2, kx, ky, T, orderMax=1, degree=1, fname='Svor2', returnf=False)
    thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, T, orderMax=2, degree=1, fname='Vor', returnf=False)
    
    PIvor = thetaSvor1[:,:,:,1] + thetaSvor2[:,:,:,2]
    T = np.sign(thetaVor[:,:,:,3]+thetaVor[:,:,:,4]) * PIvor
    
    return T

def enstrophyTransferSvor(Vor, Svor1, Svor2):
    # Enstrophy transfer (Z) using stress terms (vorticity)
    
    N = Svor1.shape[0]
    T = Svor1.shape[2]
    x, y, kx, ky = initializeGrid(N)
    
    thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, T, orderMax=2, degree=1, fname='Vor', returnf=False)
    
    Z = -(Svor1 * thetaVor[:,:,:,1] + Svor2 * thetaVor[:,:,:,2])
    
    return Z

