"""
@author: Karan Jakhar
Equation Discovery (SGS term) for 2D-FHIT and RBC
"""

"""
Requirements:
U,V, [M,N,T] data to form library of basis function
PI/Stress [M,N,T] data

Input arguments
Threshold alpha : 100 (if not known)
Tolerance: 1e-6 (if not known)
"""

import sys
import numpy as np
import scipy.io as sio
from scipy.fftpack import fft, ifft, fft2, ifft2
import math
import h5py
import os

from eqsdiscovery.rvm import RVR
from eqsdiscovery.derivative import derivative
from eqsdiscovery.derivativeLibrary import derivativeLibrary
from eqsdiscovery.denormalizeWeight import denormalizeWeight
from eqsdiscovery.util import loadmat, corr2, initializeGrid

# Input Arguments system yName filterType N Nfilter libraryDerivativeOrder thresholdAlpha tolerance nosSnapshot skipSnapshot subtractModel SAVE_DIR

# python 2DTurbulenceUV.py 2D_FHIT S Re20kNX1024nx4ny0r0p1 train1 gaussian 32 default 4 50 1e-6 1 1 0 ./
# python 2DTurbulenceUV.py RBC J NX2048NY400Ra1E6Pr100 train1 gaussian 128 default 4 50 1e-6 1 1 0 ./

system = sys.argv[1] # RBC 2D_FHIT

N = int(sys.argv[6]) # 8 16 32 64 128 1024 2048
Nfilter = sys.argv[7] # default (Nfilter = N) or 8, ...

dataType = sys.argv[3] 
# 2D_FHIT: Re20kNX1024nx4ny0r0p1 Re20kNX1024nx25ny25r0p1 Re100kNX2048nx4ny0r0p1
# RBC: NX2048NY400Ra1E6Pr100 NX2048NY400Ra4E7Pr100 NX2048NY400Ra4E7Pr7

yName = sys.argv[2] 

ICtype = sys.argv[4]  # train1, train2, ... test1, val
filterType = sys.argv[5] # gaussian box gaussian+box
dealias = 0
subtractModel = sys.argv[13] # Subtract gradient model O(\Delta^2) from residual stress before discovery

orderMax = 8
degree=1

pruneV = True  # True, False # Removing derivatives of v equal to derivatives of u (continuity equation)
galileanInvariance = True  # Removes galiliean invariant terms from the library U, V 
libraryDerivativeOrder = int(sys.argv[8]) # Derivatives of order =< libraryDerivativeOrder will be

use_h5py = 1; # use_h5py = 1 will use h5py.File to load mat, use_h5py = 0 will use sio.loadmat

# save_RESULT_DIR = 0 ; # =1 if want to save results in RESULT_DIR, else results will be saved in ./

# y = \Phi * w

# Sampling every skipSnapshot th snapshot of data from the total data
skipSnapshot = int(sys.argv[12])

# # of snapshots of data to be analyzed. Will use first nosSnapshot snapshots from the sampled data
nosSnapshot = int(sys.argv[11])

#RVM hyperparameters
# thresholdAlpha = 100
# tolerance = 1e-6
thresholdAlpha = float(sys.argv[9]) # 1 10 100 1000 100000 ...
tolerance = float(sys.argv[10]) # 1e-6 1e-4 ...

# Directory to save results
RESULT_DIR = sys.argv[14] + dataType + '/'
    
x1Name = 'U'
x2Name = 'V'

# Scaling snapshots / grid points for each grid size (N) to keep equal grid points across grid sizes (N)
scaleSnapshots = False
# Partially selecting grid points for N=1024 to keep number of grid points consistent across other Ns
# Uniformly 128 x 128 grid points for N=1024 for each snapshot
# Selecting 4 times the snapshot for N=64 to keep # grid points consistent with N=128
# Selecting 16 times the snapshot for N=32 to keep # grid points consistent with N=128

if scaleSnapshots:
    
    temp = (N*N)/(128*128)
    nosSnapshot = int(np.floor(nosSnapshot/temp))
    
    if nosSnapshot < 1:
        skipSnapshot = 1
        nosSnapshot = 1
        print('Note: nosSnapshot modified to = 1')
    else:
        if N == 2048:
            skipSnapshot = int(np.floor((100/nosSnapshot)))
        elif N == 1024:
            skipSnapshot = int(np.floor((500/nosSnapshot)))
        else:
            skipSnapshot = int(np.floor(skipSnapshot*temp))
    if skipSnapshot == 0:
        skipSnapshot = 1
        nosSnapshot = 10000
        print('Note: nosSnapshot = 10000 since skipSnapshot <1, modified skipsnapShot = ', skipSnapshot)
        
if yName == 'Svor' or yName == 'Svor1' or yName == 'Svor2':
    libraryTerms = 'UVVor'
elif yName == 'J' or yName == 'J1' or yName == 'J2':
    libraryTerms = 'UVT'
else:
    libraryTerms = 'UV' 
# 'UV' : Library containing derivatives of U, V, and, pairs of derivatives of U,V.
# 'UVVor' : Library containing derivatives of U, V, Vor, and, pairs of derivatives of U,Vor and V,Vor

if libraryTerms == 'UV':
    x3Name = ''
if libraryTerms == 'UVVor':
    x3Name = 'Vor'
if libraryTerms == 'UVT':
    x3Name = 'T'
    
if Nfilter == 'default':
    DATA_DIR = '../data/' + system + '/' + dataType + '/' + filterType + '/' + ICtype + '/NX' + str(N) + '/'

else:
    Nfilter = int(Nfilter)
    DATA_DIR = '../data/' + system + '/' + dataType + '/' + filterType + '/' + ICtype + '/NX' + str(N) + '/' + 'Nfilter' + \
    str(Nfilter) + '/'
    
if dealias == 1:
    DATA_DIR = DATA_DIR + 'dealias/'
    
########################
if system == 'RBC':
    import matlab.engine
    eng = matlab.engine.start_matlab()
    scaleSnapshots = False
    M = 400
else:
    M = N
#######################

Norg = N;
    
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

# Filter size

if system == '2D_FHIT':
    if Nfilter == 'default':
        Delta = 2*Lx/N
    else:
        Delta = 2*Lx/Nfilter
else:
    if Nfilter == 'default':
        Delta = 2*(6*pi)/N
    else:
        Delta = 2*(6*pi)/Nfilter
    
print('------------------------------------------------------------')
print('data                          =  ',dataType)
print('ICtype                        =  ',ICtype)
print('filterType                    =  ',filterType)
print('yName                         =  ',yName)
print('Library Terms                 =  ',libraryTerms)
print('Ngrid                         =  ',N)
print('Nfilter                       =  ',Nfilter)
print('libraryDerivativeOrder        =  ',libraryDerivativeOrder)
print('Threshold Alpha               =  ',thresholdAlpha)
print('Tolerance                     =  ',tolerance)
print('libraryDerivativeOrder        =  ',libraryDerivativeOrder)
print('Snapshots                     =  ',nosSnapshot)
print('skipSnapshot                  =  ',skipSnapshot)
print('Subtract Model                =  ',subtractModel)
print('dealiased data                =  ',dealias)
print('Filter Size                   =  ',Delta)
print('Data Directory                =  ',DATA_DIR)
print('Result Directory              =  ',RESULT_DIR)
print('------------------------------------------------------------')

# .mat files are column ordered while python files are row ordered
# Take a transpose 
if use_h5py == 1:
    if libraryTerms == 'UV' or libraryTerms == 'UVVor' or libraryTerms == 'UVT':
        pp = h5py.File(DATA_DIR+'UV.mat', 'r')
        Ua = np.transpose(pp['U']);
        Va = np.transpose(pp['V']);
    if libraryTerms == 'UVVor':
        pp = h5py.File(DATA_DIR+'PsiVor.mat', 'r')
        Vora = np.transpose(pp['Vor']);
    if libraryTerms == 'UVT':
        pp = h5py.File(DATA_DIR+'T.mat', 'r')
        Ta = np.transpose(pp['T']);
elif use_h5py == 0:
    if libraryTerms == 'UV' or libraryTerms == 'UVVor' or libraryTerms == 'UVT':
        pp = sio.loadmat(DATA_DIR+'UV.mat');
        Ua = pp['U'];
        Va = pp['V'];
    if libraryTerms == 'UVVor':
        pp = sio.loadmat(DATA_DIR+'PsiVor.mat');
        Vora = pp['Vor'];
    if libraryTerms == 'UVT':
        pp = sio.loadmat(DATA_DIR+'T.mat');
        Ta = pp['T'];
del pp;

if libraryTerms == 'UV' or libraryTerms == 'UVVor' or libraryTerms == 'UVT':
    
    if libraryTerms == 'UVVor':
        # Sampling first T snapshots of the data
        Vor = Vora[:,:,0:nosSnapshot]
        X3 = Vor
        del Vora
        
    if libraryTerms == 'UVT':
        # Sampling first T snapshots of the data
        T = Ta[:,:,0:nosSnapshot]
        X3 = T
        del Ta

    # Sampling first T snapshots of the data
    U = Ua[:,:,0:nosSnapshot]
    V = Va[:,:,0:nosSnapshot]
    
    del Ua
    del Va
    
    if system == 'RBC':
        
        thetaU, thetaUName = eng.derivativeLibrary_RBC(
            matlab.double(np.ndarray.tolist(U)), nosSnapshot, orderMax, degree, x1Name, nargout=2)
        thetaV, thetaVName = eng.derivativeLibrary_RBC(
            matlab.double(np.ndarray.tolist(V)), nosSnapshot, orderMax, degree, x2Name, nargout=2)
   
        thetaU = np.array(thetaU);
        thetaV = np.array(thetaV);
        thetaUName = np.array(thetaUName);
        thetaVName = np.array(thetaVName);
    
    else:
    
        thetaU, thetaUName = derivativeLibrary(U, kx, ky, nosSnapshot, orderMax=orderMax, degree=degree, fname=x1Name, returnf=not(galileanInvariance))
        thetaV, thetaVName = derivativeLibrary(V, kx, ky, nosSnapshot, orderMax=orderMax, degree=degree, fname=x2Name, returnf=not(galileanInvariance))

    # Removing white space
    thetaUName = [word.strip() for word in thetaUName]
    thetaVName = [word.strip() for word in thetaVName]

    # Pruning v derivates which are equal to u derivatives (continuity equation) i.e. all dv/dy derivatives
    if pruneV:
        pruneVIndex = np.zeros(0, dtype=int)
        for ind in range(len(thetaVName)):
            for char in thetaVName[ind]:
                if char == 'y':
                    pruneVIndex = np.append(pruneVIndex,ind)
                    break

        thetaVName = np.delete(thetaVName,pruneVIndex)
        thetaV = np.delete(thetaV,pruneVIndex, axis=3)

    print('thetaUName = ')
    print(thetaUName)

    print('thetaVName = ')
    print(thetaVName)
    
    if libraryTerms == 'UVVor' or libraryTerms == 'UVT':
        
        if system == 'RBC':
            
            thetaVor, thetaVorName = eng.derivativeLibrary_RBC(
                matlab.double(np.ndarray.tolist(X3)), nosSnapshot, orderMax, degree, x3Name, nargout=2)
            thetaVor = np.array(thetaVor);
            thetaVorName = np.array(thetaVorName);
            
        else:
    
            # Calculating derivatives
            thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, nosSnapshot, orderMax=orderMax, degree=degree, fname=x3Name, returnf=True)

        # Removing white space
        thetaVorName = [word.strip() for word in thetaVorName]
        
        if x3Name == 'Vor':
            print('thetaVorName = ')
        elif x3Name == 'T':
            print('thetaTName = ')
            
        print(thetaVorName)
        

print('******************************************************************************')
print('Derivatives Calculated')
print('******************************************************************************')

# Making library of basis functions 

if libraryTerms == 'UV' or libraryTerms == 'UVVor' or libraryTerms == 'UVT':
    shapeU = thetaU.shape
    shapeV = thetaV.shape

if libraryTerms == 'UVVor' or libraryTerms == 'UVT':
    shapeVor = thetaVor.shape

temp = np.zeros((M,N,nosSnapshot,1))

# Making Library containing pairs of U and V derivatives
if libraryTerms == 'UV':
    
    # Pairing derivatives of U (UxUxx). Also includes non-paired derivatives (Ux)
    counterU = 0
    for countU1 in range(0,shapeU[-1]):
        for countU2 in range(countU1,shapeU[-1]):
            temp[:,:,:,0] = np.multiply(thetaU[:,:,:,countU1],thetaU[:,:,:,countU2])
            tempStr = thetaUName[countU1] + thetaUName[countU2]

            if (len(tempStr)-2) <= libraryDerivativeOrder:
                if (counterU == 0):
                    thetaU2 = temp
                    thetaUName2 = tempStr
                else: 
                    thetaU2 = np.concatenate((thetaU2, temp), axis=3)
                    thetaUName2 = np.array(np.append(thetaUName2, tempStr))
                counterU = counterU + 1
                
    # Removing Duplicate terms    
    a,uniqInd,c,d = np.unique(thetaUName2, return_index=True, return_inverse=True, return_counts=True)
    uniqInd = np.sort(uniqInd)
    thetaUName2 = np.array([thetaUName2[ind] for ind in uniqInd])
    thetaU2 = np.array([thetaU2[:,:,:,ind] for ind in uniqInd])
    # axis are moved in the above operation, getting the axis back in original order
    thetaU2 = np.moveaxis(thetaU2, 0, -1) 
    
    # Pairing derivatives of V (VxVxx). Also includes non-paired derivatives (Vx)
    counterV = 0
    for countV1 in range(0,shapeV[-1]):
        for countV2 in range(countV1,shapeV[-1]):
            temp[:,:,:,0] = np.multiply(thetaV[:,:,:,countV1],thetaV[:,:,:,countV2])
            tempStr = thetaVName[countV1] + thetaVName[countV2]

            if (len(tempStr)-2) <= libraryDerivativeOrder:
                if (counterV == 0):
                    thetaV2 = temp
                    thetaVName2 = tempStr
                else: 
                    thetaV2 = np.concatenate((thetaV2, temp), axis=3)
                    thetaVName2 = np.array(np.append(thetaVName2, tempStr))
                counterV = counterV + 1
                    
    # Removing Duplicate terms    
    a,uniqInd,c,d = np.unique(thetaVName2, return_index=True, return_inverse=True, return_counts=True)
    uniqInd = np.sort(uniqInd)
    thetaVName2 = np.array([thetaVName2[ind] for ind in uniqInd])
    thetaV2 = np.array([thetaV2[:,:,:,ind] for ind in uniqInd])
    # axis are moved in the above operation, getting the axis back in original order
    thetaV2 = np.moveaxis(thetaV2, 0, -1) 
                          
    # Pairing derivatives of U and V (UxVx).      
    counterUV = 0
    thetaUV = temp
    thetaUVName = np.array('tempstr')
    for countU in range(1,shapeU[-1]):
        for countV in range(1,shapeV[-1]):
            temp[:,:,:,0] = np.multiply(thetaU[:,:,:,countU],thetaV[:,:,:,countV])
            tempStr = thetaUName[countU] + thetaVName[countV]

            if (len(tempStr)-2) <= libraryDerivativeOrder:
                
                thetaUV = np.concatenate((thetaUV, temp), axis=3)
                thetaUVName = np.array(np.append(thetaUVName, tempStr))
                counterUV = counterUV + 1

    # Removing data including constants
    thetaU2 = np.delete(thetaU2, obj=0, axis=3)
    thetaV2 = np.delete(thetaV2, obj=0, axis=3)
    thetaUV = np.delete(thetaUV, obj=0, axis=3)
    thetaUName2 = np.delete(thetaUName2, obj=0)
    thetaVName2 = np.delete(thetaVName2, obj=0)
    thetaUVName = np.delete(thetaUVName, obj=0)

    shapeU2 = thetaU2.shape
    shapeV2 = thetaV2.shape

    # Adding
    theta = np.concatenate((thetaU2,thetaV2,thetaUV),axis=3)
    thetaAllName = np.concatenate((thetaUName2,thetaVName2,thetaUVName),axis=0)
    
    print('Shape of array of all derivatives                = ', theta.shape)
    print('List of all derivatives                          = ', thetaAllName)
    
    del thetaU2
    del thetaV2
    del thetaUV
    
# Making Library containing pairs of U Vor and V Vor derivatives
if libraryTerms == 'UVVor' or libraryTerms == 'UVT':
    
    if libraryTerms == 'UVVor':
        tempsub = 4
    if libraryTerms == 'UVT':
        tempsub = 2

    # Pairing derivatives of U and Vor (UxVor). 
    counterUVor = 0
    thetaUVor = temp
    thetaUVorName = 'tempstr'
    for countUVor1 in range(1,shapeU[-1]):
        for countUVor2 in range(1,shapeVor[-1]):
            temp[:,:,:,0] = np.multiply(thetaU[:,:,:,countUVor1],thetaVor[:,:,:,countUVor2])
            tempStr = thetaUName[countUVor1] + thetaVorName[countUVor2]

            if (len(tempStr)-tempsub) <= libraryDerivativeOrder:
                
                thetaUVor = np.concatenate((thetaUVor, temp), axis=3)
                thetaUVorName = np.array(np.append(thetaUVorName, thetaUName[countUVor1] + thetaVorName[countUVor2]))
                counterUVor = counterUVor + 1

    counterVVor = 0
    thetaVVor = temp
    thetaVVorName = 'tempstr'
    for countVVor1 in range(1,shapeV[-1]):
        for countVVor2 in range(1,shapeVor[-1]):
            temp[:,:,:,0] = np.multiply(thetaV[:,:,:,countVVor1],thetaVor[:,:,:,countVVor2])
            tempStr = thetaVName[countVVor1] + thetaVorName[countVVor2]

            if (len(tempStr)-tempsub) <= libraryDerivativeOrder:
                thetaVVor = np.concatenate((thetaVVor, temp), axis=3)
                thetaVVorName = np.array(np.append(thetaVVorName, thetaVName[countVVor1] + thetaVorName[countVVor2]))
                counterVVor = counterVVor + 1  
                
    tempU = np.zeros(1,dtype=int)
    
    for countU in range(0,len(thetaUName)):
        if (len(thetaUName[countU])-2) >= libraryDerivativeOrder:
            tempU = np.append(tempU,countU)
            
    tempV = np.zeros(1,dtype=int)
    for countV in range(0,len(thetaVName)):
        if (len(thetaVName[countV])-2) >= libraryDerivativeOrder:
            tempV = np.append(tempV,countV)  
            
    tempVor = np.zeros(1,dtype=int)
    for countVor in range(0,len(thetaVorName)):
        if (len(thetaVorName[countVor])-tempsub) >= libraryDerivativeOrder:
            tempVor = np.append(tempVor,countVor)  
    
    # Removing First term including garbage value
    thetaUVor = np.delete(thetaUVor, obj=0, axis=3)
    thetaVVor = np.delete(thetaVVor, obj=0, axis=3)
    thetaUVorName = np.delete(thetaUVorName, obj=0)
    thetaVVorName = np.delete(thetaVVorName, obj=0)
            
    # Removing higher order data including constants data including constants
    thetaU = np.delete(thetaU,obj=tempU,axis=3)
    thetaUName = np.delete(thetaUName,obj=tempU)
    thetaV = np.delete(thetaV,obj=tempV,axis=3)
    thetaVName = np.delete(thetaVName,obj=tempV)
    thetaVor = np.delete(thetaVor,obj=tempVor,axis=3)
    thetaVorName = np.delete(thetaVorName,obj=tempVor)
    
    print(thetaUName)
    print(thetaVName)
    print(thetaVorName)
    print(thetaUVorName)
    print(thetaVVorName)

    # Adding
    theta = np.concatenate((thetaU,thetaV,thetaVor,thetaUVor,thetaVVor),axis=3)
    thetaAllName = np.concatenate((thetaUName,thetaVName,thetaVorName,thetaUVorName,thetaVVorName),axis=0)
    
    del thetaU
    del thetaV
    del thetaVor
    del thetaUVor
    del thetaVVor
    
print('Library Formed of # ' + str(len(thetaAllName)) + ' terms')
print('Terms of the library thetaAllName = ')
print(thetaAllName)

shapeTheta = theta.shape
theta = theta.reshape(np.prod(shapeTheta[0:3]),shapeTheta[-1])

thetaAllMean = theta.mean(axis=0)
thetaAllStd = theta.std(axis=0)

for count in range(0,theta.shape[1]):
    theta[:,count] = (theta[:,count]-thetaAllMean[count])/thetaAllStd[count]
    
print("Library of Functions Created\n")
print('Shape of Library (theta) =' , thetaAllMean.shape)

# Checking for nan, inf in theta
if (np.sum(np.isnan(theta)) > 0):
    print('NaN value encountered in theta. Exiting code')
    print('# NaN value in theta = ', np.sum(np.isnan(theta)))
    sys.exit()
       
if (np.sum(np.isinf(theta)) > 0):
    print('Inf value encountered in theta. Exiting code')
    print('# Inf value in theta = ', np.sum(np.isinf(theta)))
    sys.exit()

# ya, dataTitle = denorm(DATA_DIR, yName, datatype, skipSnapshot)

if yName == 'S' or yName == 'S1' or yName == 'S2' or yName == 'S3':
    tempMatFilename = 'S.mat'
if yName == 'Svor' or yName == 'Svor1' or yName == 'Svor2':
    tempMatFilename = 'Svor.mat'
if yName == 'J' or yName == 'J1' or yName == 'J2':
    tempMatFilename = 'J.mat'
if yName == 'PIuv' or yName == 'PIu' or yName == 'PIv':
    tempMatFilename = 'PIuv.mat'
if yName == 'PIvor':
    tempMatFilename = 'PIvor.mat'
       
# If yName = S, the RVM will run on all S1, S2, S3 in consecutive order, similarly for Svor, PIuv
# if yName = S1, the RVM will run only on S1, similariliy for Svor1, PIu
       
tempyName = [yName]
if (yName == 'S'):
    tempyName = ['S1', 'S2', 'S3']
if (yName == 'Svor'):
    tempyName = ['Svor1', 'Svor2']
if (yName == 'J'):
    tempyName = ['J1', 'J2']
if (yName == 'PIuv'):
    tempyName = ['PIu', 'PIv']

for yName in tempyName:   
       
    if use_h5py == 1:
        yData = h5py.File(DATA_DIR+tempMatFilename, 'r')
        ya = np.transpose(yData[yName])
    elif use_h5py == 0:
        yData = sio.loadmat(DATA_DIR+tempMatFilename)
        ya = yData[yName];

    y = ya[:,:,0:nosSnapshot]
    
    del yData, ya
    
    ya = np.zeros((N,N,nosSnapshot))
#     if scaleSnapshots and int(sys.argv[3]) == 1024:
#         for countx in range(len(indx)): 
#             for county in range(len(indy)):
#                 ya[countx,county,:] = y[indx[countx],indy[county],:]
#         y = ya
#         del ya
        
    if subtractModel == 'GM2':

        thetaU, thetaUName = derivativeLibrary(U, kx, ky, nosSnapshot, orderMax=1, degree=1, fname=x1Name, returnf=False)
        thetaV, thetaVName = derivativeLibrary(V, kx, ky, nosSnapshot, orderMax=1, degree=1, fname=x2Name, returnf=False)
        
        if yName == 'Svor1' or yName == 'Svor2':
            thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, nosSnapshot, orderMax=1, degree=1, fname=x2Name, returnf=False)
        
        if yName == 'S1':
            yS1G2 = (Delta*Delta/12)*(thetaU[:,:,:,1]*thetaU[:,:,:,1]+thetaU[:,:,:,2]*thetaU[:,:,:,2])
            y = y - yS1G2
            del yS1G2
        if yName == 'S2':
            yS2G2 = (Delta*Delta/12)*(thetaU[:,:,:,1]*thetaV[:,:,:,1]+thetaU[:,:,:,2]*thetaV[:,:,:,2])
            y = y - yS2G2
            del yS2G2
        if yName == 'S3':
            yS3G2 = (Delta*Delta/12)*(thetaV[:,:,:,1]*thetaV[:,:,:,1]+thetaV[:,:,:,2]*thetaV[:,:,:,2])
            y = y - yS3G2
            del yS3G2

        if yName == 'Svor1':
            ySvor1G2 = (Delta*Delta/12)*(thetaU[:,:,:,1]*thetaVor[:,:,:,1]+thetaU[:,:,:,2]*thetaVor[:,:,:,2])
            y = y - ySvor1G2
            del ySvor1G2
            
        if yName == 'Svor2':
            ySvor2G2 = (Delta*Delta/12)*(thetaV[:,:,:,1]*thetaVor[:,:,:,1]+thetaV[:,:,:,2]*thetaVor[:,:,:,2])
            y = y - ySvor2G2
            del ySvor2G2
            
        del thetaU, thetaV
        
    if subtractModel == 'GM4':

        thetaU, thetaUName = derivativeLibrary(U, kx, ky, nosSnapshot, orderMax=3, degree=1, fname=x1Name, returnf=False)
        thetaV, thetaVName = derivativeLibrary(V, kx, ky, nosSnapshot, orderMax=3, degree=1, fname=x2Name, returnf=False)
        
        if yName == 'Svor1' or yName == 'Svor2':
            thetaVor, thetaVorName = derivativeLibrary(Vor, kx, ky, nosSnapshot, orderMax=3, degree=1, fname=x2Name, returnf=False)
        
        Delta2 = Delta*Delta
        Delta4 = Delta*Delta*Delta*Delta
                
        if yName == 'S1':
            yS1G4 = (Delta2/12)*(thetaU[:,:,:,1]*thetaU[:,:,:,1]+thetaU[:,:,:,2]*thetaU[:,:,:,2]) - (
                10*Delta4/1728*(thetaU[:,:,:,1]*thetaU[:,:,:,9]+thetaU[:,:,:,1]*thetaU[:,:,:,8])) + (
                Delta4/288*(thetaU[:,:,:,3]*thetaU[:,:,:,3] + thetaU[:,:,:,4]*thetaU[:,:,:,4])) - (
                10*Delta4/3456*(thetaU[:,:,:,3]*thetaU[:,:,:,4]))
            y = y - yS1G4
            del yS1G4
            
        if yName == 'S2':
            yS2G4 = (Delta2/12)*(thetaU[:,:,:,1]*thetaV[:,:,:,1]+thetaU[:,:,:,2]*thetaV[:,:,:,2]) - (
                5*Delta4/1728*(thetaU[:,:,:,1]*(thetaV[:,:,:,9]+thetaV[:,:,:,8]) + thetaV[:,:,:,1]*(
                thetaU[:,:,:,8]+thetaU[:,:,:,9]))) + (
                Delta4/288*(thetaU[:,:,:,3]*thetaV[:,:,:,3] + thetaU[:,:,:,4]*thetaV[:,:,:,4])) - (
                5*Delta4/3456*(thetaU[:,:,:,3]*thetaV[:,:,:,4] + thetaU[:,:,:,4]*thetaV[:,:,:,3]))
            y = y - yS2G4
            del yS2G4
            
        if yName == 'S3':
            yS3G4 = (Delta2/12)*(thetaV[:,:,:,1]*thetaV[:,:,:,1]+thetaV[:,:,:,2]*thetaV[:,:,:,2]) - (
                10*Delta4/1728*(thetaV[:,:,:,1]*thetaV[:,:,:,9]+thetaV[:,:,:,1]*thetaV[:,:,:,8])) + (
                Delta4/288*(thetaV[:,:,:,3]*thetaV[:,:,:,3] + thetaV[:,:,:,4]*thetaV[:,:,:,4])) - (
                10*Delta4/3456*(thetaV[:,:,:,3]*thetaV[:,:,:,4]))
            y = y - yS3G4
            del yS3G4
            
        if yName == 'Svor1':
            ySvor1G4 = (Delta2/12)*(thetaU[:,:,:,1]*thetaVor[:,:,:,1]+thetaU[:,:,:,2]*thetaVor[:,:,:,2]) - (
                5*Delta4/1728*(thetaU[:,:,:,1]*(thetaVor[:,:,:,9]+thetaVor[:,:,:,8]) + thetaVor[:,:,:,1]*(
                thetaU[:,:,:,8]+thetaU[:,:,:,9]))) + (
                Delta4/288*(thetaU[:,:,:,3]*thetaVor[:,:,:,3] + thetaU[:,:,:,4]*thetaVor[:,:,:,4])) - (
                5*Delta4/3456*(thetaU[:,:,:,3]*thetaVor[:,:,:,4] + thetaU[:,:,:,4]*thetaVor[:,:,:,3]))
            y = y - ySvor1G4
            del ySvor1G4
            
        if yName == 'Svor2':
            ySvor2G4 = (Delta2/12)*(thetaV[:,:,:,1]*thetaVor[:,:,:,1]+thetaV[:,:,:,2]*thetaVor[:,:,:,2]) - (
                5*Delta4/1728*(thetaV[:,:,:,1]*(thetaVor[:,:,:,9]+thetaVor[:,:,:,8]) + thetaVor[:,:,:,1]*(
                thetaV[:,:,:,8]+thetaV[:,:,:,9]))) + (
                Delta4/288*(thetaV[:,:,:,3]*thetaVor[:,:,:,3] + thetaV[:,:,:,4]*thetaVor[:,:,:,4])) - (
                5*Delta4/3456*(thetaV[:,:,:,3]*thetaVor[:,:,:,4] + thetaV[:,:,:,4]*thetaVor[:,:,:,3]))
            y = y - ySvor2G4
            del ySvor2G4
            
        del thetaU, thetaV

    yUnNorm = y.flatten(order='C')
    yMean = yUnNorm.mean(axis=0)
    yStd = yUnNorm.std(axis=0)
    del y
    y = (yUnNorm-yMean)/yStd
    del yUnNorm
    
    print(y.shape)
    print(theta.shape)
    print(thetaAllName.shape)

    """RVM"""

    print('----------------------------------------')
    print("Starting RVM on", yName)
    print('Threshold Alpha    =  ',str(thresholdAlpha) )
    print('Tolerance          =  ', tolerance)
    print('========================================')

    clf = RVR(threshold_alpha= thresholdAlpha, tol=tolerance, verbose=True, standardise=True)
    fitted = clf.fit(theta    , y     , thetaAllName        )
    scoreR2 = clf.score_R2(theta,y) 
    MSE     = clf.score_MSE(theta,y) 
    weightMean  = clf.m_#[0]
    alpha   = clf.alpha_

    # Eq 5 in the paper
    i_s = clf.beta_ * np.dot( clf.phi.T, clf.phi ) +  np.diag( clf.alpha_ )
    # Eq 5 in the paper
    sigma_ = np.linalg.inv(i_s)
    weightStd = np.sqrt(np.diag(sigma_))
    # Eq 6 in the paper| m_ is the mu in the paper
    #weightMean = clf.beta_ * np.dot( sigma_, np.dot( clf.phi.T, clf.t ) ) 
    # Error bar

    retainedInd = np.where(clf.retained_)[0]
    thetaName = [thetaAllName[v] for v in retainedInd]
    thetaMean = np.array([thetaAllMean[v] for v in retainedInd])
    thetaStd = np.array([thetaAllStd[v] for v in retainedInd])

    errorBar = np.zeros(len(sigma_), dtype=float)
    for count in range(0,len(sigma_)):
        errorBar[count] = (weightStd[count]/(weightMean[count]))**2

    # Scaling discovered weights

    weightMeanScaled, constant = denormalizeWeight(weightMean, yMean, yStd, thetaMean, thetaStd)

    weightMeanScaledDelta = np.zeros(len(thetaName))
    for count in range(0,len(thetaName)):
        tempU = str.count(thetaName[count],'U')
        tempV = str.count(thetaName[count],'V')
        
        if libraryTerms == 'UVVor':
            tempVor = str.count(thetaName[count],'Vor')
            tempV = tempV - tempVor
            orderTerm = len(thetaName[count]) - tempU - tempV - 3*tempVor
        elif libraryTerms == 'UVT':
            tempT = str.count(thetaName[count],'T')
            orderTerm = len(thetaName[count]) - tempU - tempV - tempT
        else:
            orderTerm = len(thetaName[count]) - tempU - tempV
            
        weightMeanScaledDelta[count] = weightMeanScaled[count]/(Delta**orderTerm)
        
    print('------------------------------------------------------------')
    print('Threshold Alpha        =  ', str(thresholdAlpha) )
    print('Tolerance              =  ', tolerance)
    print('Bases Retained         =  ', clf.n_retained       )
    print('Model R2 score         =  ', clf.score_R2(theta,y)  )
    print('Model MSE              =  ', clf.score_MSE(theta,y) )
    print('Error Bar              =  ', np.sum(errorBar))
    print('------------------------------------------------------------')
    print('Weights Normalized     =  ', weightMean)
    print('Std of bases           =  ', weightStd)
    print('Error Bar Total        =  ', errorBar)
    print('============================================================')
    print('Weights Scaled         =  ', weightMeanScaled)
    print('Delta                  =  ', Delta)
    print('Weights Scaled (Delta) =  ', weightMeanScaledDelta)
    print('============================================================')

    tempstr = yName + x1Name + x2Name + x3Name
    tempstr2 = '_' + 'NX' + str(Norg)
    
    if Nfilter == 'default':
        tempstr3 = ''
    else:
        tempstr3 = '_' + str(int(Nfilter))
            
    tempstr4 = '_O' + str(int(libraryDerivativeOrder))
    tempstr5 = '_T' + str(int(nosSnapshot))
    
    if (dealias == 1):
        tempstr6 = '_dealias'
    else:
        tempstr6 = ''
    
    if (subtractModel == '0'):
        tempstr7 = ''
    else:
        tempstr7 = '_' + subtractModel

    filename = tempstr + '_' + filterType + tempstr2 + tempstr3 + tempstr4 + '_' + str(int(thresholdAlpha)) + '_' + str(tolerance) + tempstr5 + tempstr6 + tempstr7
    
    try:
        os.mkdir(RESULT_DIR)
    except OSError as error:
        print(error)

    mdict = {'dataType':dataType, 'filterType':filterType, 'ICtype':ICtype, 'Ngrid':N,
             'thresholdAlpha':thresholdAlpha, 'tolerance':tolerance, 'Nfilter':Nfilter,
             'yName':yName, 'yMean':yMean, 'yStd':yStd, 'thetaMean':thetaMean,
             'thetaStd':thetaStd, 'weightMean':weightMean, 'weightStd':weightStd,
             'weightMeanScaled':weightMeanScaled, 'constant':constant,
             'weightMeanScaledDelta':weightMeanScaledDelta, 'Delta':Delta,
             'errorBar':errorBar, 'MSE':MSE, 'alpha':alpha, 'R2':scoreR2,
             'thetaName':thetaName, 'thetaAllName':thetaAllName, 'skipSnapshot':skipSnapshot,
             'libraryDerivativeOrder':libraryDerivativeOrder, 'libraryTerms':libraryTerms, 'noOfBases':clf.n_retained,
             'pruneV':pruneV, 'galileanInvariance':galileanInvariance, 'subtractModel':subtractModel, 'dealias':dealias,
             'retainedInd':retainedInd, 'nosSnapshot':nosSnapshot}
    sio.savemat(RESULT_DIR+filename+'.mat', mdict)
    # To load .mat file use sio.loadmat(file_name.mat, simplify_cells=True)
    # Each term thetaAllName and other strings would be saved as characters of equal length. Empty spaces would be used to make all the characters of equal lengths. Use the following to remove the empty space in these.
    # thetaAllName = [word.strip() for word in thetaAllName]

    np.savez(RESULT_DIR+filename , dataType=dataType, filterType=filterType, ICtype=ICtype, Ngrid=N,
             thresholdAlpha=thresholdAlpha, tolerance=tolerance, Nfilter=Nfilter,
             yName=yName, yMean=yMean, yStd=yStd, thetaMean=thetaMean,
             thetaStd=thetaStd, weightMean=weightMean, weightStd=weightStd,
             weightMeanScaled=weightMeanScaled, constant=constant,
             weightMeanScaledDelta=weightMeanScaledDelta, Delta=Delta,
             errorBar=errorBar, MSE=MSE, alpha=alpha, R2=scoreR2,
             thetaName = thetaName, thetaAllName=thetaAllName,
             libraryDerivativeOrder = libraryDerivativeOrder, libraryTerms=libraryTerms, noOfBases=clf.n_retained,
             pruneV=pruneV, galileanInvariance=galileanInvariance, subtractModel=subtractModel, dealias=dealias,
             retainedInd=retainedInd, nosSnapshot=nosSnapshot)
    print(filename)
