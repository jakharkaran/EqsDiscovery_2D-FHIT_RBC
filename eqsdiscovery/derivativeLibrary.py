"""
@author: Karan Jakhar

Creates library of of basis functions upto degree n and order 8
"""

import sys
import numpy as np
import scipy.io as sio
from scipy.fftpack import fft, ifft, fft2, ifft2
import math
from eqsdiscovery.derivative import derivative

def derivativeLibrary(f, kx, ky, T, orderMax, degree, fname, returnf):
    
    # f: 2D matrix: containing function f
    # T: int: # of snapshots
    # orderMax: int: Maximum order of drivatives in the output
    # degree: int: Maximum degree of derivatives
    # fname: string: Name of function f
    # returnf: boolean: return function f along with its derivatives 'True: returns f. False: doesn't return f'
    # returnEmpty: boolean: return the first element (axis=3) as ones
    
    if orderMax > 8:
        print('**Error** Input derivative order is greater than 8. Derivatives upto order 8 are supported. ')
        sys.exit()
    
    M = len(kx)
    N = len(ky)
    
    if (orderMax >= 1):
        fx = np.zeros((M,N,T))
        fy = np.zeros((M,N,T))

        for t in range(0,T):
            fx[:,:,t],fy[:,:,t],fname1= derivative(f[:,:,t],kx,ky,order=1,fname=fname)
            
    if (orderMax >= 2):
        fxx = np.zeros((M,N,T))
        fyy = np.zeros((M,N,T))
        fxy = np.zeros((M,N,T))
        for t in range(0,T):
            fxx[:,:,t],fyy[:,:,t], fxy[:,:,t], fname2= derivative(f[:,:,t],kx,ky,order=2,fname=fname)
    
    if (orderMax >= 3):
        fxxx = np.zeros((M,N,T))
        fyyy = np.zeros((M,N,T))
        fxxy = np.zeros((M,N,T))
        fxyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxx[:,:,t], fyyy[:,:,t], fxxy[:,:,t], fxyy[:,:,t], fname3 = derivative(f[:,:,t],kx,ky,order=3,fname=fname)
            
    if (orderMax >= 4):
        fxxxx = np.zeros((M,N,T))
        fyyyy = np.zeros((M,N,T))
        fxxxy = np.zeros((M,N,T))
        fxxyy = np.zeros((M,N,T))
        fxyyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxxx[:,:,t], fyyyy[:,:,t], fxxxy[:,:,t], fxxyy[:,:,t], fxyyy[:,:,t], fname4 =\
            derivative(f[:,:,t],kx,ky,order=4,fname=fname)
            
    if (orderMax >= 5):
        fxxxxx = np.zeros((M,N,T))
        fyyyyy = np.zeros((M,N,T))
        fxxxxy = np.zeros((M,N,T))
        fxxxyy = np.zeros((M,N,T))
        fxxyyy = np.zeros((M,N,T))
        fxyyyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxxxx[:,:,t], fyyyyy[:,:,t], fxxxxy[:,:,t], fxxxyy[:,:,t], fxxyyy[:,:,t], fxyyyy[:,:,t], fname5 =\
            derivative(f[:,:,t],kx,ky,order=5,fname=fname)
            
    if (orderMax >= 6):
        fxxxxxx = np.zeros((M,N,T))
        fyyyyyy = np.zeros((M,N,T))
        fxxxxxy = np.zeros((M,N,T))
        fxxxxyy = np.zeros((M,N,T))
        fxxxyyy = np.zeros((M,N,T))
        fxxyyyy = np.zeros((M,N,T))
        fxyyyyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxxxxx[:,:,t], fyyyyyy[:,:,t], fxxxxxy[:,:,t], fxxxxyy[:,:,t], fxxxyyy[:,:,t], fxxyyyy[:,:,t], fxyyyyy[:,:,t], fname6 = derivative(f[:,:,t],kx,ky,order=6,fname=fname)            
            
    if (orderMax >= 7):
        fxxxxxxx = np.zeros((M,N,T))
        fyyyyyyy = np.zeros((M,N,T))
        fxxxxxxy = np.zeros((M,N,T))
        fxxxxxyy = np.zeros((M,N,T))
        fxxxxyyy = np.zeros((M,N,T))
        fxxxyyyy = np.zeros((M,N,T))
        fxxyyyyy = np.zeros((M,N,T))
        fxyyyyyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxxxxxx[:,:,t], fyyyyyyy[:,:,t], fxxxxxxy[:,:,t], fxxxxxyy[:,:,t], fxxxxyyy[:,:,t], fxxxyyyy[:,:,t], fxxyyyyy[:,:,t], fxyyyyyy[:,:,t], fname7 = derivative(f[:,:,t],kx,ky,order=7,fname=fname)
            
    if (orderMax >= 8):
        fxxxxxxxx = np.zeros((M,N,T))
        fyyyyyyyy = np.zeros((M,N,T))
        fxxxxxxxy = np.zeros((M,N,T))
        fxxxxxxyy = np.zeros((M,N,T))
        fxxxxxyyy = np.zeros((M,N,T))
        fxxxxyyyy = np.zeros((M,N,T))
        fxxxyyyyy = np.zeros((M,N,T))
        fxxyyyyyy = np.zeros((M,N,T))
        fxyyyyyyy = np.zeros((M,N,T))
        
        for t in range(0,T):
            fxxxxxxxx[:,:,t], fyyyyyyyy[:,:,t], fxxxxxxxy[:,:,t], fxxxxxxyy[:,:,t], fxxxxxyyy[:,:,t], fxxxxyyyy[:,:,t], fxxxyyyyy[:,:,t], fxxyyyyyy[:,:,t], fxyyyyyyy[:,:,t], fname8 = derivative(f[:,:,t],kx,ky,order=8,fname=fname)
            
        
    countf = 0
    
    if (orderMax == 0):
        farray = np.zeros((M,N,T,(1)+1))
        fnameArray = [None] * (1+1)
    
    if (orderMax == 1):
        farray = np.zeros((M,N,T,(3)+1))
        fnameArray = [None] * (3+1)
    
    if (orderMax == 2):
        farray = np.zeros((M,N,T,(6)+1))
        fnameArray = [None] * (6+1)
    
    if (orderMax == 3):
        farray = np.zeros((M,N,T,(10)+1))
        fnameArray = [None] * (10+1)
        
    if (orderMax == 4):
        farray = np.zeros((M,N,T,(15)+1))
        fnameArray = [None] * (15+1)
        
    if (orderMax == 5):
        farray = np.zeros((M,N,T,(21)+1))
        fnameArray = [None] * (21+1)
        
    if (orderMax == 6):
        farray = np.zeros((M,N,T,(28)+1))
        fnameArray = [None] * (28+1)
        
    if (orderMax == 7):
        farray = np.zeros((M,N,T,(36)+1))
        fnameArray = [None] * (36+1)
        
    if (orderMax == 8):
        farray = np.zeros((M,N,T,(45)+1))
        fnameArray = [None] * (45+1)
    
    if (degree >= 1):
        
        if (orderMax >= 0):
                      
            farray[:,:,:,countf] = np.ones(f.shape)
            fnameArray[countf] = ' ' 
            countf = countf+1
                      
            farray[:,:,:,countf] = f
            fnameArray[countf] = fname
            countf = countf+1

        if (orderMax >= 1):
                      
            farray[:,:,:,countf] = fx; 
            fnameArray[countf] = fname + 'x' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fy; 
            fnameArray[countf] = fname + 'y' 
            countf = countf+1

        if (orderMax >= 2):
                      
            farray[:,:,:,countf] = fxx; 
            fnameArray[countf] = fname + 'xx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyy; 
            fnameArray[countf] = fname + 'yy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxy; 
            fnameArray[countf] = fname + 'xy' 
            countf = countf+1

        if (orderMax >= 3):
                      
            farray[:,:,:,countf] = fxxx; 
            fnameArray[countf] = fname + 'xxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyy; 
            fnameArray[countf] = fname + 'yyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxy; 
            fnameArray[countf] = fname + 'xxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxyy; 
            fnameArray[countf] = fname + 'xyy' 
            countf = countf+1
            
        if (orderMax >= 4):
                      
            farray[:,:,:,countf] = fxxxx; 
            fnameArray[countf] = fname + 'xxxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyyy; 
            fnameArray[countf] = fname + 'yyyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxy; 
            fnameArray[countf] = fname + 'xxxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxyy; 
            fnameArray[countf] = fname + 'xxyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxyyy; 
            fnameArray[countf] = fname + 'xyyy' 
            countf = countf+1
            
        if (orderMax >= 5):
                      
            farray[:,:,:,countf] = fxxxxx; 
            fnameArray[countf] = fname + 'xxxxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyyyy; 
            fnameArray[countf] = fname + 'yyyyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxy; 
            fnameArray[countf] = fname + 'xxxxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxyy; 
            fnameArray[countf] = fname + 'xxxyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxyyy; 
            fnameArray[countf] = fname + 'xxyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxyyyy; 
            fnameArray[countf] = fname + 'xyyyy' 
            countf = countf+1
            
        if (orderMax >= 6):
                      
            farray[:,:,:,countf] = fxxxxxx; 
            fnameArray[countf] = fname + 'xxxxxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyyyyy; 
            fnameArray[countf] = fname + 'yyyyyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxxy; 
            fnameArray[countf] = fname + 'xxxxxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxyy; 
            fnameArray[countf] = fname + 'xxxxyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxyyy; 
            fnameArray[countf] = fname + 'xxxyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxyyyy; 
            fnameArray[countf] = fname + 'xxyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxyyyyy; 
            fnameArray[countf] = fname + 'xyyyyy' 
            countf = countf+1
            
        if (orderMax >= 7):
                      
            farray[:,:,:,countf] = fxxxxxxx; 
            fnameArray[countf] = fname + 'xxxxxxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyyyyyy; 
            fnameArray[countf] = fname + 'yyyyyyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxxxy; 
            fnameArray[countf] = fname + 'xxxxxxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxxyy; 
            fnameArray[countf] = fname + 'xxxxxyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxxyyy; 
            fnameArray[countf] = fname + 'xxxxyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxyyyy; 
            fnameArray[countf] = fname + 'xxxyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxyyyyy; 
            fnameArray[countf] = fname + 'xxyyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxyyyyyy; 
            fnameArray[countf] = fname + 'xyyyyyy' 
            countf = countf+1
            
        if (orderMax >= 8):
                      
            farray[:,:,:,countf] = fxxxxxxxx; 
            fnameArray[countf] = fname + 'xxxxxxxx' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fyyyyyyyy; 
            fnameArray[countf] = fname + 'yyyyyyyy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxxxxy; 
            fnameArray[countf] = fname + 'xxxxxxxy' 
            countf = countf+1
                      
            farray[:,:,:,countf] = fxxxxxxyy; 
            fnameArray[countf] = fname + 'xxxxxxyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxxxyyy; 
            fnameArray[countf] = fname + 'xxxxxyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxxyyyy; 
            fnameArray[countf] = fname + 'xxxxyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxxyyyyy; 
            fnameArray[countf] = fname + 'xxxyyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxxyyyyyy; 
            fnameArray[countf] = fname + 'xxyyyyyy' 
            countf = countf+1
            
            farray[:,:,:,countf] = fxyyyyyyy; 
            fnameArray[countf] = fname + 'xyyyyyyy' 
            countf = countf+1
            
    if not(returnf):
        farray = np.delete(farray,obj=1,axis=3)
        fnameArray = np.delete(fnameArray,obj=1)
        
# Higher degree derivatives
 
    if (degree >= 2):
        
        tempShape = farray.shape[-1]
        temp = np.zeros((M,N,T,1))
                
        for degree in range(2,degree+1):
            for count in range (1,tempShape):
                temp[:,:,:,0] = np.power(farray[:,:,:,count],degree)
                farray = np.concatenate((farray, temp), axis=3)
                fnameArray = np.append(fnameArray, fnameArray[count]*degree )
                
    return farray, fnameArray  
