"""
This fuction normalizes a vector or columns of a 2D matrix

Input:
y: Vector or a 2D matrix

returns:
yNorm: Normalized vector or a 2D matrix y
yMean: Mean of a vector or columns of a 2D matrix 
yStd: Standard Deviation of a vector or columns of a 2D matrix

"""
import sys
import numpy as np

def normalize(y):
    
    if (len(y.shape)==1):
        # Normalizing a vector
        
        yMean = np.mean(y)
        yStd = np.std(y)
        
        yNorm = (y-yMean)/yStd
                
    elif (len(y.shape)==2):
        # Normalizing a matrix along columns
        
        yMean = np.mean(y,axis=0)
        yStd = np.std(y,axis=0)
        
        yNorm = np.zeros(y.shape)
        for count in range(0,y.shape[1]):
            yNorm[:,count] = (y[:,count]-yMean[count])/yStd[count]
                                    
    elif (len(y.shape)>2):
        sys.exit("Error: Please input a vector or 2D matrix for normalization")
        
    return yNorm, yMean, yStd    