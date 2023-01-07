"""
@author: Karan Jakhar

**WARNING**
This function de-normalizes the output weight (~w) from RVM
"""
"""
This fuction de-normalizes normalized weight ~w

~y = ~X.~w

~w = RVM(~y,~X)
~w is output from RVM

~y and ~w are normalized vectors
~X is normalized 2D matrix (column are normalized)

Normalizstion process for vectos or columns of the matrix

~y = (y - mean(y))/std(y)
~X = (X - mean(X, axis=0))/(std(X,axis=0))

The de-normalized equation would be of the following form

y = X.w + c

y: nX1 vector
X: nXd matrix
w: dX1 vector
c: dX1 vector

w = yStd/Xstd * ~w
c = -(yStd/Xstd * ~w)*XMean + yMean
  = -w * XMean + yMean

Input:
~w (weightMeanNorm): nx1 vector
mean(y) (yMean): scalar
std(y) (yStd): scalar
mean(X) (XMean): dX1 vector

returns:
weightMean: dX1 matrix
c: dX1 matrix
"""

"""
Example: d = 2

~y = ~X.~w
y = X[:,0]*w[0] + X[:,1]*w[1] + np.sum(c)

"""
import numpy as np

def denormalizeWeight(wNorm, yMean, yStd, XMean, XStd):
    
    wMean = yStd*np.multiply(wNorm,1/XStd)
    c = -np.multiply(wMean, XMean) + yMean
        
    return wMean, c   