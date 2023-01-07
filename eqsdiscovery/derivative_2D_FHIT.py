import numpy as np
from scipy.fftpack import fft2, ifft2
import math

def derivative_2D_FHIT(T,order,Tname):
    
## Calculate spatial derivatives for 2D_FHIT
# Derivatives are calculated in spectral space
# Boundary conditions are periodic in x and y spatial dimensions
# Length of domain 2*pi

# Tdash, 'Txxy'] = derivative_2D_FHIT(T, [2, 1], 'T')

## Input
# T: Input flow field: Square Matrix NxN
# order [orderX, orderY]: Array of order of derivatives in x and y spatial dimensions: [Interger (>=0), Integer (>=0)] 
# Tname: Name of the input flow field: Character 'T', 'Psi', 'tau', 'U',...

## Output
# Tdash: derivative of the flow field T: Square Matrix NxN
# TnameOut: Name of the output derivative of flow field: Characters

    orderX = order[0]
    orderY = order[1]

    # Validating user input
    if orderX < 0 or orderY < 0:
        raise ValueError("Order of derivatives must be 0 or positive")
    elif orderX == 0 and orderY == 0:
        raise ValueError("Both order of derivatives are 0, atleast one of them should be positive")

    # orderX 0,1,2,3,4...
    # orderY,0,1,2,3,4...

    Ngrid = T.shape[0];

    # Initialize grid
    
    pi = math.pi
    Lx = 2*pi

    kx = np.zeros(Ngrid)
    kx[0:int(Ngrid/2)+1]=2*pi*np.arange(0,Ngrid/2+1)/float(Lx)
    kx[int(Ngrid/2)+1:Ngrid]=-kx[int(Ngrid/2)-1:0:-1]

    # Making meshgrid 
    [Ky,Kx] = np.meshgrid(kx,kx)

    # Calculating derivatives in spectral space

    T_hat = fft2(T);
    Tdash_hat = np.multiply(np.multiply(((1j*Kx)**orderX), ((1j*Ky)**orderY)), T_hat);
    Tdash = np.real(ifft2(Tdash_hat)) ;

    # Naming the variable
    TnameOut = Tname;
    if orderX > 0:
        for count in range(1,orderX+1):
            TnameOut = TnameOut + 'x';

    if orderY > 0:
        for count in range(1,orderY+1):
            TnameOut = TnameOut + 'y';

    return Tdash, TnameOut