"""
@author: Karan Jakhar

Caompute derivatives for 2D_FHIT
"""

import numpy as np
from scipy.fftpack import fft, ifft, fft2, ifft2
import math
import sys

## Taking derivative     
def derivative(f, kx, ky, order, fname):
    
    if order > 8:
        print('**Error** Derivatives upto order 8 are supported, input order greater than 8')
        sys.exit()
    
    M = len(kx)
    N = len(ky)
    
    [kky,kkx] = np.meshgrid(kx,ky) 
        
    if (order>=2):
        kk2xx = np.multiply(kkx,kkx)
        kk2yy = np.multiply(kky,kky)
        kk2xy = np.multiply(kkx,kky)
        
    if (order>=3):
        kk3xxx = np.multiply(kk2xx,kkx)
        kk3yyy = np.multiply(kk2yy,kky)
        kk3xxy = np.multiply(kk2xx,kky)
        kk3xyy = np.multiply(kk2xy,kky)
        
    if (order>=4):
        kk4xxxx = np.multiply(kk3xxx,kkx)
        kk4yyyy = np.multiply(kk3yyy,kky)
        kk4xxxy = np.multiply(kk3xxx,kky)
        kk4xxyy = np.multiply(kk3xxy,kky)
        kk4xyyy = np.multiply(kk3xyy,kky)
        
    if (order>=5):
        kk5xxxxx = np.multiply(kk4xxxx,kkx)
        kk5yyyyy = np.multiply(kk4yyyy,kky)
        kk5xxxxy = np.multiply(kk4xxxx,kky)
        kk5xxxyy = np.multiply(kk4xxxy,kky)
        kk5xxyyy = np.multiply(kk4xxyy,kky)
        kk5xyyyy = np.multiply(kk4xyyy,kky)
        
    if (order>=6):
        kk6xxxxxx = np.multiply(kk5xxxxx,kkx)
        kk6yyyyyy = np.multiply(kk5yyyyy,kky)
        kk6xxxxxy = np.multiply(kk5xxxxx,kky)
        kk6xxxxyy = np.multiply(kk5xxxxy,kky)
        kk6xxxyyy = np.multiply(kk5xxxyy,kky)
        kk6xxyyyy = np.multiply(kk5xxyyy,kky)
        kk6xyyyyy = np.multiply(kk5xyyyy,kky)
        
    if (order>=7):
        kk7xxxxxxx = np.multiply(kk6xxxxxx,kkx)
        kk7yyyyyyy = np.multiply(kk6yyyyyy,kky)
        kk7xxxxxxy = np.multiply(kk6xxxxxx,kky)
        kk7xxxxxyy = np.multiply(kk6xxxxxy,kky)
        kk7xxxxyyy = np.multiply(kk6xxxxyy,kky)
        kk7xxxyyyy = np.multiply(kk6xxxyyy,kky)
        kk7xxyyyyy = np.multiply(kk6xxyyyy,kky)
        kk7xyyyyyy = np.multiply(kk6xyyyyy,kky)
        
    if (order>=8):
        kk8xxxxxxxx = np.multiply(kk7xxxxxxx,kkx)
        kk8yyyyyyyy = np.multiply(kk7yyyyyyy,kky)
        kk8xxxxxxxy = np.multiply(kk7xxxxxxx,kky)
        kk8xxxxxxyy = np.multiply(kk7xxxxxxy,kky)
        kk8xxxxxyyy = np.multiply(kk7xxxxxyy,kky)
        kk8xxxxyyyy = np.multiply(kk7xxxxyyy,kky)
        kk8xxxyyyyy = np.multiply(kk7xxxyyyy,kky)
        kk8xxyyyyyy = np.multiply(kk7xxyyyyy,kky)
        kk8xyyyyyyy = np.multiply(kk7xyyyyyy,kky)
        

    fhat = fft2(f)

    fhat[0,0]=0

    if (order==1):
            
        fxhat = 1j*np.multiply(kkx,fhat)
        fyhat = 1j*np.multiply(kky,fhat)
                    
        fx = np.real(ifft2(fxhat))
        fy = np.real(ifft2(fyhat))

    if (order==2):
            
        fxxhat = -1*np.multiply(kk2xx,fhat)
        fyyhat = -1*np.multiply(kk2yy,fhat)
        fxyhat = -1*np.multiply(kk2xy,fhat)

        fxx = np.real(ifft2(fxxhat))
        fyy = np.real(ifft2(fyyhat))
        fxy = np.real(ifft2(fxyhat))
    
    if (order==3):
            
        fxxxhat = -1j*np.multiply(kk3xxx,fhat)
        fyyyhat = -1j*np.multiply(kk3yyy,fhat)
        fxxyhat = -1j*np.multiply(kk3xxy,fhat)
        fxyyhat = -1j*np.multiply(kk3xyy,fhat)
    
        fxxx = np.real(ifft2(fxxxhat))
        fyyy = np.real(ifft2(fyyyhat))
        fxxy = np.real(ifft2(fxxyhat))
        fxyy = np.real(ifft2(fxyyhat))
        
    if (order==4):
            
        fxxxxhat = 1*np.multiply(kk4xxxx,fhat)
        fyyyyhat = 1*np.multiply(kk4yyyy,fhat)
        fxxxyhat = 1*np.multiply(kk4xxxy,fhat)
        fxxyyhat = 1*np.multiply(kk4xxyy,fhat)
        fxyyyhat = 1*np.multiply(kk4xyyy,fhat)
    
        fxxxx = np.real(ifft2(fxxxxhat))
        fyyyy = np.real(ifft2(fyyyyhat))
        fxxxy = np.real(ifft2(fxxxyhat))
        fxxyy = np.real(ifft2(fxxyyhat))
        fxyyy = np.real(ifft2(fxyyyhat))
        
    if (order==5):
            
        fxxxxxhat = 1j*np.multiply(kk5xxxxx,fhat)
        fyyyyyhat = 1j*np.multiply(kk5yyyyy,fhat)
        fxxxxyhat = 1j*np.multiply(kk5xxxxy,fhat)
        fxxxyyhat = 1j*np.multiply(kk5xxxyy,fhat)
        fxxyyyhat = 1j*np.multiply(kk5xxyyy,fhat)
        fxyyyyhat = 1j*np.multiply(kk5xyyyy,fhat)
    
        fxxxxx = np.real(ifft2(fxxxxxhat))
        fyyyyy = np.real(ifft2(fyyyyyhat))
        fxxxxy = np.real(ifft2(fxxxxyhat))
        fxxxyy = np.real(ifft2(fxxxyyhat))
        fxxyyy = np.real(ifft2(fxxyyyhat))
        fxyyyy = np.real(ifft2(fxyyyyhat))
        
    if (order==6):
            
        fxxxxxxhat = -1*np.multiply(kk6xxxxxx,fhat)
        fyyyyyyhat = -1*np.multiply(kk6yyyyyy,fhat)
        fxxxxxyhat = -1*np.multiply(kk6xxxxxy,fhat)
        fxxxxyyhat = -1*np.multiply(kk6xxxxyy,fhat)
        fxxxyyyhat = -1*np.multiply(kk6xxxyyy,fhat)
        fxxyyyyhat = -1*np.multiply(kk6xxyyyy,fhat)
        fxyyyyyhat = -1*np.multiply(kk6xyyyyy,fhat)
    
        fxxxxxx = np.real(ifft2(fxxxxxxhat))
        fyyyyyy = np.real(ifft2(fyyyyyyhat))
        fxxxxxy = np.real(ifft2(fxxxxxyhat))
        fxxxxyy = np.real(ifft2(fxxxxyyhat))
        fxxxyyy = np.real(ifft2(fxxxyyyhat))
        fxxyyyy = np.real(ifft2(fxxyyyyhat))
        fxyyyyy = np.real(ifft2(fxyyyyyhat))
        
    if (order==7):
            
        fxxxxxxxhat = -1j*np.multiply(kk7xxxxxxx,fhat)
        fyyyyyyyhat = -1j*np.multiply(kk7yyyyyyy,fhat)
        fxxxxxxyhat = -1j*np.multiply(kk7xxxxxxy,fhat)
        fxxxxxyyhat = -1j*np.multiply(kk7xxxxxyy,fhat)
        fxxxxyyyhat = -1j*np.multiply(kk7xxxxyyy,fhat)
        fxxxyyyyhat = -1j*np.multiply(kk7xxxyyyy,fhat)
        fxxyyyyyhat = -1j*np.multiply(kk7xxyyyyy,fhat)
        fxyyyyyyhat = -1j*np.multiply(kk7xyyyyyy,fhat)
    
        fxxxxxxx = np.real(ifft2(fxxxxxxxhat))
        fyyyyyyy = np.real(ifft2(fyyyyyyyhat))
        fxxxxxxy = np.real(ifft2(fxxxxxxyhat))
        fxxxxxyy = np.real(ifft2(fxxxxxyyhat))
        fxxxxyyy = np.real(ifft2(fxxxxyyyhat))
        fxxxyyyy = np.real(ifft2(fxxxyyyyhat))
        fxxyyyyy = np.real(ifft2(fxxyyyyyhat))
        fxyyyyyy = np.real(ifft2(fxyyyyyyhat))
        
    if (order==8):
            
        fxxxxxxxxhat = 1*np.multiply(kk8xxxxxxxx,fhat)
        fyyyyyyyyhat = 1*np.multiply(kk8yyyyyyyy,fhat)
        fxxxxxxxyhat = 1*np.multiply(kk8xxxxxxxy,fhat)
        fxxxxxxyyhat = 1*np.multiply(kk8xxxxxxyy,fhat)
        fxxxxxyyyhat = 1*np.multiply(kk8xxxxxyyy,fhat)
        fxxxxyyyyhat = 1*np.multiply(kk8xxxxyyyy,fhat)
        fxxxyyyyyhat = 1*np.multiply(kk8xxxyyyyy,fhat)
        fxxyyyyyyhat = 1*np.multiply(kk8xxyyyyyy,fhat)
        fxyyyyyyyhat = 1*np.multiply(kk8xyyyyyyy,fhat)
    
        fxxxxxxxx = np.real(ifft2(fxxxxxxxxhat))
        fyyyyyyyy = np.real(ifft2(fyyyyyyyyhat))
        fxxxxxxxy = np.real(ifft2(fxxxxxxxyhat))
        fxxxxxxyy = np.real(ifft2(fxxxxxxyyhat))
        fxxxxxyyy = np.real(ifft2(fxxxxxyyyhat))
        fxxxxyyyy = np.real(ifft2(fxxxxyyyyhat))
        fxxxyyyyy = np.real(ifft2(fxxxyyyyyhat))
        fxxyyyyyy = np.real(ifft2(fxxyyyyyyhat))
        fxyyyyyyy = np.real(ifft2(fxyyyyyyyhat))


    if (order==1):
        fnameOut = [fname+'x', fname+'y']
        return fx, fy, fnameOut
    
    if (order==2):
        fnameOut = [fname+'xx', fname+'yy', fname+'xy']
        return fxx, fyy, fxy, fnameOut
    
    if (order==3):
        fnameOut = [fname+'xxx', fname+'yyy', fname+'xxy', fname+'xyy']
        return fxxx, fyyy, fxxy, fxyy, fnameOut
    
    if (order==4):
        fnameOut = [fname+'xxxx', fname+'yyyy', fname+'xxxy', fname+'xxyy', fname+'xyyy']
        return fxxxx, fyyyy, fxxxy, fxxyy, fxyyy, fnameOut
    
    if (order==5):
        fnameOut = [fname+'xxxxx', fname+'yyyyy', fname+'xxxxy', fname+'xxxyy', fname+'xxyyy', fname+'xyyyy']
        return fxxxxx, fyyyyy, fxxxxy, fxxxyy, fxxyyy, fxyyyy, fnameOut
    
    if (order==6):
        fnameOut = [fname+'xxxxxx', fname+'yyyyyy', fname+'xxxxxy', fname+'xxxxyy', fname+'xxxyyy', fname+'xxyyyy', fname+'xyyyyy']
        return fxxxxxx, fyyyyyy, fxxxxxy, fxxxxyy, fxxxyyy, fxxyyyy, fxyyyyy, fnameOut
    
    if (order==7):
        fnameOut = [fname+'xxxxxxx', fname+'yyyyyyy', fname+'xxxxxxy', fname+'xxxxxyy', fname+'xxxxyyy', fname+'xxxyyyy', fname+'xxyyyyy', fname+'xyyyyyy']
        return fxxxxxxx, fyyyyyyy, fxxxxxxy, fxxxxxyy, fxxxxyyy, fxxxyyyy, fxxyyyyy, fxyyyyyy, fnameOut
    
    if (order==8):
        fnameOut = [fname+'xxxxxxxx', fname+'yyyyyyyy', fname+'xxxxxxxy', fname+'xxxxxxyy', fname+'xxxxxyyy', fname+'xxxxyyyy', fname+'xxxyyyyy', fname+'xxyyyyyy', fname+'xyyyyyyy']
        return fxxxxxxxx, fyyyyyyyy, fxxxxxxxy, fxxxxxxyy, fxxxxxyyy, fxxxxyyyy, fxxxyyyyy, fxxyyyyyy, fxyyyyyyy, fnameOut