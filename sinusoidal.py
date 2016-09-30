from __future__ import division

import numpy as np
from numpy import sin,cos, pi
import matplotlib.pyplot as plt

import scipy.optimize

def sin_profile(min_val,max_val,psiMin,psiMax):
    A = max_val - min_val
    delta_psiN = psiMax - psiMin
    k= 2*pi/delta_psiN
    return lambda x: max_val - A*(1-cos(k*(x-psiMin)))/2.0

if __name__ == "__main__":
    xped=1.0
    w = 2
    dxdr = 0.1
    xsol = xped - dxdr*w
    
    psiN0 = 1.0
    psiN1 = psiN0 + w*2
    psiN=np.linspace(psiN0,psiN1)
    x=sin_profile(xsol,xped,psiN0,psiN1)
    plt.plot(psiN,x(psiN))
    plt.show()
