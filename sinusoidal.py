from __future__ import division

import numpy as np
from numpy import sin,cos, pi
import matplotlib.pyplot as plt
from diff_matrix import diff_matrix

import scipy.optimize

def sin_profile(min_val,max_val,psiMin,psiMax):
    A = max_val - min_val
    delta_psiN = psiMax - psiMin
    k= 2*pi/delta_psiN
    return lambda x: max_val - A*(1-cos(k*(x-psiMin)))/2.0

def ddx_sin_profile(min_val,max_val,psiMin,psiMax):
    A = max_val - min_val
    delta_psiN = psiMax - psiMin
    k= 2*pi/delta_psiN
    return lambda x: -0.5*k*A*sin(k*(x-psiMin))

def parameter_wrapper(xPed,xPedGrad,psiMinPed,psiMaxPed):
    max_val = xPed
    w = psiMaxPed - psiMinPed
    min_val = xPed + xPedGrad*w
    psiMin = psiMinPed
    psiMax = psiMin + 2*w
    return (min_val,max_val,psiMin,psiMax)

def generate_sin_profile(xPed,xPedGrad,psiMinPed,psiMaxPed):
    (min_val,max_val,psiMin,psiMax) = parameter_wrapper(xPed,xPedGrad,psiMinPed,psiMaxPed)
    return (sin_profile(min_val,max_val,psiMin,psiMax),ddx_sin_profile(min_val,max_val,psiMin,psiMax))

if __name__ == "__main__":
    xped=1.0
    w = 2
    dxdr = 0.1
    xsol = xped - dxdr*w
    
    psiN0 = 1.0
    psiN1 = psiN0 + w*2
    psiN=np.linspace(psiN0,psiN1)
    D=diff_matrix(psiN[0],psiN[-1],len(psiN))

    x=sin_profile(xsol,xped,psiN0,psiN1)
    dxdpsiN = ddx_sin_profile(xsol,xped,psiN0,psiN1)
    numerical_dxdpsiN = np.dot(D,x(psiN))

    plt.subplot(2, 1, 1)
    plt.hold(True)
    plt.plot(psiN,x(psiN))

    plt.subplot(2, 1, 2)
    plt.hold(True)
    plt.plot(psiN,dxdpsiN(psiN))
    plt.plot(psiN,numerical_dxdpsiN)

    
    plt.show()
