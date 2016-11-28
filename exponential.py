from __future__ import division

import numpy 
from numpy import exp,sqrt
from diff_matrix import diff_matrix
import matplotlib.pyplot as plt

#import scipy.optimize

def constant_delta_T(T0,psiN0,Z,B):
    #T profile is different since rho_p depends on T
    T =lambda psiN: T0*(1 -  numpy.fabs(Z)*B*(psiN - psiN0))**2
    return T

def constant_delta(eta0,psiN0,Z,Beta,BT):
    #this formula assumes a T profile calculated by the above T!
    # B should be the same in both
    eta =lambda psiN: eta0 * (1 -  numpy.fabs(Z)*B*(psiN - psiN0))**(2*numpy.fabs(Z)*Beta)
    return eta
    
if __name__ == "__main__":
    
    Z = 1
    RHat =1
    mHat = 1
    Delta = 0.0006  #based on He papar
    psiAHat = 0.0117  # based on He paper
    # for T
    T0 = 1
    psiN0 = 0.9
    delta = 0.1
    A = delta *psiAHat/(sqrt(mHat)*RHat*Delta)
    B = A/(2*sqrt(T0))
    
    T = constant_delta_T(T0,psiN0,Z,B)

    # for eta
    eta0 = 1
    psiN0 = 0.9
    delta = 0.1
    A = delta *psiAHat/(sqrt(mHat)*RHat*Delta)
    Beta = A/(2*sqrt(T0))
    C = A/(2*sqrt(T0)*sqrt(eta0))

    eta = constant_delta(eta0,psiN0,Z,Beta,B)

    ###### GRID ###################
    
    psiN = numpy.linspace(-2,1,100)
    D=diff_matrix(psiN[0],psiN[-1],len(psiN))

    ###### PLOT T #################

    plt.subplot(3, 1, 1)
    plt.plot(psiN,T(psiN))
    plt.ylabel(r"$T$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(3, 1, 2)
    dTdpsiN = numpy.dot(D,T(psiN))
    plt.plot(psiN,dTdpsiN)
    plt.ylabel(r"$dT/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.subplot(3, 1, 3)
    deltaT = -Delta/(Z*psiAHat)*RHat*sqrt(mHat)*dTdpsiN/sqrt(T(psiN))
    plt.plot(psiN,deltaT)
    plt.ylabel(r"$\delta_T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,0.2])
    plt.savefig("exponentialT.pdf")
    plt.show()

    ############# plot ETA ####################
    plt.subplot(3, 1, 1)
    plt.plot(psiN,eta(psiN))
    plt.ylabel(r"$\eta$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(3, 1, 2)
    detadpsiN = numpy.dot(D,eta(psiN))
    plt.plot(psiN,detadpsiN)
    plt.ylabel(r"$d\eta/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.subplot(3, 1, 3)
    deltaEta = -Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN))*detadpsiN/eta(psiN)
    plt.plot(psiN,deltaEta)
    plt.ylabel(r"$\delta_\eta$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,0.2])
    plt.savefig("exponentialEta.pdf")
    plt.show()
