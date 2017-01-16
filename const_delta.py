from __future__ import division

import numpy 
from numpy import exp,sqrt,log
from diff_matrix import diff_matrix
import matplotlib.pyplot as plt

#import scipy.optimize

def constant_delta_T(T0,psiN0,Z,B):
    #T profile is different since rho_p depends on T
    return lambda psiN: T0*(1 -  numpy.fabs(Z)*B*(psiN - psiN0))**2

def constant_delta_dTdpsiN(T0,psiN0,Z,B):
    #T profile is different since rho_p depends on T
    return lambda psiN: 2*T0*(1 -  numpy.fabs(Z)*B*(psiN - psiN0))*(-numpy.fabs(Z)*B)

def constant_delta_X(eta0,psiN0,Z,deltaEta,BT,deltaT):
    #this formula assumes a T profile calculated by the above T!
    # BT should be the same in both
    return lambda psiN: eta0 * (1 -  numpy.fabs(Z)*BT*(psiN - psiN0))**(2*deltaEta/deltaT)

def constant_delta_dXdpsiN(eta0,psiN0,Z,deltaEta,BT,deltaT):
    #this formula assumes a T profile calculated by the above T!
    # BT should be the same in both
    return lambda psiN: -numpy.fabs(Z)*BT *(2*deltaEta/deltaT)*eta0*(1 -  numpy.fabs(Z)*BT*(psiN - psiN0))**(2.0*deltaEta/deltaT-1.0)

def single_ion_all_eta_Phi(etai,etae,Ti,Te,Zi):
    return lambda psiN: (Ti(psiN)*Te(psiN)/(Ti(psiN) + Zi*Te(psiN))) * log(Zi*etai(psiN)/etae(psiN))

def single_ion_all_eta_dPhidpsiN(etai,etae,Ti,Te,ddx_etai,ddx_etae,ddx_Ti,ddx_Te,Zi):
    return lambda psiN: (
        (1/(Ti(psiN) + Zi*Te(psiN))) * (
        (ddx_Ti(psiN) * Te(psiN) + Ti(psiN) * ddx_Te(psiN) - (ddx_Ti(psiN) + Zi*ddx_Te(psiN))*Ti(psiN)*Te(psiN)/(Ti(psiN) + Zi*Te(psiN))) * log(Zi*etai(psiN)/etae(psiN))
        + Ti(psiN)*Te(psiN) * (ddx_etai(psiN)/etai(psiN) - ddx_etae(psiN)/etae(psiN)))
    )


def single_ion_n_e_eta_i_Phi(etai,ne,Te,Ti,Zi):
    return lambda psiN: - (Ti(psiN)/Zi) * log(ne(psiN)/etai(psiN))


if __name__ == "__main__":

    ##### TUNEABLES ##############

    deltaT = 0.1
    deltaEta = 0.06
    deltaT_e = 0.01
    deltaEta_e = 0.02

    ###### GRID ###################
    
    psiN = numpy.linspace(0.8,1.099,100)
    D=diff_matrix(psiN[0],psiN[-1],len(psiN))

    ##### Constants ###############
    
    Z = 1
    Zi = Z
    Ze = -1
    RHat =1.3
    mHat = 1.0
    mHat_e = 0.000271933691768
    Delta = 0.0006  #based on He papar
    psiAHat = 0.0117  # based on He paper

    psiN0 = 0.9 # based on He paper
    
    # for T
    T0 = 1
    A = deltaT *psiAHat/(sqrt(mHat)*RHat*Delta)
    B = A/(2*sqrt(T0))
    
    T = constant_delta_T(T0,psiN0,Z,B)
    Ti=T
    dTdpsiN = constant_delta_dTdpsiN(T0,psiN0,Z,B)

    Ae = deltaT_e *psiAHat/(sqrt(mHat_e)*RHat*Delta)
    Be = Ae/(2*sqrt(T0))
    Te = constant_delta_T(T0,psiN0,Ze,Be)
    dTedpsiN = constant_delta_dTdpsiN(T0,psiN0,Ze,Be)

    ###### PLOT T #################

    plt.subplot(3, 1, 1)
    plt.plot(psiN,T(psiN))
    plt.plot(psiN,Te(psiN))
    plt.ylabel(r"$T$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(3, 1, 2)
    #num_dTdpsiN = numpy.dot(D,T(psiN))
    #num_dTedpsiN = numpy.dot(D,Te(psiN))
    #plt.plot(psiN,num_dTdpsiN)
    #plt.plot(psiN,num_dTedpsiN)
    plt.plot(psiN,dTdpsiN(psiN))
    plt.plot(psiN,dTedpsiN(psiN))
    plt.ylabel(r"$dT/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.subplot(3, 1, 3)
    deltaT = -Delta/(Z*psiAHat)*RHat*sqrt(mHat)*dTdpsiN(psiN)/sqrt(T(psiN))
    deltaTe = -Delta/(Z*psiAHat)*RHat*sqrt(mHat_e)*dTedpsiN(psiN)/sqrt(Te(psiN))
    plt.plot(psiN,deltaT)
    plt.plot(psiN,deltaTe)
    plt.ylabel(r"$\delta_T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,0.2])
    plt.savefig("exponentialT.pdf")
    plt.show()
    
    # for eta
    eta0 = 1.0
    A = deltaEta *psiAHat/(sqrt(mHat)*RHat*Delta)
    Beta = A/(2*sqrt(T0))
    eta = constant_delta_X(eta0,psiN0,Z,deltaEta,B,deltaT)
    detadpsiN = constant_delta_dXdpsiN(eta0,psiN0,Z,deltaEta,B,deltaT)

    eta0e= 1.0
    Ae = deltaEta_e *psiAHat/(sqrt(mHat_e)*RHat*Delta)
    Betae = Ae/(2*sqrt(T0))
    etae = constant_delta_X(eta0,psiN0,Ze,deltaEta_e,Be,deltaT_e)
    detaedpsiN = constant_delta_dXdpsiN(eta0,psiN0,Ze,deltaEta_e,Be,deltaT_e)

    ############# PLOT ETA ####################
    plt.subplot(3, 1, 1)
    plt.plot(psiN,eta(psiN))
    plt.plot(psiN,etae(psiN))
    plt.ylabel(r"$\eta$")
    plt.xlabel(r"$\psi_N$")


    plt.subplot(3, 1, 2)
    #num_detadpsiN = numpy.dot(D,eta(psiN))
    #num_detaedpsiN = numpy.dot(D,etae(psiN))
    #plt.plot(psiN,num_detadpsiN)
    #plt.plot(psiN,num_detaedpsiN)
    plt.plot(psiN,detadpsiN(psiN))
    plt.plot(psiN,detaedpsiN(psiN))
    plt.ylabel(r"$d\eta/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.subplot(3, 1, 3)
    deltaEta = -Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN))*detadpsiN(psiN)/eta(psiN)
    deltaEtae = -Delta/(numpy.fabs(Ze)*psiAHat)*RHat*sqrt(mHat_e*Te(psiN))*detaedpsiN(psiN)/etae(psiN)
    plt.plot(psiN,deltaEta)
    plt.plot(psiN,deltaEtae)
    plt.ylabel(r"$\delta_\eta$")
    plt.xlabel(r"$\psi_N$")
    plt.savefig("exponentialEta.pdf")
    plt.show()


    # for Phi
    
    Phi = single_ion_all_eta_Phi(eta,etae,Ti,Te,Zi)

    ############# PLOT Phi ####################
    plt.subplot(2, 1, 1)
    plt.plot(psiN,Phi(psiN))
    plt.ylabel(r"$\Phi$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(2, 1, 2)
    dPhidpsiN = numpy.dot(D,Phi(psiN))
    plt.plot(psiN,dPhidpsiN)
    plt.ylabel(r"$d\Phi/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.savefig("exponentialPhi.pdf")
    plt.show()

    # for n
    n = lambda psiN: eta(psiN)*exp(-Z*Phi(psiN)/T(psiN))
    n_e = lambda psiN: etae(psiN)*exp(-Ze*Phi(psiN)/Te(psiN))

    
    ############# PLOT n ####################
    plt.subplot(3, 1, 1)
    plt.plot(psiN,n(psiN))
    plt.plot(psiN,n_e(psiN))
    plt.ylabel(r"$n$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(3, 1, 2)
    dndpsiN = numpy.dot(D,n(psiN))
    dn_edpsiN = numpy.dot(D,n_e(psiN))
    plt.plot(psiN,dndpsiN)
    plt.plot(psiN,dn_edpsiN)
    plt.ylabel(r"$dn/d\psi_N$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(3, 1, 3)
    deltaN = -Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN))*dndpsiN/n(psiN)
    mHat_e = 0.000271933691768
    deltaN_e = -Delta/(Z*psiAHat)*RHat*sqrt(mHat_e*Te(psiN))*dn_edpsiN/n_e(psiN)
    plt.plot(psiN,deltaN)
    plt.plot(psiN,deltaN_e)
    #plt.ylim([0,0.2])
    plt.ylabel(r"$\delta_n$")
    plt.xlabel(r"$\psi_N$")

    
    plt.savefig("exponentialN.pdf")
    plt.show()


    ############# plot higher derivatives ##########
    """d2etadpsiN2 = numpy.dot(D,num_detadpsiN)
    d2TdpsiN2 = numpy.dot(D,num_dTdpsiN)
    
    delta2Eta = (Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN)))**2*d2etadpsiN2/eta(psiN)
    delta2T = (Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN)))**2*d2TdpsiN2/T(psiN)

    delta1Eta = (Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN)))**1*num_detadpsiN/eta(psiN)
    delta1T = (Delta/(Z*psiAHat)*RHat*sqrt(mHat*T(psiN)))**1*num_dTdpsiN/T(psiN)
    

    print delta2Eta
    print delta1Eta*(delta1Eta - delta1T/2)
    
    print delta2T
    print delta1T*(delta1T - delta1T/2)
    
    plt.plot(psiN,delta2Eta)
    plt.plot(psiN,delta2T)
    plt.ylabel(r"$\delta^2$")
    plt.xlabel(r"$\psi_N$")

    plt.savefig("delta2.pdf")
    plt.show()"""
    
