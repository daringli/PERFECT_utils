from const_delta import constant_delta_T, constant_delta_dTdpsiN, constant_delta_X, constant_delta_dXdpsiN, single_ion_all_eta_Phi
from mtanh import mtanh_transition, ddx_mtanh_transition
from diff_matrix import diff_matrix

import matplotlib.pyplot as plt
import numpy
from numpy import exp,sqrt,log



if __name__ == "__main__":

    only_mtanh=False
    ##### TUNEABLES ##############

    deltaTi = 0.1
    deltaEtai = 0.1
    deltaTe = 0.005
    deltaEtae = 0.05

    eta0 = 1.0
    T0 = 1.0

    psiN0 = 0.9 # based on He paper

    p = psiN0 + 0.02
    w = 0.02
    
    if only_mtanh:
        #Make true to get ordinary mtanh profiles
        deltaTi = 0.00000000000001
        deltaEtai = 0.00000000000001
        deltaTe = 0.00000000000001
        deltaEtae = 0.00000000000001
        p = psiN0 + 0.05
        w = 0.05
        eta0 = 1.0
        T0 = 1.0
    


    ##### GRID #######
    
    psiN = numpy.linspace(0.8,1.099,100)
    D=diff_matrix(psiN[0],psiN[-1],len(psiN))

    #### constants
    
    Zi = 1
    Ze = -1
    RHat =1.3
    mHati = 1.0
    mHate = 0.000271933691768
    Delta = 0.0006  #based on He papar
    psiAHat = 0.0117  # based on He paper


    #lazy core profiles
    T0c=T0
    eta0c=eta0
    a = 1.0
    core_Ti = lambda psiN: T0c - a*(psiN-psiN0)
    core_Te = lambda psiN: T0c - a*(psiN-psiN0)
    core_etai = lambda psiN: eta0c - a*(psiN-psiN0)
    core_etae = lambda psiN: eta0c - a*(psiN-psiN0)

    core_dTidpsiN = lambda psiN: -a + 0.0*psiN
    core_dTedpsiN = lambda psiN: -a + 0.0*psiN
    core_detaidpsiN = lambda psiN: -a + 0.0*psiN
    core_detaedpsiN = lambda psiN: -a + 0.0*psiN

    #constant delta Ti

    if only_mtanh:
        T0=T0c/5
        eta0=eta0c/5
    
    ATi = deltaTi *psiAHat/(sqrt(mHati)*RHat*Delta)
    BTi = ATi/(2*sqrt(T0))
    Ti = constant_delta_T(T0,psiN0,Zi,BTi)
    dTidpsiN = constant_delta_dTdpsiN(T0,psiN0,Zi,BTi)
    #constant delta Te
    ATe = deltaTe *psiAHat/(sqrt(mHate)*RHat*Delta)
    BTe = ATe/(2*sqrt(T0))
    Te = constant_delta_T(T0,psiN0,Ze,BTe)
    dTedpsiN = constant_delta_dTdpsiN(T0,psiN0,Ze,BTe)

    #constant delta etai
    etai = constant_delta_X(eta0,psiN0,Zi,deltaEtai,BTi,deltaTi)
    detaidpsiN = constant_delta_dXdpsiN(eta0,psiN0,Zi,deltaEtai,BTi,deltaTi)
    #constant delta etae
    etae = constant_delta_X(eta0,psiN0,Ze,deltaEtae,BTe,deltaTe)
    detaedpsiN = constant_delta_dXdpsiN(eta0,psiN0,Ze,deltaEtae,BTe,deltaTe)    

    # full profiles
    full_Ti = mtanh_transition(Ti,core_Ti,w,p)
    full_dTidpsiN = ddx_mtanh_transition(Ti,core_Ti,dTidpsiN,core_dTidpsiN,w,p)
    full_Te = mtanh_transition(Te,core_Te,w,p)
    full_dTedpsiN = ddx_mtanh_transition(Te,core_Te,dTedpsiN,core_dTedpsiN,w,p)
    full_etai = mtanh_transition(etai,core_etai,w,p)
    full_detaidpsiN = ddx_mtanh_transition(etai,core_etai,detaidpsiN,core_detaidpsiN,w,p)
    full_etae = mtanh_transition(etae,core_etae,w,p)
    full_detaedpsiN = ddx_mtanh_transition(etae,core_etae,detaedpsiN,core_detaedpsiN,w,p)

    # potential
    full_Phi = single_ion_all_eta_Phi(full_etai,full_etae,full_Ti,full_Te,Zi)

    # densities
    full_ni = lambda psiN: full_etai(psiN)*exp(-Zi*full_Phi(psiN)/full_Ti(psiN))
    full_ne = lambda psiN: full_etae(psiN)*exp(-Ze*full_Phi(psiN)/full_Te(psiN))


    ###### PLOT T #################

    plt.subplot(3, 1, 1)
    plt.plot(psiN,full_Ti(psiN))
    plt.plot(psiN,full_Te(psiN))
    plt.plot(psiN,Ti(psiN),ls=":",color="darkblue")
    plt.plot(psiN,Te(psiN),ls=":",color="darkgreen")
    plt.ylabel(r"$T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1.2])

    plt.subplot(3, 1, 2)
    num_full_dTidpsiN = numpy.dot(D,full_Ti(psiN))
    num_full_dTedpsiN = numpy.dot(D,full_Te(psiN))
    #plt.plot(psiN,num_dTdpsiN)
    #plt.plot(psiN,num_dTedpsiN)
    plt.plot(psiN,full_dTidpsiN(psiN))
    plt.plot(psiN,full_dTedpsiN(psiN))
    plt.plot(psiN,dTidpsiN(psiN),ls=":",color="darkblue")
    plt.plot(psiN,dTedpsiN(psiN),ls=":",color="darkgreen")
    plt.ylabel(r"$dT/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-12,0])
    
    plt.subplot(3, 1, 3)
    deltaT = -Delta/(Zi*psiAHat)*RHat*sqrt(mHati)*full_dTidpsiN(psiN)/sqrt(full_Ti(psiN))
    deltaTe = -Delta/(Zi*psiAHat)*RHat*sqrt(mHate)*full_dTedpsiN(psiN)/sqrt(full_Te(psiN))
    plt.plot(psiN,deltaT)
    plt.plot(psiN,deltaTe)
    plt.ylabel(r"$\delta_T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,0.2])
    plt.savefig("mtanh_constDeltaT.pdf")
    plt.show()
    

    ############# PLOT ETA ####################
    plt.subplot(3, 1, 1)
    plt.plot(psiN,full_etai(psiN))
    plt.plot(psiN,full_etae(psiN))
    plt.plot(psiN,etai(psiN),ls=":",color="darkblue")
    plt.plot(psiN,etae(psiN),ls=":",color="darkgreen")
    plt.ylabel(r"$\eta$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1.2])


    plt.subplot(3, 1, 2)
    plt.plot(psiN,full_detaidpsiN(psiN))
    plt.plot(psiN,full_detaedpsiN(psiN))
    num_full_detaidpsiN= numpy.dot(D,full_etai(psiN))
    num_full_detaedpsiN= numpy.dot(D,full_etae(psiN))
    plt.plot(psiN,full_detaidpsiN(psiN),color="blue")
    plt.plot(psiN,full_detaedpsiN(psiN),color="green")
    #plt.plot(psiN,num_full_detaidpsiN,ls=":",color="blue")
    #plt.plot(psiN,num_full_detaedpsiN,ls=":",color="green")
    plt.plot(psiN,detaidpsiN(psiN),ls=":",color="darkblue")
    plt.plot(psiN,detaedpsiN(psiN),ls=":",color="darkgreen")
    plt.ylabel(r"$d\eta/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-30,0])
    
    plt.subplot(3, 1, 3)
    deltaEta = -Delta/(Zi*psiAHat)*RHat*sqrt(mHati*Ti(psiN))*full_detaidpsiN(psiN)/full_etai(psiN)
    deltaEtae = -Delta/(numpy.fabs(Ze)*psiAHat)*RHat*sqrt(mHate*Te(psiN))*full_detaedpsiN(psiN)/full_etae(psiN)
    
        
    plt.plot(psiN,deltaEta)
    plt.plot(psiN,deltaEtae)
    plt.ylabel(r"$\delta_\eta$")
    plt.xlabel(r"$\psi_N$")
    plt.savefig("mtanh_constDeltaEta.pdf")
    plt.show()

    ############# PLOT Phi ####################
    plt.subplot(2, 1, 1)
    plt.plot(psiN,full_Phi(psiN))
    plt.ylabel(r"$\Phi$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(2, 1, 2)
    num_full_dPhidpsiN = numpy.dot(D,full_Phi(psiN))
    plt.plot(psiN,num_full_dPhidpsiN)
    plt.ylabel(r"$d\Phi/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    
    plt.savefig("mtanh_constDeltaPhi.pdf.pdf")
    plt.show()

    ############# PLOT n ####################
    plt.subplot(3, 1, 1)
    plt.plot(psiN,full_ni(psiN))
    plt.plot(psiN,full_ne(psiN))
    plt.ylabel(r"$n$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1.2])

    plt.subplot(3, 1, 2)
    num_full_dnidpsiN = numpy.dot(D,full_ni(psiN))
    num_full_dnedpsiN = numpy.dot(D,full_ne(psiN))
    plt.plot(psiN,num_full_dnidpsiN)
    plt.plot(psiN,num_full_dnedpsiN)
    plt.ylabel(r"$dn/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-30,0])

    plt.subplot(3, 1, 3)
    deltaNi = -Delta/(numpy.fabs(Zi)*psiAHat)*RHat*sqrt(mHati*full_Ti(psiN))*num_full_dnidpsiN/full_ni(psiN)
    mHat_e = 0.000271933691768
    deltaNe = -Delta/(numpy.fabs(Ze)*psiAHat)*RHat*sqrt(mHate*full_Te(psiN))*num_full_dnedpsiN/full_ne(psiN)
    plt.plot(psiN,deltaNi)
    plt.plot(psiN,deltaNe)
    #plt.ylim([0,0.2])
    plt.ylabel(r"$\delta_n$")
    plt.xlabel(r"$\psi_N$")

    
    plt.savefig("mtanh_constDeltaN.pdf")
    plt.show()


    ###### PLOT higher derivatives ##########
    num_full_d2etaidpsiN2 = numpy.dot(D,num_full_detaidpsiN)
    num_full_d2etaedpsiN2 = numpy.dot(D,num_full_detaedpsiN)
    
    num_full_d2TidpsiN2 = numpy.dot(D,num_full_dTidpsiN)
    num_full_d2TedpsiN2 = numpy.dot(D,num_full_dTedpsiN)

    delta2Etai = (Delta/(Zi*psiAHat)*RHat*sqrt(mHati*full_Ti(psiN)))**2*num_full_d2etaidpsiN2/full_etai(psiN)
    delta2Ti = (Delta/(Zi*psiAHat)*RHat*sqrt(mHati*full_Ti(psiN)))**2*num_full_d2TidpsiN2/full_Ti(psiN)

    delta2Etae = (Delta/(Ze*psiAHat)*RHat*sqrt(mHate*full_Te(psiN)))**2*num_full_d2etaedpsiN2/full_etae(psiN)
    delta2Te = (Delta/(Ze*psiAHat)*RHat*sqrt(mHate*full_Te(psiN)))**2*num_full_d2TedpsiN2/full_Te(psiN)

    plt.plot(psiN,delta2Etai)
    plt.plot(psiN,delta2Ti)
    plt.plot(psiN,delta2Etae)
    plt.plot(psiN,delta2Te)
    
    plt.ylabel(r"$\delta^2$")
    plt.xlabel(r"$\psi_N$")

    plt.savefig("delta2.pdf")
    plt.show()
                    
