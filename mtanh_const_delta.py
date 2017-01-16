from __future__ import division

from const_delta import constant_delta_T, constant_delta_dTdpsiN, constant_delta_X, constant_delta_dXdpsiN, single_ion_all_eta_Phi, single_ion_all_eta_dPhidpsiN
from mtanh import mtanh_transition, ddx_mtanh_transition
from diff_matrix import diff_matrix

import matplotlib.pyplot as plt
import numpy
from numpy import exp,sqrt,log
import scipy.optimize


def eta_parameter_wrapper(Zi,Ze,mi,me,T0i,T0e,deltaTi,deltaTe,ATi,ATe,Ti,Te,n_ped,dndr_core,n_LCFS,ped_width,ped_pos):
    #eta and n overlap in the core:
    eta0 = n_ped
    eta0_core = eta0
    detadx_core = dndr_core

    #n_LFCS = n_ped + dndr_ped * ped_width # what we want
    print n_LCFS
    x_LCFS = ped_pos + ped_width # where we want it
    print x_LCFS
    deltaEtai = 0.1
    def f(deltaEtae): # that which gives it
        etae = mtanh_const_delta_X(Ze,eta0,deltaEtae,T0e,ATe,deltaTe,eta0_core,detadx_core,ped_width,ped_pos)
        etai = mtanh_const_delta_X(Zi,eta0,deltaEtai,T0i,ATi,deltaTi,eta0_core,detadx_core,ped_width,ped_pos)
        Phi = single_ion_all_eta_Phi(etai,etae,Ti,Te,Zi)
        return etae(x_LCFS)*exp(-Ze*Phi(x_LCFS)/Te(x_LCFS)) - n_LCFS
    #first_guess = -detadx_core/eta0_core # how we view it
    first_guess = 0.02 # how we view it
    deltaEtae = scipy.optimize.fsolve(f,first_guess)[0] # how it is
    print "delta_eta_e:" + str(deltaEtae)
    #deltaEtai = sqrt(mi)/sqrt(me) * deltaEtae
    etae = mtanh_const_delta_X(Ze,eta0,deltaEtae,T0e,ATe,deltaTe,eta0_core,detadx_core,ped_width,ped_pos)
    etai = mtanh_const_delta_X(Zi,eta0,deltaEtai,T0i,ATi,deltaTi,eta0_core,detadx_core,ped_width,ped_pos)
    ddx_etae = mtanh_const_delta_X_ddpsiN(Ze,eta0,deltaEtae,T0e,ATe,deltaTe,eta0_core,detadx_core,ped_width,ped_pos)
    ddx_etai = mtanh_const_delta_X_ddpsiN(Zi,eta0,deltaEtai,T0i,ATi,deltaTi,eta0_core,detadx_core,ped_width,ped_pos)
    return (etae,etai,ddx_etae,ddx_etai)

def T_parameter_wrapper(Z,T0,T_LCFS,ddx_T_core,ped_width,ped_pos):
    T0c = T0
    x_LCFS = ped_pos + ped_width # where we want it
    def f(AT): # that which gives it
        T=mtanh_const_delta_T(Z,T0,AT,T0c,ddx_T_core,ped_width,ped_pos)
        return T(x_LCFS) - T_LCFS
    first_guess = 0.1 # how we view it
    AT = scipy.optimize.fsolve(f,first_guess)[0] # how it is

    T=mtanh_const_delta_T(Z,T0,AT,T0c,ddx_T_core,ped_width,ped_pos)
    ddx_T=mtanh_const_delta_T_ddpsiN(Z,T0,AT,T0c,ddx_T_core,ped_width,ped_pos)
    return (T,ddx_T,AT)
    
def generate_mtanh_const_delta_T(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos,dpsi_ds=1):
    (a_ped,a_sol,a_etb,a_delta,a_slope,b_slope) = parameter_wrapper(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos)
    return (m2tanh_profile(a_ped,a_sol,a_etb,a_delta,dXdr_core,dXdr_sol),lambda x: dpsi_ds*ddx_m2tanh_profile(a_ped,a_sol,a_etb,a_delta,dXdr_core,dXdr_sol)(x))
    
def ped_pos_to_middle(ped_pos,ped_width):
    return ped_pos + ped_width/5

def ped_width_to_trans_width(ped_width):
    return ped_width/5.5

def mtanh_const_delta_T(Z,T0,AT,T0_core,dTdx_core,ped_width,ped_pos):
    BT = AT/(2*sqrt(T0))
    p = ped_pos_to_middle(ped_pos,ped_width)
    w = ped_width_to_trans_width(ped_width)
    T_SOL = constant_delta_T(T0,ped_pos,Z,BT)
    T_core = lambda psiN: T0_core + dTdx_core*(psiN-ped_pos)
    return mtanh_transition(T_SOL,T_core,w,p)

def mtanh_const_delta_T_ddpsiN(Z,T0,AT,T0_core,dTdx_core,ped_width,ped_pos):
    BT = AT/(2*sqrt(T0))
    p = ped_pos_to_middle(ped_pos,ped_width)
    w = ped_width_to_trans_width(ped_width)
    T_SOL = constant_delta_T(T0,ped_pos,Z,BT)
    dT_SOL_dpsiN = constant_delta_dTdpsiN(T0,ped_pos,Z,BT)
    T_core = lambda psiN: T0_core + dTdx_core*(psiN-ped_pos)
    dT_core_dpsiN = lambda psiN:  dTdx_core + 0.0*psiN
    return ddx_mtanh_transition(T_SOL,T_core,dT_SOL_dpsiN,dT_core_dpsiN,w, p)

def mtanh_const_delta_X(Z,eta0,deltaEta,T0,AT,deltaT,eta0_core,detadx_core,ped_width,ped_pos):
    BT = AT/(2*sqrt(T0))
    p = ped_pos_to_middle(ped_pos,ped_width)
    w = ped_width_to_trans_width(ped_width)
    eta_SOL = constant_delta_X(eta0,ped_pos,Z,deltaEta,BT,deltaT)
    eta_core = lambda psiN: eta0_core + detadx_core*(psiN-ped_pos)
    return mtanh_transition(eta_SOL,eta_core,w,p)

def mtanh_const_delta_X_ddpsiN(Z,eta0,deltaEta,T0,AT,deltaT,eta0_core,detadx_core,ped_width,ped_pos):
    BT = AT/(2*sqrt(T0))
    p = ped_pos_to_middle(ped_pos,ped_width)
    w = ped_width_to_trans_width(ped_width)
    eta_SOL = constant_delta_X(eta0,ped_pos,Z,deltaEta,BT,deltaT)
    deta_SOL_dpsiN = constant_delta_dXdpsiN(eta0,ped_pos,Z,deltaEta,BT,deltaT)
    eta_core = lambda psiN: eta0_core + detadx_core*(psiN-ped_pos)
    deta_core_dpsiN = lambda psiN: detadx_core +0.0*psiN
    return ddx_mtanh_transition(eta_SOL,eta_core,deta_SOL_dpsiN,deta_core_dpsiN,w, p)    


if __name__ == "__main__":

    
    old=False
    pure_AUG=False
    if pure_AUG==False:
        #RBar = 1.51486662398
        #BBar = 3.00850351298
        RHat = 1 #FSA(RHat_ped), RBar=1m
        Delta = 0.0014188872143349514  #
        psiAHat = 0.024908932380078564  # based on John's ASDEX work. RBar=1.51486662398 m, BBar=3.00850351298 T
        
        mHati = 1.0
        mHate = 0.000271933691768
        Zi=1
        Ze=-1
 
        T0ce = 0.42
        T0e=T0ce
        Te_LCFS = 0.1
        ddx_Te_core = -1.4
        T0ci = 0.42
        T0i=T0ci
        Ti_LCFS = 0.35
        ddx_Ti_core = -1.2

        psiN0 = 0.94776 # based on AUG #27963 T_e
        w = 0.05224 # based on AUG #27963 T_e

        n0 = 0.70 # 10^-20 m^{-3}
        dndpsiN_core = -0.40616 # 10^-20 m^{-3}
        n_LCFS = 0.12 # 10^-20 m^{-3}
    else:
        #RBar = 1
        #BBar = 1
        RHat = 1.51486662398 #FSA(RHat_ped), RBar=1m
        Delta = 0.0064665523149661508  #
        psiAHat = 0.113522  # based on John's ASDEX work. RBar=1m, BBar=1T
        
        mHati = 1.0
        mHate = 0.000271933691768
        Zi=1
        Ze=-1
        
        #deltaTi = 0.1
        #deltaTe = 0.005
        #ATi = deltaTi *psiAHat/(sqrt(mHati)*RHat*Delta)
        #ATe = deltaTe *psiAHat/(sqrt(mHate)*RHat*Delta)
        
        T0ce = 0.42
        T0e=T0ce
        Te_LCFS = 0.1
        ddx_Te_core = -1.4
        T0ci = 0.38
        T0i=T0ci
        Ti_LCFS = 0.18
        ddx_Ti_core = -1.4

        psiN0 = 0.94776 # based on AUG #27963 T_e
        w = 0.05224 # based on AUG #27963 T_e

        n0 = 0.70 # 10^-20 m^{-3}
        dndpsiN_core = -0.40616 # 10^-20 m^{-3}
        n_LCFS = 0.12 # 10^-20 m^{-3}

    (Ti,ddx_Ti,ATi) = T_parameter_wrapper(Zi,T0i,Ti_LCFS,ddx_Ti_core,w,psiN0)
    (Te,ddx_Te,ATe) = T_parameter_wrapper(Ze,T0e,Te_LCFS,ddx_Te_core,w,psiN0)
    deltaTi = ATi*(sqrt(mHati)*RHat*Delta)/psiAHat 
    deltaTe = ATe*(sqrt(mHate)*RHat*Delta)/psiAHat 

    #Ti=mtanh_const_delta_T(Zi,T0i,ATi,T0ci,dTidx_core,w,psiN0)
    #Te=mtanh_const_delta_T(Ze,T0e,ATe,T0ce,dTedx_core,w,psiN0)
    #ddx_Ti=mtanh_const_delta_T_ddpsiN(Zi,T0i,ATi,T0ci,ddx_Ti_core,w,psiN0)
    #ddx_Te=mtanh_const_delta_T_ddpsiN(Ze,T0e,ATe,T0ce,ddx_Te_core,w,psiN0)

    (etae,etai,ddx_etae,ddx_etai) = eta_parameter_wrapper(Zi,Ze,mHati,mHate,T0i,T0e,deltaTi,deltaTe,ATi,ATe,Ti,Te,n0,dndpsiN_core,n_LCFS,w,psiN0)
    Phi = single_ion_all_eta_Phi(etai,etae,Ti,Te,Zi)
    ddx_Phi = single_ion_all_eta_dPhidpsiN(etai,etae,Ti,Te,ddx_etai,ddx_etae,ddx_Ti,ddx_Te,Zi)
    ne = lambda psiN:  etae(psiN)*exp(-Ze*Phi(psiN)/Te(psiN))
    ni = lambda psiN:  etai(psiN)*exp(-Zi*Phi(psiN)/Ti(psiN))

    psiN = numpy.linspace(0.64,1.0,200)
    D=diff_matrix(psiN[0],psiN[-1],len(psiN))
    plt.subplot(3, 1, 1)
    plt.plot(psiN,ni(psiN),color="red")
    plt.plot(psiN,ne(psiN),color="darkblue")
    plt.ylabel(r"$n$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1])

    plt.subplot(3, 1, 2)
    num_dnidpsiN = numpy.dot(D,ni(psiN))
    num_dnedpsiN = numpy.dot(D,ne(psiN))
    plt.plot(psiN,num_dnidpsiN,color="red")
    plt.plot(psiN,num_dnedpsiN,color="darkblue")
    plt.ylabel(r"$dn/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-30,0])

    plt.subplot(3, 1, 3)
    deltaNi = -Delta/(numpy.fabs(Zi)*psiAHat)*RHat*sqrt(mHati*Ti(psiN))*num_dnidpsiN/ni(psiN)
    mHat_e = 0.000271933691768
    deltaNe = -Delta/(numpy.fabs(Ze)*psiAHat)*RHat*sqrt(mHate*Te(psiN))*num_dnedpsiN/ne(psiN)
    plt.plot(psiN,deltaNi,color="red")
    plt.plot(psiN,deltaNe,color="darkblue")
    #plt.ylim([0,0.2])
    plt.ylabel(r"$\delta_n$")
    plt.xlabel(r"$\psi_N$")

    plt.savefig("mtanh_constDeltaN.pdf")
    plt.show()

    plt.subplot(2, 1, 1)
    plt.plot(psiN,Phi(psiN),color='k')
    plt.ylabel(r"$\Phi$")
    plt.xlabel(r"$\psi_N$")

    plt.subplot(2, 1, 2)
    num_ddx_Phi = numpy.dot(D,Phi(psiN))
    plt.plot(psiN,ddx_Phi(psiN),color='k')
    plt.plot(psiN,num_ddx_Phi,color='r',ls=":") # overlaps with analytical
    plt.ylabel(r"$d\Phi/d\psi_N$")
    plt.xlabel(r"$\psi_N$")

    plt.savefig("mtanh_constDeltaPhi.pdf")
    plt.show()

    plt.subplot(3, 1, 1)
    plt.plot(psiN,etai(psiN),color="red")
    plt.plot(psiN,etae(psiN),color="darkblue")


    plt.ylabel(r"$\eta$")
    plt.xlabel(r"$\psi_N$")
    #plt.ylim([0,1.2])
    plt.ylim([0,2])

    plt.subplot(3, 1, 2)
    plt.plot(psiN,ddx_etai(psiN),color="red")
    plt.plot(psiN,ddx_etae(psiN),color="darkblue")
    plt.ylabel(r"$d\eta/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-40,0])

    plt.subplot(3, 1, 3)
    deltaEtai = -Delta/(Zi*psiAHat)*RHat*sqrt(mHati*Ti(psiN))*ddx_etai(psiN)/etai(psiN)
    deltaEtae = -Delta/(numpy.fabs(Ze)*psiAHat)*RHat*sqrt(mHate*Te(psiN))*ddx_etae(psiN)/etae(psiN)


    plt.plot(psiN,deltaEtai,color="red")
    plt.plot(psiN,deltaEtae,color="darkblue")
    plt.ylabel(r"$\delta_\eta$")
    plt.xlabel(r"$\psi_N$")
    plt.savefig("mtanh_constDeltaEta.pdf")
    plt.show()

    plt.subplot(3, 1, 1)
    plt.plot(psiN,Ti(psiN),color="red")
    plt.plot(psiN,Te(psiN),color="darkblue")
    plt.ylabel(r"$T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1])

    plt.subplot(3, 1, 2)
    plt.plot(psiN,ddx_Ti(psiN),color="red")
    plt.plot(psiN,ddx_Te(psiN),color="darkblue")
    plt.ylabel(r"$dT/d\psi_N$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([-15,0])

    plt.subplot(3, 1, 3)
    deltaTi = -Delta/(Zi*psiAHat)*RHat*sqrt(mHati)*ddx_Ti(psiN)/sqrt(Ti(psiN))
    deltaTe = -Delta/(Zi*psiAHat)*RHat*sqrt(mHate)*ddx_Te(psiN)/sqrt(Te(psiN))
    plt.plot(psiN,deltaTi,color="red")
    plt.plot(psiN,deltaTe,color="darkblue")    
    plt.ylabel(r"$\delta_T$")
    plt.xlabel(r"$\psi_N$")
    plt.ylim([0,1])
    plt.savefig("mtanh_constDeltaT.pdf")
    plt.show()

    if old:
        only_mtanh=False
        ##### TUNEABLES ##############

        deltaTi = 0.1
        deltaEtai = 0.1
        deltaTe = 0.005
        deltaEtae = 0.05

        eta0c = 1.0
        T0c = 1.0
        T0=T0c
        eta0=eta0c

        psiN0 = 0.9 # based on He paper

        w = 0.1/2
        p = psiN0 + w


        if only_mtanh:
            #Make true to get ordinary mtanh profiles
            deltaTi = 0.00000000000001
            deltaEtai = 0.00000000000001
            deltaTe = 0.00000000000001
            deltaEtae = 0.00000000000001
            w = 0.05
            p = psiN0 + w/2
            eta0c = 1.0
            T0c = 1.0
            T0=T0c/5
            eta0=eta0c/5



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
        full_dPhidpsiN = single_ion_all_eta_dPhidpsiN(full_etai,full_etae,full_Ti,full_Te,full_detaidpsiN,full_detaedpsiN,full_dTidpsiN,full_dTedpsiN,Zi)

        # densities
        full_ni = lambda psiN: full_etai(psiN)*exp(-Zi*full_Phi(psiN)/full_Ti(psiN))
        full_ne = lambda psiN: full_etae(psiN)*exp(-Ze*full_Phi(psiN)/full_Te(psiN))


        ############# PLOT ETA ####################
        plt.subplot(3, 1, 1)
        plt.plot(psiN,full_etai(psiN))
        plt.plot(psiN,full_etae(psiN))
        plt.plot(psiN,etai(psiN),ls=":",color="darkblue")
        plt.plot(psiN,etae(psiN),ls=":",color="darkgreen")


        etai2=mtanh_const_delta_X(Zi,eta0,deltaEtai,T0i,ATi,deltaTi,eta0c,-a,w,psiN0)
        etae2=mtanh_const_delta_X(Ze,eta0,deltaEtae,T0e,ATe,deltaTe,eta0c,-a,w,psiN0)

        #etae2=mtanh_const_delta_X(Ze,eta0,deltaEtae,ATe,deltaTe,eta0c,-a/100.0,w,psiN0)
        plt.plot(psiN,etai2(psiN),ls="--",color="darkred")
        plt.plot(psiN,etae2(psiN),ls="--",color="magenta")

        plt.ylabel(r"$\eta$")
        plt.xlabel(r"$\psi_N$")
        #plt.ylim([0,1.2])
        plt.ylim([0,2])


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

        etai2=mtanh_const_delta_X_ddpsiN(Zi,eta0,deltaEtai,T0i,ATi,deltaTi,eta0c,-a,w,psiN0)
        etae2=mtanh_const_delta_X_ddpsiN(Ze,eta0,deltaEtae,T0e,ATe,deltaTe,eta0c,-a,w,psiN0)
        plt.plot(psiN,etai2(psiN),ls="--",color="darkred")
        plt.plot(psiN,etae2(psiN),ls="--",color="magenta")

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

        ###### PLOT T #################

        plt.subplot(3, 1, 1)
        plt.plot(psiN,full_Ti(psiN))
        plt.plot(psiN,full_Te(psiN))
        plt.plot(psiN,Ti(psiN),ls=":",color="darkblue")
        plt.plot(psiN,Te(psiN),ls=":",color="darkgreen")

        Ti2=mtanh_const_delta_T(Zi,T0,ATi,T0c,-a,w,psiN0)
        Te2=mtanh_const_delta_T(Ze,T0,ATe,T0c,-a,w,psiN0)
        plt.plot(psiN,Ti2(psiN),ls="--",color="darkred")
        plt.plot(psiN,Te2(psiN),ls="--",color="magenta")

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

        Ti2=mtanh_const_delta_T_ddpsiN(Zi,T0,ATi,T0c,-a,w,psiN0)
        Te2=mtanh_const_delta_T_ddpsiN(Ze,T0,ATe,T0c,-a,w,psiN0)
        plt.plot(psiN,Ti2(psiN),ls="--",color="darkred")
        plt.plot(psiN,Te2(psiN),ls="--",color="magenta")

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




        ############# PLOT Phi ####################
        plt.subplot(2, 1, 1)
        plt.plot(psiN,full_Phi(psiN))
        plt.ylabel(r"$\Phi$")
        plt.xlabel(r"$\psi_N$")

        plt.subplot(2, 1, 2)
        num_full_dPhidpsiN = numpy.dot(D,full_Phi(psiN))
        plt.plot(psiN,full_dPhidpsiN(psiN))
        #plt.plot(psiN,num_full_dPhidpsiN) # overlaps with analytical
        plt.ylabel(r"$d\Phi/d\psi_N$")
        plt.xlabel(r"$\psi_N$")

        plt.savefig("mtanh_constDeltaPhi.pdf")
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

