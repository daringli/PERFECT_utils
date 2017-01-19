import scipy.optimize
from numpy import exp

from mtanh import mtanh_profile, ddx_mtanh_profile
from mtanh_const_delta import single_ion_all_eta_Phi

from mtanh import mtanh_transition, ddx_mtanh_transition
from const_delta import constant_delta_X, constant_delta_dXdpsiN


def n_mtanh_const_delta_from_n_LCFS(n_ped,n_LCFS,ddx_n_core,ped_width,ped_pos,Z,BT,deltaT):
    # fit a transition point (TP) to get the n_LCFS we want
    n_c = lambda psiN: n_ped + ddx_n_core * (psiN - ped_pos)
    ddx_n_c = lambda psiN: ddx_n_core + 0.0 * psiN
    x_LCFS = ped_pos + ped_width

    deltaN = 0.10
    n_SOLp = n_LCFS*(4./7.0) #based on AUG #27963
    x_SOLp = 1.06 #based on AUG #27963
    #TW = ped_width/2.0 #transition width
    def f(TP):
        # create a constant delta SOL profile
        TW= TP - ped_pos
        

        n_SOL = constant_delta_X(n_SOLp,x_SOLp,Z,deltaN,BT,deltaT)
        
        # create transition
        n = mtanh_transition(n_SOL,n_c,TW,TP)
        return n(x_LCFS) - n_LCFS
    TP0 = x_LCFS
    TP = scipy.optimize.fsolve(f,TP0)[0]
    TW= TP - ped_pos
    print "TP: " + str(TP)
    n_SOL = constant_delta_X(n_SOLp,x_SOLp,Z,deltaN,BT,deltaT)
    ddx_n_SOL = constant_delta_dXdpsiN(n_SOLp,x_SOLp,Z,deltaN,BT,deltaT)
    n = mtanh_transition(n_SOL,n_c,TW,TP)
    ddx_n = ddx_mtanh_transition(n_SOL,n_c,ddx_n_SOL,ddx_n_c,TW,TP)
    return (n,ddx_n,n_SOL,ddx_n_SOL)
    
def n_mtanh_from_n_LCFS(n_ped,n_LCFS,ddx_n_core,ped_width,ped_pos):
    a_delta = ped_width/4.0
    a_ped = n_ped + ddx_n_core * ped_width/2.0
    a_etb = ped_pos + ped_width/2.0
    x_LCFS = ped_pos + ped_width
    

    # fit a_sol to get the n_LCFS we want for
    # our mtanh density profile
    def f(a_sol):
        a_slope = -4 * a_delta * ddx_n_core / (a_ped - a_sol)
        n = mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
        return n(x_LCFS) - n_LCFS
    a_sol0=n_LCFS
    a_sol = scipy.optimize.fsolve(f,a_sol0)[0]
    print "a_sol:" + str(a_sol)
    a_slope = -4 * a_delta * ddx_n_core / (a_ped - a_sol)
    n = mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    ddx_n = ddx_mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    return (n,ddx_n)

if __name__=="__main__":
    import matplotlib.pyplot as plt
    import numpy
    n_ped = 1.0
    n_LCFS = n_ped/8.0
    ddx_n_core = -0.1
    ped_width = 0.06
    ped_pos = 0.94
    (n,ddx_n) = n_mtanh_from_n_LCFS(n_ped,n_LCFS,ddx_n_core,ped_width,ped_pos)

    Z=1
    deltaT=0.123142398652 #from AUG #27963 #27963
    BT=1.66786235477 #from AUG
    (n2,ddx_n2) = n_mtanh_const_delta_from_n_LCFS(n_ped,n_LCFS,ddx_n_core,ped_width,ped_pos,Z,BT,deltaT)
    x=numpy.linspace(0.8,1.2)
    plt.subplot(3,1,1)
    plt.plot(numpy.sqrt(x),n(x))
    plt.plot(numpy.sqrt(x),n2(x))
    plt.subplot(3,1,2)
    plt.plot(numpy.sqrt(x),ddx_n(x))
    plt.plot(numpy.sqrt(x),ddx_n2(x))
    plt.subplot(3,1,3)
    plt.plot(numpy.sqrt(x),ddx_n(x)/n(x))
    plt.plot(numpy.sqrt(x),ddx_n2(x)/n2(x))
    plt.show()
