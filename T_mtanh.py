import scipy.optimize
from numpy import exp,sqrt

from mtanh import mtanh_profile, ddx_mtanh_profile
from mtanh_const_delta import single_ion_all_eta_Phi

from mtanh import mtanh_transition, ddx_mtanh_transition
from const_delta import constant_delta_T, constant_delta_dTdpsiN

def T_mtanh_const_delta_from_T_LCFS(T_ped,T_LCFS,ddx_T_core,ped_width,ped_pos,Z,AT):
    # fit a transition point (TP) to get the n_LCFS we want
    T_c = lambda psiN: T_ped + ddx_T_core * (psiN - ped_pos)
    ddx_T_c = lambda psiN: ddx_T_core + 0.0 * psiN
    x_LCFS = ped_pos + ped_width

    BT=AT/(2*sqrt(T_ped))
    T_SOLp = T_LCFS*(1./3.0) #based on AUG #27963
    x_SOLp = 1.06 #based on AUG #27963
    #TW = ped_width/2.0 #transition width
    def f(TP):
        # create a constant delta SOL profile
        TW= TP - ped_pos
        

        T_SOL = constant_delta_T(T_SOLp,x_SOLp,Z,BT)
        
        # create transition
        T = mtanh_transition(T_SOL,T_c,TW,TP)
        return T(x_LCFS) - T_LCFS
    TP0 = x_LCFS
    TP = scipy.optimize.fsolve(f,TP0)[0]
    TW= TP - ped_pos
    print "T TP: " + str(TP)
    T_SOL = constant_delta_T(T_SOLp,x_SOLp,Z,BT)
    ddx_T_SOL = constant_delta_dTdpsiN(T_SOLp,x_SOLp,Z,BT)
    T = mtanh_transition(T_SOL,T_c,TW,TP)
    ddx_T = ddx_mtanh_transition(T_SOL,T_c,ddx_T_SOL,ddx_T_c,TW,TP)

    return (T,ddx_T,AT,T_SOL,ddx_T_SOL)
