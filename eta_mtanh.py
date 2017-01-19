from mtanh import mtanh_profile, ddx_mtanh_profile
from mtanh_const_delta import single_ion_all_eta_Phi
import scipy.optimize
from numpy import exp


def eta_from_n_LCFS(n_ped,n_LCFS,ddx_n_core,Zi,Ze,Ti,Te,ped_width,ped_pos):
    a_delta = ped_width/4.0
    #a_delta = ped_width/20.0
    a_ped = n_ped + ddx_n_core * ped_width/2.0
    a_etb = ped_pos + ped_width/2.0
    x_LCFS = ped_pos + ped_width
    

    # nice and flat ion eta mtanh
    a_soli = n_ped + 3*ddx_n_core * ped_width
    a_slopei = -4 * a_delta * ddx_n_core / (a_ped - a_soli)
    etai = mtanh_profile(a_ped,a_soli,a_etb,a_delta,a_slopei)
    ddx_etai = ddx_mtanh_profile(a_ped,a_soli,a_etb,a_delta,a_slopei)

    # nice and flat ion eta
    #etai = lambda psiN: n_ped + ddx_n_core * (psiN - ped_pos)
    #ddx_etai = lambda psiN: ddx_n_core + 0 * psiN
    a_delta = a_delta/4
    # fit a_sol for the eta_e mtanh profile to get the n_sol we want
    def f(a_sol):
        a_slope = -4 * a_delta * ddx_n_core / (a_ped - a_sol)
        etae = mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
        Phi = single_ion_all_eta_Phi(etai,etae,Ti,Te,Zi)
        ne = lambda psiN:  etae(psiN)*exp(-Ze*Phi(psiN)/Te(psiN))
        ni = lambda psiN:  etai(psiN)*exp(-Zi*Phi(psiN)/Ti(psiN))
        return ne(x_LCFS) - n_LCFS
    a_sol0=0.0
    a_sol = scipy.optimize.fsolve(f,a_sol0)[0]
    print "a_sol:" + str(a_sol)
    a_slope = -4 * a_delta * ddx_n_core / (a_ped - a_sol)
    etae = mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    ddx_etae = ddx_mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    return (etae,etai,ddx_etae,ddx_etai)

