from __future__ import division

import numpy as np
from numpy import exp
from diff_matrix import diff_matrix
import matplotlib.pyplot as plt

from bezier_transition import derivative_bezier_transition
import scipy.optimize


def mtanh_wikipedia(a,b,c,d):
    #returns a modified tanh function with parameters a,b,c,d
    # https://en.wikipedia.org/wiki/Modified_hyperbolic_tangent
    # NOTE: not the actual mtanh used in the fusion community???
    return lambda x: (exp(a*x) - exp(-b*x))/(exp(c*x)+exp(-d*x))

def mtanh(a):
    # From:
    # Deconvolution of Thomson scattering temperature profiles
    # R. Scannell, et al.
    # Rev. Sci. Instrum. 82, 053501 (2011)
    # http://dx.doi.org/10.1063/1.3581230
    return lambda r: ((1+a*r)*exp(r)-exp(-r))/(exp(r) + exp(-r))

def m2tanh(a,b):
    return lambda r: (exp(r) - exp(-r))/(exp(r) + exp(-r)) + (a*r*exp(r) + b*r*exp(-r))/(exp(r) + exp(-r))

def m2tanh_old(a,b):
    return lambda r: ((1+a*r)*exp(r)-(1-b*r)*exp(-r))/(exp(r) + exp(-r))

def ddx_mtanh(a):
    #return lambda r: (a*r*exp(2*r) + 2+3*a*r + 2*exp(-2*r))*((exp(r) + exp(-r))**(-2))
    return lambda r: (exp(2*r)/((exp(2*r)+1)**2))*(a*(2*r + exp(2*r) +1) + 4)

def ddx_m2tanh_old(a,b):
    return lambda r: (exp(2*r)*(a*(2*r+exp(2*r)+1)+4)+b*(exp(2*r)*(1-2*r)+1))/(exp(2*r)+1)**2

def ddx_m2tanh(a,b):
    return lambda r: 4/((exp(r) + exp(-r))**2) + (a*(exp(2*r) + 1 + 2*r) + b*(exp(-2*r) + 1 - 2*r))/((exp(r) + exp(-r))**2)

def mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope):
    # From:
    # Deconvolution of Thomson scattering temperature profiles
    # R. Scannell, et al.
    # Rev. Sci. Instrum. 82, 053501 (2011)
    # http://dx.doi.org/10.1063/1.3581230
    return lambda r: a_sol + ((a_ped - a_sol)/2)*(1 + mtanh(a_slope)((a_etb-r)/(2*a_delta)))

def m2tanh_profile(a_ped,a_sol,a_etb,a_delta,ddx_P_core,ddx_P_sol):
    return lambda r: a_sol + (((a_ped - a_sol)/2) + ((a_ped - a_sol)/2)*(exp(((a_etb-r)/(2*a_delta))) - exp(-((a_etb-r)/(2*a_delta))))/(exp(((a_etb-r)/(2*a_delta))) + exp(-((a_etb-r)/(2*a_delta))))
                              - (ddx_P_core*(a_etb-r)*exp(((a_etb-r)/(2*a_delta))) + ddx_P_sol*(a_etb-r)*exp(-((a_etb-r)/(2*a_delta))))/(exp(((a_etb-r)/(2*a_delta))) + exp(-((a_etb-r)/(2*a_delta)))))

                                                  
def m2tanh_profile_old(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope):
    return lambda r: a_sol + ((a_ped - a_sol)/2)*(1 + m2tanh(a_slope,b_slope)((a_etb-r)/(2*a_delta)))

def ddx_mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope):
    return lambda r: -1/(4*a_delta)*(a_ped - a_sol)*ddx_mtanh(a_slope)((a_etb-r)/(2*a_delta))

def ddx_m2tanh_profile_old(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope):
    return lambda r: -1/(4*a_delta)*(a_ped - a_sol)*ddx_m2tanh_old(a_slope,b_slope)((a_etb-r)/(2*a_delta))

def ddx_m2tanh_profile(a_ped,a_sol,a_etb,a_delta,ddx_P_core,ddx_P_sol):
    #return lambda r: -1/(4*a_delta)*(a_ped - a_sol)*ddx_m2tanh(a_slope,b_slope)
    #a_slope = -ddx_P_core *(4*a_delta)/(a_ped - a_sol)
    #b_slope = -ddx_P_sol *(4*a_delta)/(a_ped - a_sol)
    return lambda r: -(a_ped - a_sol)/(4*a_delta)*(4/((exp(((a_etb-r)/(2*a_delta)))+ exp(-((a_etb-r)/(2*a_delta))))**2)) + ((ddx_P_core*(exp(2*((a_etb-r)/(2*a_delta))) + 1 + 2*((a_etb-r)/(2*a_delta))) + ddx_P_sol*(exp(-2*((a_etb-r)/(2*a_delta))) + 1 - 2*((a_etb-r)/(2*a_delta))))/((exp(((a_etb-r)/(2*a_delta))) + exp(-((a_etb-r)/(2*a_delta))))**2))
    
def parameter_wrapper(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos):
    #takes parameters from spline profile generation and maps them to
    #modified modified arctangens hyperbolicus profiles
    a_delta = ped_width/4.0
    a_ped = X_ped + dXdr_core * ped_width/2.0
    X_sol = X_ped + dXdr_ped * ped_width
    a_sol = X_ped + (dXdr_ped  - dXdr_sol/2.0) * ped_width
    print "a_ped: " + repr(a_ped)
    print "a_sol: " + repr(a_sol)
    try:
        a_slope = -4 * a_delta * dXdr_core / (a_ped - a_sol)
        b_slope = -4 * a_delta * dXdr_sol / (a_ped - a_sol)
    except ZeroDivisionError:
        print "parameter_wrapper: WARNING: a_ped = a_sol, thus a_slope and b_slope are undefined. This should NOT be a problem for the new profile generation script"
        a_slope = np.divide(-4 * a_delta * dXdr_core,(a_ped - a_sol))
        b_slope = np.divide(-4 * a_delta * dXdr_sol,(a_ped - a_sol))
    a_etb = ped_pos + ped_width/2.0
    return (a_ped,a_sol,a_etb,a_delta,a_slope,b_slope)

def generate_m2tanh_profile(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos,dpsi_ds=1):
    (a_ped,a_sol,a_etb,a_delta,a_slope,b_slope) = parameter_wrapper(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos)
    return (m2tanh_profile(a_ped,a_sol,a_etb,a_delta,dXdr_core,dXdr_sol),lambda x: dpsi_ds*ddx_m2tanh_profile(a_ped,a_sol,a_etb,a_delta,dXdr_core,dXdr_sol)(x))

def extrapolate_m2tanh_sections(P,dPdx,x_core_stop,x_ped_start,x_ped_stop,x_sol_start,dpsi_ds=1):
    #for some parts of the profile generation, we need to extrapolate certain
    #segments of the profiles
    #this creates linearly extrapolated functions for the core, pedestal
    #and buffer region
    # x_core_stop: x where core stops
    # x_ped_start : x where pedestal starts
    # x_ped_stop : x where pedestal stops
    # x_sol_start: x where sol starts
    
    slope_core_stop = dpsi_ds * dPdx(x_core_stop)
    slope_ped_start = dpsi_ds * dPdx(x_ped_start)
    slope_ped_stop = dpsi_ds * dPdx(x_ped_stop)
    slope_sol_start = dpsi_ds * dPdx(x_sol_start)
    
    value_core_stop = P(x_core_stop)
    value_ped_start = P(x_ped_start)
    value_ped_stop = P(x_ped_stop)
    value_sol_start = P(x_sol_start)

    #create extrapolations
    after_core = lambda x: value_core_stop + slope_core_stop*(x-x_core_stop)
    before_ped = lambda x: value_ped_start + slope_ped_start*(x-x_ped_start)
    after_ped = lambda x: value_ped_stop + slope_ped_stop*(x-x_ped_stop)
    before_sol = lambda x: value_sol_start + slope_sol_start*(x-x_sol_start)

    ddx_after_core = lambda x: slope_core_stop
    ddx_before_ped = lambda x: slope_ped_start
    ddx_after_ped = lambda x: slope_ped_stop
    ddx_before_sol = lambda x: slope_sol_start

    #NOTE: SMOOTHNESS OF TRANSITION A PROBLEM.
    # mtanh transition function
    pair_list = [[0,0]] 
    core_funclist = [P,after_core]
    core_derivlist = [dPdx,ddx_after_core]
    core_pointlist = [x_core_stop]
    (core,ddx_core) = derivative_bezier_transition(core_funclist,core_derivlist,core_pointlist,pair_list)

    pair_list = [[0,0],[0,0]]
    ped_funclist = [before_ped,P,after_ped]
    ped_derivlist = [ddx_before_ped,dPdx,ddx_after_ped]
    ped_pointlist = [x_ped_start,x_ped_stop]
    (ped,ddx_ped) = derivative_bezier_transition(ped_funclist,ped_derivlist,ped_pointlist,pair_list)

    pair_list = [[0,0]]
    sol_funclist = [before_sol,P]
    sol_derivlist = [ddx_before_sol,dPdx]
    sol_pointlist = [x_sol_start]
    (sol,ddx_sol) = derivative_bezier_transition(sol_funclist,sol_derivlist,sol_pointlist,pair_list)

    return (core,ped,sol,ddx_core,ddx_ped,ddx_sol)

def match_heat_flux_proxy(T_ped,dTdr_core_0,dTdr_ped,dTdr_sol,ped_width,ped_pos,n_a,n_b,r_a,r_b):
    # find a dTdr_core to match heat flux proxy
    # at r_a and r_b, assuming
    # dTdr_core_0: initial guess
    def f(a):
        (T,dTdr)=generate_m2tanh_profile(T_ped,a,dTdr_ped,dTdr_sol,ped_width,ped_pos)
        (T_a,dTdr_a,T_b,dTdr_b) = (T(r_a),dTdr(r_a),T(r_b),dTdr(r_b))
        #n_a = n(r_a)
        #n_b = n(r_b)
        return T_a**(3./2.)*dTdr_a*n_a - T_b**(3./2.)*dTdr_b*n_b
    dTdr_core=scipy.optimize.fsolve(f,dTdr_core_0)[0]
    return generate_m2tanh_profile(T_ped,dTdr_core,dTdr_ped,dTdr_sol,ped_width,ped_pos)

def mtanh_transition(f,g,a=1):
    return lambda r: (f(r)*exp(r/a) + g(r)*exp(-r/a))/(exp(r/a)+exp(-r/a))

def ddx_mtanh_transition(f,g,ddx_f,ddx_g,a=1):
    return lambda r: ((2/a)*(f(r) - g(r)) + ddx_f(r)*(1+exp(2*r/a)) + ddx_g(r)*(1+exp(-2*r/a)))/((exp(r/a)+exp(-r/a))**2.0)


def generate_m3tanh_profile(a_ped,a_sol,a_etb,dXdr_core,dXdr_sol,a_delta,dpsi_ds=lambda r : 1 + 0*r):
    a=2*a_delta
    if a_ped == a_sol:
        print "a_ped = a_sol"
        f = lambda r : a_sol - dXdr_sol *(r - a_etb)
        ddx_f = lambda r : - dXdr_sol + 0*r
        g = lambda r : a_ped - dXdr_core *(r - a_etb)
        ddx_g = lambda r : - dXdr_core + 0*r
    return (mtanh_transition(f,g,a),lambda x: dpsi_ds(x)*ddx_mtanh_transition(f,g,ddx_f,ddx_g,a)(x))


if __name__ == "__main__":
    (a_ped,a_sol,a_etb,a_delta,a_slope) = (2,1,0,1,0.1)
    P=mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    dPdx = ddx_mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)

    xlo = -10
    xhi = 10
    xlims = [xlo,xhi] 
    x=np.linspace(xlo,xhi,100)
    D=diff_matrix(x[0],x[-1],len(x))

    #plt.subplot(2, 1, 1)
    #plt.plot(x, P(x))
    #plt.ylabel(r'$P$')
    #plt.xlim(xlims)

    #plt.subplot(2, 1, 2)
    #plt.plot(x, dPdx(x))
    #plt.hold(True)
    #plt.plot(x, np.dot(D,P(x)))
    #plt.ylabel(r'$dP/dx$')
    #plt.xlim(xlims)

    #plt.show()

    #(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope) = (2,1,0,1,0.1,0)
    #P=m2tanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope)
    #dPdx = ddx_m2tanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope)
    #(core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(P,dPdx,-2,-2,2,2)
    xlo = 0.84
    xhi = 1.11
    xlims = [xlo,xhi] 
    x=np.linspace(xlo,xhi,100)

    psiMinPed = 0.94927395957
    psiMaxPed = psiMinPed + 0.0338173602865
    #Tpeds[e_index],TCoreGrads[e_index],TpedGrads[e_index],TSOLGrads[e_index],psiMaxPed-psiMinPed,psiMinPed
    X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos =0.9, -1.77423664921, -17.7423664921, -1.77423664921, 0.0338173602865, 0.94927395957
    P,dPdx = generate_m2tanh_profile(X_ped,dXdr_core,dXdr_ped,dXdr_sol,ped_width,ped_pos)
    (core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(P,dPdx,xlo,ped_pos+ped_width/2,ped_pos+ped_width/2,xhi)
    a_ped = X_ped + dXdr_core * ped_width/2.0
    a_sol = X_ped + (dXdr_ped  - dXdr_sol/2.0) * ped_width
    X_sol = X_ped + dXdr_ped * ped_width
    #ped = lambda x: X_ped + dXdr_ped * (x - ped_pos)

    flist = [core,ped,sol]
    ddx_flist = [ddx_core,ddx_ped,ddx_sol]
    pointlist = [psiMinPed,psiMaxPed]
    offset = (psiMaxPed-psiMinPed)*0.2
    pairList=[[offset,offset],[offset,offset]]
    P2,dP2dx = derivative_bezier_transition(flist,ddx_flist,pointlist,pairList)
    THatPre =(lambda psiN: (X_ped + dXdr_core*(psiN-psiMinPed)))
    THatPed =(lambda psiN: (X_ped + dXdr_ped*(psiN-psiMinPed)))
    THatAft =(lambda psiN: (X_ped + dXdr_ped*(psiMaxPed-psiMinPed) + dXdr_sol*(psiN-psiMaxPed)))
    dTHatPredpsi = (lambda psiN: dXdr_core + 0*psiN)
    dTHatPeddpsi = (lambda psiN: dXdr_ped + 0*psiN)
    dTHatAftdpsi = (lambda psiN: dXdr_sol + 0*psiN)
    Tlist=[THatPre,THatPed,THatAft]
    dTdpsiList = [dTHatPredpsi,dTHatPeddpsi,dTHatAftdpsi]
    P3,dP3dx = derivative_bezier_transition(Tlist,dTdpsiList,pointlist,pairList)

    plt.title("linear extrapolation")
    plt.subplot(2, 1, 1)

    ax=plt.gca()
    ax.axvline(x=ped_pos+ped_width/2,color='k',linestyle=':')
    ax.axhline(y=a_ped,color='k',linestyle=':')
    ax.axhline(y=a_sol,color='k',linestyle=':')
    
    plt.plot(x, P(x))
    plt.hold(True)
    plt.plot(x, P2(x))
    plt.plot(x, P3(x))
    plt.plot(x, core(x))
    plt.plot(x, sol(x))
    plt.plot(x, ped(x))
    plt.ylabel(r'$P$')
    plt.xlim(xlims)
    plt.ylim([0,1.1])
    ax.plot(ped_pos,X_ped,'o')
    ax.plot(ped_pos+ped_width,X_sol,'o')

    plt.subplot(2, 1, 2)
    plt.plot(x, dPdx(x))
    plt.hold(True)
    plt.plot(x, dP2dx(x))
    plt.plot(x, dP3dx(x))
    #plt.plot(x, np.dot(D,P(x)))
    plt.plot(x, ddx_core(x))
    #plt.plot(x, ddx_ped(x))
    plt.plot(x, ddx_sol(x))
    plt.ylabel(r'$dP/dx$')
    plt.xlim(xlims)
    plt.show()

    plt.title("heat flux proxy matching")
    xlo = -10
    xhi = 10
    xlims = [xlo,xhi] 
    x=np.linspace(xlo,xhi,100)
    
    n = lambda x: 6.0 - (3.0/20)*(x + 10.0)
    print "n_a: " + str(n(xlo))
    print "n_b: " + str(n(xhi))
    (P,dPdx) = match_heat_flux_proxy(4,-0.1,-0.1,-0.2,1,0,n(xlo),n(xhi),xlo,xhi)
    
    plt.subplot(3, 1, 1)
    plt.plot(x, P(x))
    plt.xlim(xlims)
    
    plt.subplot(3, 1, 2)
    plt.plot(x, dPdx(x))
    plt.xlim(xlims)


    Q = lambda x : n(x)*P(x)**(3.0/2.0)*dPdx(x)
    plt.subplot(3, 1, 3)
    plt.plot(x, Q(x))
    plt.xlim(xlims)

    print "Q_a: " + str(Q(xlo))
    print "Q_b: " + str(Q(xhi))
    print "Q_b - Q_a: " + str.format('{0:.5e}', Q(xhi)-Q(xlo))
    
    plt.show()
    

    # plt.title("mtanh transition")
    
    # a = 0.1
    # b = 1.0
    # n = lambda x: 6.0 - (3.0/20)*(x + 10.0)
    # ddx_n = lambda x: - (3.0/20) + 0*x
    # T = lambda x: 3-a*x
    # ddx_T = lambda x: -a + 0*x
    # Phi_c = 8.0
    # f = lambda r: n(r)*exp(Phi_c/T(x))
    # g = lambda r: 1+b*(r - 2)
    # ddx_f = lambda r: (ddx_n(r) - n(r)*Phi_c/(T(r)**2)*ddx_T(r))*exp(Phi_c/T(r))
    # ddx_g = lambda r: b + r*0
    # h = mtanh_transition(f,g)
    # #h2 = m2tanh(a,b)
    # ddx_h = ddx_mtanh_transition(f,g,ddx_f,ddx_g)
    # #ddx_h2 = ddx_m2tanh(a,b)
    
    # plt.subplot(3, 1, 1)
    # plt.hold(True)
    # #plt.plot(x, h2(x))
    # plt.plot(x, h(x))
    # plt.plot(x, f(x))
    # plt.plot(x, g(x))
    # plt.xlim(xlims)

    # plt.subplot(3, 1, 2)
    # plt.hold(True)
    # #plt.plot(x, ddx_h2(x))
    # plt.plot(x, ddx_h(x))
    # plt.plot(x, ddx_f(x))
    # plt.plot(x, ddx_g(x))
    # plt.xlim(xlims)
    
    # plt.show()

    plt.title("old vs new")

    a,b = 1,2
    P1 = m2tanh_old(a,b)
    P2 = m2tanh(a,b)
    ddx_P1 = ddx_m2tanh_old(a,b)
    ddx_P2 = ddx_m2tanh(a,b)

    plt.subplot(4, 1, 1)
    plt.hold(True)
    plt.plot(x, P1(x))
    plt.plot(x, P2(x))
    plt.xlim(xlims)

    plt.subplot(4, 1, 2)
    plt.hold(True)
    plt.plot(x, ddx_P1(x))
    plt.plot(x, ddx_P2(x))
    plt.xlim(xlims)

    a_ped = 2.35
    a_sol = 1.23
    a_etb = 0.31
    a_delta = 0.78
    ddx_P_core = -0.11
    ddx_P_sol = -0.25
    a_slope = -ddx_P_core*4*a_delta/(a_ped - a_sol)
    b_slope = -ddx_P_sol*4*a_delta/(a_ped - a_sol)

    P3 = m2tanh_profile_old(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope)
    P4 = m2tanh_profile(a_ped,a_sol,a_etb,a_delta,ddx_P_core,ddx_P_sol)
    ddx_P3 = ddx_m2tanh_profile_old(a_ped,a_sol,a_etb,a_delta,a_slope,b_slope)
    ddx_P4 = ddx_m2tanh_profile(a_ped,a_sol,a_etb,a_delta,ddx_P_core,ddx_P_sol)

    plt.subplot(4, 1, 3)
    plt.hold(True)
    plt.plot(x, P3(x))
    plt.plot(x, P4(x))
    
    plt.subplot(4, 1, 4)
    plt.hold(True)
    plt.plot(x, ddx_P3(x))
    plt.plot(x, ddx_P4(x))
    
    plt.show()
