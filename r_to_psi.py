import numpy
import matplotlib.pyplot as plt
import scipy.integrate,scipy.misc
from toroidal_flux_q import func_q,func_polyfit,func_polyfit_derivative,toroidal_flux_derivative

def dpsidr(q,dPsiTdr):
    return lambda r: (1/(2*numpy.pi*q(r)))*dPsiTdr(r)

def func_psi(q,dPsiTdr):
    integrand=dpsidr(q,dPsiTdr)
    def integral(x):
        i,err=scipy.integrate.quad(integrand,0,x)
        return i
    return integral
    #print scipy.integrate.quad(integrand,0,0.6)

def dpsiNdr(psi,r):
    from diff_matrix import diff_matrix
    psi_array=numpy.array([psi(ri) for ri in r])
    psiN_array=psi_array/psi_array[-1]
    rmin=r[0]
    rmax=r[-1]
    N=len(r)
    diff_matrix=diff_matrix(rmin,rmax,N)
    dpsiNdr_array=numpy.dot(diff_matrix,psiN_array)
    return dpsiNdr_array

def dpsiNdr_analytic(q,dPsiTdr,psi,r_array):
    return dpsidr(q,dPsiTdr)(r_array)/psi(r_array[-1])

def drdPsiN_at_psiN_one(q,dPsiTdr,psi,r_ped):
    #assumes r_ped is where psiN=1
    #which could be a good assumption since dpsiN/dr decreases at in the pedestal and hence psi_N is not that sensitive to the exact choice of r? (still greater than 1, though...)
    if type(r_ped) != list:
        r_ped=numpy.array([r_ped])
        
    dpsiNdr=dpsiNdr_analytic(q,dPsiTdr,psi,r_ped)
    return 1/dpsiNdr


def psi_pedestal_width(deltar,q,dPsiTdr,psi,r_ped):
    drdpsi=drdPsiN_at_psiN_one(q,dPsiTdr,psi,r_ped)[0]
    return deltar/drdpsi

if __name__=='__main__':
    y=numpy.array([1.23,1.37,1.6])
    x=numpy.array([0,0.5,0.9])
    #flat profile
    #y=numpy.array([1.73,1.73,1.73])
    #x=numpy.array([0,0.5,0.9])
    deg=2
    kappa=func_polyfit(x,y,deg)
    dkappadrho=func_polyfit_derivative(x,y,deg)
    q_of_rho=func_q(1.0,3.5,0.6,1.6)
    
    #flat profile
    #q=func_q(3.5,3.5,0.6,3.5)
    #does not actually matter in psiN since its normalized away
    #would matter if we calculated psi(r)
    BT=lambda x: 2.93 #Value for the average BT, taken from BBar in T for JETLike profile
    dBTdr=lambda x: 0 
    
    
    N=100
    rmin=0
    rmax=0.7
    q=lambda r:q_of_rho(r/rmax)
    dPsiTdr=toroidal_flux_derivative(rmax,kappa,dkappadrho,BT,dBTdr)
    psi=func_psi(q,dPsiTdr)
    r=numpy.linspace(rmin,rmax,N)    
    dpsiN_array=dpsiNdr(psi,r)
    dpsiN_analarray=dpsiNdr_analytic(q,dPsiTdr,psi,r)
    plt.plot(r,dpsiN_array)
    plt.plot(r,dpsiN_analarray)
    #plt.plot(r,1/dpsiN_array)
    plt.title("dpsi_N/dr")
    plt.show()
    r=r[N/2:-1]
    oneOverdpsiN_array=1/dpsiN_analarray[N/2:-1]
    plt.plot(r,oneOverdpsiN_array)
    plt.title("dr/dpsi_N")
    plt.show()
    print drdPsiN_at_psiN_one(q,dPsiTdr,psi,r[-1])
