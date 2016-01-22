import numpy
import matplotlib.pyplot as plt
import scipy.optimize
from latexify import latexify


def func_q(q0,q95,x,qx):
    #from "Safety factor profile shape and the ideal magnetohydrodynamic stability of n = 1 internal modes"
    def to_solve(nu):
        a=((q95/q0)**nu-1)*0.95**(-2*nu)
        return q0*(1+a*x**(2*nu))**(1/nu)-qx 
    nu=2.0
    nu0=2.0
    nu=scipy.optimize.fsolve(to_solve,2)
    if nu>10:
        print "warning, high nu!"
    def q(rho):
        a=lambda nu: ((q95/q0)**nu-1)*0.95**(-2*nu)
        return q0*(1+a(nu)*rho**(2*nu))**(1/nu)
    return q

def func_q_exp(q0,q95):
    #failed attempt to describe q by e^(-a(x/(x-1)))
    #captures flatness near 0, but too sudden increase by r/a=1
    def q(x):
        a=((100-95)/95.0)*numpy.log(q95/q0)
        print a
        return q0*numpy.exp(-a*(x/(x-1)))
    return q
    

def func_polyfit(x,y,deg):
    p=numpy.polyfit(x,y,deg)
    deg=len(p)-1
    def temp(x):
        out=0
        for i in range(len(p)):
            out=out + p[i]*x**(deg-i)
        return out
    return temp

def func_polyfit_derivative(x,y,deg):
    p=numpy.polyfit(x,y,deg)
    deg=len(p)-1
    def temp(x):
        out=0
        for i in range(len(p)-1):
            out=out + (deg-i)*p[i]*x**(deg-i-1)
        return out
    return temp

def toroidal_flux_derivative(a,kappa,dkappadrho,BT,dBTdr):
    return lambda r :numpy.pi*((1.0/a)*dkappadrho(r/a)*BT(r)*r**2+2*r*kappa(r/a)*BT(r)+dBTdr(r)*kappa(r/a)*r**2)


if __name__=='__main__':
    latexify()
    y=numpy.array([1.23,1.37,1.6])
    x=numpy.array([0,0.5,0.9])
    rho=numpy.linspace(0,0.95)
    deg=2
    kappa=func_polyfit(x,y,deg)
    #dkappadr=func_polyfit_derivative(x,y,deg)
    q=func_q(1.0,3.5,0.6,1.6)
    #q=func_q(4.0,4.0,0.6,4.0)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_ylim([0,4])
    plt.plot(rho,kappa(rho),'--')
    #plt.plot(rho,dkappadr(rho))
    
    ax.set_xlabel(r"$\rho=r/a$",fontsize=12)
    ax.set_ylabel("$\kappa$, $q$",fontsize=12)
    fig.subplots_adjust(bottom=0.20)

    #ax = fig.add_subplot(2,1,2)
    plt.plot(rho,q(rho))
    #ax.set_xlabel("r/a")
    #ax.set_ylabel("q")
    #print q(rho)
    #print q(0.8)
    plt.savefig("q_kappa.pdf")
