from __future__ import division

from profile import Profile
from generator import Generator
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from findInProfileList import extractProfilesFromList
from numpy import log,exp
import numpy

import matplotlib.pyplot as plt


class ProfilePhiFromEtaNAdvanced(Generator):
    """Profile generator that takes eta and n profile of some species and gives Phi."""
    
    def __init__(self,eta_Zs,eta_species,n_Zs,n_species,DeltaOver2omega=1.0):
        # solves for Phi from etas and quasineutrality
        self.DeltaOver2omega = DeltaOver2omega
        self.eta_Zs = eta_Zs
        self.eta_species = eta_species
        self.n_Zs = n_Zs
        self.n_species = n_species
        self.profile = "Phi"
        self.species = "none"
        
    def generate(self,profileList=[]):
        #extract needed profiles
        eta_species_wantedProfiles = ["T","ddx_T","eta","ddx_eta"]
        n_species_wantedProfiles = ["n","ddx_n"]
        
        Ts = []
        ddx_Ts = []
        etas = []
        ddx_etas = []
        ns = []
        ddx_ns = []

        profileList = list(profileList)

        for s in self.eta_species:
            [T,ddx_T,eta,ddx_eta] = extractProfilesFromList(profileList,eta_species_wantedProfiles,s)
            Ts.append(T)
            ddx_Ts.append(ddx_T)
            etas.append(eta)
            ddx_etas.append(ddx_eta)

        for s in self.n_species:
            [n,ddx_n] = extractProfilesFromList(profileList,n_species_wantedProfiles, s)
            ns.append(n)
            ddx_ns.append(ddx_n)

        #print ns[0](0.75054)
        #print ns[0].generator_dict["y"][1]

        _x = numpy.linspace(0.8,1.1)
        #plt.plot(_x,ns[0](_x))
        #plt.ylabel("during Phi")
        #plt.show()

        Neta = len(etas)
        Nn = len(ns)
        def X(x):
            # X = exp(-Phi) 
            def f(y):
                return numpy.sum([self.eta_Zs[i]*etas[i](x) * exp(-y*self.eta_Zs[i]/Ts[i](x)) for i in range(Neta)]) + numpy.sum([self.n_Zs[i]*ns[i](x) for i in range(Nn)])
                
            def ddy_f(y):
                ret = numpy.sum([-(self.eta_Zs[i]**2 * etas[i](x)/Ts[i](x)) * exp(-y*self.eta_Zs[i]/Ts[i](x)) for i in range(Neta)])
                return numpy.array([ret]) # wants row of Jacobian for all variables. Does not work the same with a scalar for 1D. Bug?
            
            x0s = [0.0,1.0,-1.0]
            tol = 1e-10
            for x0 in x0s:
                sol = fsolve(f,x0,fprime = ddy_f)[0]
                #sol = fsolve(f,x0)[0]
                if numpy.fabs(f(sol)) < tol:
                    return sol
            
            raise ValueError("PhiFromEtaNAdvanced: fsolve did not converge to tolerance '" + str(tol) + "' for initial guesses '" + str(x0) + "'!")
    
        Phi = numpy.vectorize(X)
        ddx_Phi = numpy.vectorize(lambda x : (1.0/(numpy.sum([self.eta_Zs[i]**2 * (etas[i](x)/Ts[i](x)) * exp(-self.eta_Zs[i]*Phi(x)/Ts[i](x)) for i in range(Neta)]))) * (numpy.sum([self.n_Zs[i]*ddx_ns[i](x) for i in range(Nn)]) + numpy.sum([(self.eta_Zs[i]*ddx_etas[i](x) + self.eta_Zs[i]**2 * (etas[i](x)*ddx_Ts[i](x)/Ts[i](x)**2) * Phi(x)) * exp(-self.eta_Zs[i] * Phi(x)/Ts[i](x)) for i in range(Nn)])))
        
        p = Profile(Phi)
        ddx_p = Profile(ddx_Phi)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'PhiFromEtaN'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'PhiFromEtaN'
        ddx_p.generator_dict = vars(self)
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    from profile3Linear import Profile3Linear

    width=0.03
    XPedStop = 1.0
    XPedStart = XPedStop - width
    xPed = 0.97
    XStart = xPed - width
    XStop =  XPedStop #xPed + 2 * width
    
    x=numpy.linspace(XStart,XStop,50)

    nPed=0.4
    ddx_nCore=-0.1*nPed/width
    ddx_nPed=-nPed/width
    ddx_nSOL=-0.05*nPed/width

    etaPed=0.4
    ddx_etaCore=-0.1*nPed/width
    ddx_etaPed=-0.1*nPed/width
    ddx_etaSOL=-0.05*nPed/width

    TPed=0.9
    ddx_TCore=-0.1*TPed/width
    ddx_TPed=-0.1*TPed/width
    ddx_TSOL=-0.05*TPed/width
    
    # bezier curves
    genn = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_nCore,ddx_YSOL=ddx_nSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=nPed,ddx_YPed=ddx_nPed,transWidthToPedWidth=0.2,profile="n")
    geneta = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_etaCore,ddx_YSOL=ddx_etaSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=etaPed,ddx_YPed=ddx_etaPed,transWidthToPedWidth=0.2,profile="eta")
    genT = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_TCore,ddx_YSOL=ddx_TSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=TPed,ddx_YPed=ddx_TPed,transWidthToPedWidth=0.2,profile="T")
    #genn = Profile3Linear(XStart,XStop,XPedStart,XPedStop,nPed,ddx_nCore,ddx_nPed,ddx_nSOL,width,transWidthToPedWidth=0.2,profile="n")
    #geneta = Profile3Linear(XStart,XStop,XPedStart,XPedStop,etaPed,ddx_etaCore,ddx_etaPed,ddx_etaSOL,width,transWidthToPedWidth=0.2,profile="eta")
    #genT = Profile3Linear(XStart,XStop,XPedStart,XPedStop,TPed,ddx_TCore,ddx_TPed,ddx_TSOL,width,transWidthToPedWidth=0.2,profile="T")
    
    (n, ddx_n) = genn.generate()
    (eta, ddx_eta) = geneta.generate()
    (T, ddx_T) = genT.generate()
    profiles = [n,ddx_n,eta,ddx_eta,T,ddx_T]
    
    # potential
    Z = 1
    Delta = 0.0006  #based on He papar
    omega = Delta/2
    DeltaOver2omega=Delta/(2*omega)
    genPhi = ProfilePhiFromEtaNAdvanced([1.0],[''],[-1.0],[''],DeltaOver2omega)
    (Phi,ddx_Phi) = genPhi.generate(profiles)
    profiles.append(Phi)
    profiles.append(ddx_Phi)

    Z1=1
    Z2=-1
    Nn=1
    Neta=1
    def X2(x):
        # X = exp(-Phi) 
        def f2(y):
            return numpy.sum([Z1*eta(x) * exp(-y*Z1/T(x)) for i in range(Neta)]) + numpy.sum([Z2*n(x) for i in range(Nn)])
        print "f2 " + str(f2(x))
        return fsolve(f2,2.2)[0]

    #print X2(x)
    
    print Phi(0.98)
    print X2(0.98)

    
    
    # plt.plot(x,n(x))
    # plt.plot(x,eta(x))
    
    # plt.show()

    plt.plot(x,Phi(x))
    plt.show()

    
    plt.plot(x,ddx_Phi(x))
    
    num_ddx_Phi = (Phi(x[1:]) - Phi(x[:-1]))/(x[1:] - x[:-1])
    plt.plot(x[1:],num_ddx_Phi)
    
    plt.show()

    
    ni = lambda x : eta(x) * exp(-Phi(x)/T(x))
    plt.plot(x,ni(x))
    plt.plot(x,n(x))
    plt.show()
