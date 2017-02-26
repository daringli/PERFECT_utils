from profile import Profile
from profileSplineFromLinear import ProfileSplineFromLinear
from scipy.interpolate import UnivariateSpline
from bezier_transition import derivative_bezier_transition
from mtanh import generate_m2tanh_profile
from const_delta import constant_delta_T,constant_delta_dTdpsiN
from numpy import sqrt

class ProfileConstDeltaT(object):
    """Profile generator that takes a deltaT, etc and gives a const-delta T"""
    def __init__(self,T0,X0,Z,deltaT,psiAHat,mHat,RHat,Delta,species=''):
        self.species = species
        self.profile = "T"
        self.psiAHat = psiAHat
        self.Delta = Delta
        self.deltaT = deltaT
        self.Z = Z
        self.mHat = mHat
        self.RHat = RHat
        self.T0 = T0
        self.X0 = X0
        self.A = self.deltaT * self.psiAHat/(sqrt(self.mHat)*self.RHat*self.Delta)
        self.B = self.A/(2*sqrt(self.T0))
        
    def generate(self,profileList=[]):
        self.f = constant_delta_T(self.T0,self.X0,self.Z,self.B)
        self.ddx_f = constant_delta_dTdpsiN(self.T0,self.X0,self.Z,self.B)
        p = Profile(self.f)
        ddx_p = Profile(self.ddx_f)
        p.species = self.species
        p.profile = self.profile
        p.generator = 'constDeltaT'
        ddx_p.species = self.species
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.generator = 'constDeltaT'
        
        return (p,ddx_p)
        
        

            
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    
    # pedestal profile
    deltaT = 0.1
    Z = 1
    RHat =1.3
    mHat = 1.0
    Delta = 0.0006  #based on He papar
    psiAHat = 0.0117  # based on He paper
    T0 = 1.0
    X0 = 0.9
    x=numpy.linspace(0.8,1.099,100)
    # bezier 1
    gen1 = ProfileConstDeltaT(T0,X0,Z,deltaT,psiAHat,mHat,RHat,Delta,species='')
    (f1, ddx_f1) = gen1.generate()
    
    plt.subplot(2,1,1)
    plt.plot(X0,T0,'o')
    plt.plot(x,f1(x))
    plt.subplot(2,1,2)
    plt.plot(x,ddx_f1(x))
    plt.show()
