from profile import Profile
from profileSplineFromLinear import ProfileSplineFromLinear
from scipy.interpolate import UnivariateSpline
from bezier_transition import derivative_bezier_transition
from mtanh import generate_m2tanh_profile
from const_delta import constant_delta_X,constant_delta_dXdpsiN
from numpy import sqrt

from findInProfileList import findInProfileList
from profileConstDeltaT import ProfileConstDeltaT

class ProfileConstDeltaY(object):
    """Profile generator that takes a deltaY, etc and gives a const-delta Y"""
    def __init__(self, Y0, X0, Z, deltaY, BT, deltaT, species='', profile=''):
        self.Y0 = Y0
        self.X0 = X0
        self.Z = Z
        self.deltaY = deltaY
        self.BT = BT
        self.deltaT = deltaT
        self.species = species
        self.profile = profile
       
    def generate(self,profileList=[]):
        # check if T is a constant delta T
        TIndex = findInProfileList(profileList,"T",self.species,'constDeltaT')
                    
        self.f = constant_delta_X(self.Y0,self.X0,self.Z,self.deltaY,self.BT,self.deltaT)
        self.ddx_f = constant_delta_dXdpsiN(self.Y0,self.X0,self.Z,self.deltaY,self.BT,self.deltaT)
        p = Profile(self.f)
        ddx_p = Profile(self.ddx_f)
        p.species = self.species
        p.profile = self.profile
        p.generator = 'constDeltaEta'
        ddx_p.species = self.species
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.generator = 'constDeltaEta'
        
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
    genT = ProfileConstDeltaT(T0,X0,Z,deltaT,psiAHat,mHat,RHat,Delta,species='')
    eta0=0.4
    deltaEta = 0.1
    
    (T, ddx_T) = genT.generate()
    gen1 = ProfileConstDeltaY(eta0, X0, Z, deltaEta, genT.B, genT.deltaT, species='', profile='')
    (f1, ddx_f1) = gen1.generate([T])
    
    plt.subplot(2,1,1)
    plt.plot(X0,eta0,'o')
    plt.plot(x,f1(x))
    plt.subplot(2,1,2)
    plt.plot(x,ddx_f1(x))
    plt.show()
