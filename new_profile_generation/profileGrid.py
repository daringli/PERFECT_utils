from profile import Profile
from profileSplineFromLinear import ProfileSplineFromLinear
from scipy.interpolate import UnivariateSpline
from bezier_transition import derivative_bezier_transition
from mtanh import generate_m2tanh_profile
from setPedestalParams import setPedestalParams,setPedestalStartStop

from nonuniformGrid import generate_nonuniform_grid_Int_arctan


class ProfileGrid(object):
    """Generates a non-uniform grid"""
    def __init__(self,s,psiAHat,XPedStart=None,XPedStop=None,width=None,gridType = "uniform",**kwargs):
        self.profile="grid"
        self.species="none"
        self.Npsi=len(s)
        
        # calculate
        (XPedStart,XPedStop,width) = setPedestalStartStop(XPedStart,XPedStop,width)
        if gridType == "uniform":
            self.psi = lambda x: x*psiAHat
            self.dds_psi = lambda x:  psiAHat + 0*x #trivial mapping 
        elif gridType == "IntArctan":
            a = kwargs["transitionLength"]
            b = kwargs["pedestalGridDensity"]
            c=1

            # by design of non-uniform grid, s[0]=psiN[0], s[-1]=psiN[-1]
            # for the physical psiN. Not true for indices between
            if "leftShift" in kwargs.keys():
                self.leftshift = (XPedStop-XPedStart)*kwargs["leftShift"]
            else:
                self.leftshift = 0

            if "rightShift" in kwargs.keys():
                self.rightshift = (XPedStop-XPedStart)*kwargs["rightShift"]
            else:
                self.rightshift = 0

            X1= XPedStart - self.leftshift
            X2= XPedStop + self.rightshift

            self.psi,self.dds_psi = generate_nonuniform_grid_Int_arctan(s,psiAHat,a,b,c,X1,X2)
        else:
            "generate_compatible_profiles: ERROR: unrecognized grid type!"
            raise ValueError('Unrecognized non-uniform grid option')

        
    
    def generate(self,profileList=[]):
        p = Profile(self.psi)
        ddx_p = Profile(self.dds_psi)
        p.species = self.species
        p.profile = self.profile
        p.generator = 'grid'
        ddx_p.species = self.species
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.generator = 'grid'

        return (p,ddx_p)
        
        

            
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt


    XStart = 0.9
    XPedStart = 0.95
    XPedStop = 1.0
    XStop = 1.05
    x = numpy.linspace(XStart,XStop,50)

    psiAHat = 0.1
    s = numpy.linspace(XStart,XStop,50)
    
    #uniform
    gen1 = ProfileGrid(s,psiAHat,XPedStart,XPedStop,gridType = "uniform")
    psi1,psiHat1 = gen1.generate()
    
    #nonuniform
    kw={"gridType":"IntArctan","transitionLength":0.025,"pedestalGridDensity":0.7,"leftShift":0.025*20,"rightShift":0.025*20}
    gen2 = ProfileGrid(s,psiAHat,XPedStart,XPedStop,**kw)
    psi2,psiHat2 = gen2.generate()
    
    plt.subplot(2,1,1)
    plt.plot(x,psiHat1(x))

    plt.subplot(2,1,2)
    plt.plot(x,psiHat2(x))
    
    
    plt.show()
