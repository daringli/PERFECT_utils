from profile import Profile
from profileSplineFromLinear import ProfileSplineFromLinear
from scipy.interpolate import UnivariateSpline
from bezierTransition import derivativeBezierTransition
from mtanh import generate_m2tanh_profile
from setPedestalParams import setPedestalParams,setPedestalStartStop


class Profile3Linear(ProfileSplineFromLinear):
    """Profile generator that takes 3 linear splines and implements a smooth transition between them."""
    def __init__(self,XStart,XStop,ddx_YCore,ddx_YSOL,width=None,XPedStart=None,XPedStop=None,YPed=None,ddx_YPed=None,YLCFS=None,transWidth=None,transWidthToPedWidth=0.2,kind="bezier",profile='',species=''):
        # calculate
        (self.XPedStart,self.XPedStop,self.width) = setPedestalStartStop(XPedStart,XPedStop,width)
        (self.YPed,self.ddx_YPed, self.width, self.YLCFS) = setPedestalParams(YPed,ddx_YPed, width,YLCFS)

        self.kind = kind
        if transWidth is None:
            self.transWidth= width * transWidthToPedWidth
        else:
            self.transWidth = transWidth
        self.XStart = XStart
        self.XStop = XStop
        self.ddx_YCore = ddx_YCore 
        self.ddx_YSOL = ddx_YSOL
        self.transXs = [float("-inf"),self.XPedStart,self.XPedStop,float("inf")]
        self.ddx_Ys = [self.ddx_YCore,self.ddx_YPed,self.ddx_YSOL]
        super(Profile3Linear,self).__init__(self.transXs,self.ddx_Ys,self.XPedStart,self.YPed,profile=profile,species=species)

    def generate(self,profileList=[]):
        if self.kind is "bezier":
            self.f,self.ddx_f = derivativeBezierTransition(self.flist,self.ddx_flist,self.x[1:-1],[[self.transWidth,self.transWidth],[self.transWidth,self.transWidth]])
            print "y skeleton: " + str(self.y)
            print "x skeleton: " + str(self.x)
        elif self.kind is "mtanh":
            self.f,self.ddx_f = generate_m2tanh_profile(self.YPed,self.ddx_YCore,self.ddx_YPed, self.ddx_YSOL, self.width, self.XPedStart)
        p = Profile(self.f)
        ddx_p = Profile(self.ddx_f)
        p.species = self.species
        p.profile = self.profile
        p.generator = '3Linear'
        ddx_p.species = self.species
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.generator = '3Linear'
        return (p,ddx_p)
        
        

            
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    
    # pedestal profile
    width=0.03
    YPed=0.9
    ddx_YCore=-0.1*YPed/width
    ddx_YPed=-YPed/width
    #print dTpeddx_e
    ddx_YSOL=-0.05*YPed/width

    XPedStop = 1

    XPedStart = XPedStop - width
    XStart = XPedStart - width
    XStop =  XPedStart + 2 * width
    x=numpy.linspace(XStart,XStop,50)
    # bezier 1
    gen1 = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_YCore,ddx_YSOL=ddx_YSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=YPed,ddx_YPed=ddx_YPed,transWidthToPedWidth=0.2)
    (f1, ddx_f1) = gen1.generate()
    
    # bezier 2
    gen2 = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_YCore,ddx_YSOL=ddx_YSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=YPed,ddx_YPed=ddx_YPed,transWidthToPedWidth=0.4)
    (f2, ddx_f2) = gen2.generate()
    
    # mtanh
    gen3 = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_YCore,ddx_YSOL=ddx_YSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=YPed,ddx_YPed=ddx_YPed,kind="mtanh")
    (f3, ddx_f3) = gen3.generate()

    # plot
    plt.subplot(2,1,1)
    plt.plot(x,f1(x))
    plt.plot(x,f2(x))
    plt.plot(x,f3(x))

    plt.subplot(2,1,2)
    plt.plot(x,ddx_f1(x))
    plt.plot(x,ddx_f2(x))
    plt.plot(x,ddx_f3(x))
    
    plt.show()
