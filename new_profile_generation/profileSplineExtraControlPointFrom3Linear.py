from profile import Profile
from profileSplineFromLinear import ProfileSplineFromLinear
from scipy.interpolate import UnivariateSpline

class ProfileFrom3Linear(ProfileSplineFromLinear):
    """Profile generator that takes 3 linear splines and implements a smooth transition between them."""
    def __init__(self,XStart,XStop,XPedStart,XPedStop,YPed,ddx_YCore,ddx_YPed,ddx_YSOL,width,transWidth=None,profile='',species=''):
        XPedStop = XPedStart + width
        if transWidth is None:
            transWidth= width * 0.2
        transXs = [XStart,XPedStart-transWidth/2,XPedStart,XPedStart+transWidth/2,XPedStop-transWidth/2,XPedStop,XPedStop+transWidth/2,XStop]
        ddx_Ys = [ddx_YCore,ddx_YCore,ddx_YPed,ddx_YPed,ddx_YPed,ddx_YSOL,ddx_YSOL]
        w = [1.0,1.0,0.2,1.0,1.0,0.2,1.0,1.0]
        
        
        
        super(ProfileFrom3Linear,self).__init__(transXs,ddx_Ys,XPedStart,YPed,w=w,profile=profile,species=species,kind=(2,0.08))

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

    xPed = 0.97

    XStart = xPed - width
    XStop =  xPed + 2 * width
    XPedStart = xPed
    XPedStop = xPed + width
    x=numpy.linspace(XStart,XStop,50)
    # linear spline
    gen1 = ProfileFrom3Linear(XStart,XStop,XPedStart,XPedStop,YPed,ddx_YCore,ddx_YPed,ddx_YSOL,width)
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    #plt.plot(x,ddx_f1(x))
    
    plt.show()
