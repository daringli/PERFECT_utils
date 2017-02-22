from profile import Profile
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import UnivariateSpline


class ProfileGeneratorFromLinear(object):
    """Profile generator that takes linear splines and implements a smooth transition between them.
    transXs: list of X at which unsmoothened spline transitions
    ddx_Ys: list of slope in each interval (transXs[i],transXs[i+1]), should have length len(transXs)-1
    transXStartStops: list of tuples (x_i1,x_i2) that determine the transition interval at each transX: (transX-x_i1, tranX + x_i2)
    X0,Y0: value of spline at one point X0"""
    
    
    def __init__(self,transXs,ddx_Ys,X0,Y0,profile='',species=''):
        self.profile=profile
        self.species=species
        self.X0 = X0
        self.Y0 = Y0
        
        

    def generate(self,profileList=[]):
        if self.type is "smoothingSpline":
            f = UnivariateSpline(self.x, self.y, k=self.order,s=self.smoothing)
            ddx_f = f.derivative()
        else: 
            print "ProfileGeneratorFromData : error: unsupported interpolation scheme!"
            raise ValueError('Unsupported interpolation scheme!')
        p = Profile(f)
        ddx_p = Profile(ddx_f)
        p.profile = self.profile
        p.species = self.species
        ddx_p.profile = self.profile
        ddx_p.species = self.species
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    
    # data
    xp = numpy.linspace(0,1,8)
    yp=numpy.sin(xp*5.0)
    ddx_yp = 5*numpy.cos(xp*5.0)
    plt.plot(xp,yp,'o')
    plt.plot(xp,ddx_yp,'D')

    x=numpy.linspace(0,1,50)
    # cubic interpolating spline
    gen1 = ProfileGeneratorFromData(xp,yp)
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    plt.plot(x,ddx_f1(x))
    # linear interpolating spline
    gen2 = ProfileGeneratorFromData(xp,yp,1)
    (f2,ddx_f2) = gen2.generate()
    plt.plot(x,f2(x))
    plt.plot(x,ddx_f2(x))
    # crazy interpolating spline
    gen3 = ProfileGeneratorFromData(xp,yp,5,profile="test",species="dog")
    (f3,ddx_f3) = gen3.generate()
    plt.plot(x,f3(x))
    plt.plot(x,ddx_f3(x))
    print f3.profile
    print f3.species
    plt.show()
