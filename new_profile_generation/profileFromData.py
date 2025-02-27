from profile import Profile
from generator import Generator
from scipy.interpolate import UnivariateSpline
import numpy

class ProfileFromData(Generator):
    """Profile generator that takes data y sampled at points x and interpolates/fits using the kind of interpolation specified by kind."""
    
    def __init__(self,x,y,kind=3,profile='',species='',ddx_y=None):
        self.x = x
        self.y = y
        self.ddx_y = ddx_y
        self.profile=profile
        self.species=species
        if type(kind) is int:
            self.type = "smoothingSpline"
            self.order = kind
            self.smoothing = 0
        
    def generate(self,profileList=[]):
        if self.type is "smoothingSpline":
            if numpy.any(self.x[1:] <= self.x[:-1]):
                raise ValueError("x vector to profileFromData must be monotonic for scipy.interpolate.UnivariateSpline In profile: " + self.profile + "_" + self.species + ", offending index: " + str(numpy.where(numpy.diff(self.x)<=0)))
            f = UnivariateSpline(self.x, self.y, k=self.order,s=self.smoothing)
            if self.ddx_y is None:
                ddx_f = f.derivative()
            else:
                ddx_f = UnivariateSpline(self.x, self.ddx_y, k=self.order,s=self.smoothing)
        else: 
            print "ProfileFromData : error: unsupported interpolation scheme!"
            raise ValueError('Unsupported interpolation scheme!')
        p = Profile(f)
        ddx_p = Profile(ddx_f)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'fromData'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'fromData'
        ddx_p.generator_dict = vars(self)
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
    gen1 = ProfileFromData(xp,yp)
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    plt.plot(x,ddx_f1(x))
    # linear interpolating spline
    gen2 = ProfileFromData(xp,yp,1)
    (f2,ddx_f2) = gen2.generate()
    plt.plot(x,f2(x))
    plt.plot(x,ddx_f2(x))
    # crazy interpolating spline
    gen3 = ProfileFromData(xp,yp,5,profile="test",species="dog")
    (f3,ddx_f3) = gen3.generate()
    plt.plot(x,f3(x))
    plt.plot(x,ddx_f3(x))
    print f3.profile
    print f3.species

    # externally set derivative to something (here, xp)
    gen1 = ProfileFromData(xp,yp,ddx_y=xp)
    (f4, ddx_f4) = gen1.generate()
    plt.plot(x,f4(x))
    plt.plot(x,ddx_f4(x))
    
    plt.show()
