from profile import Profile
from generator import Generator
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import numpy

class ProfileExtendFromData(Generator):
    """Profile generator that takes data y sampled at points x and interpolates/fits using the kind of interpolation specified by kind. 
After interpolation, it extrends the profile outside a range (typically outside the last-closed flux-surface to create a numerical buffer zone).
Extension is designed to be smooth and to let local theory be valid at boundary"""
    
    def __init__(self,x,y,interval,kind=3,profile='',species='',extrapolation=2.0,offset=1.0,a=0.05,b=-0.15,C=7.0):
        self.x = x
        self.y = y
        self.profile=profile
        self.species=species
        self.x0 = interval[1]

        if type(kind) is int:
            self.type = "smoothingSpline"
            self.order = kind
            self.smoothing = 0

        if type(extrapolation) is float:
            self.fraction = extrapolation
            self.extrapolation = "derivativeFraction"
            self.offset=offset

        elif extrapolation == "istvan":
            self.extrapolation = extrapolation
            self.a=a
            self.b=b
            self.C=C
            
        
    def generate(self,profileList=[]):
        if self.type is "smoothingSpline":
            if numpy.any(self.x[1:] <= self.x[:-1]):
                raise ValueError("x vector to profileFromData must be monotonic for scipy.interpolate.UnivariateSpline")
            f = UnivariateSpline(self.x, self.y, k=self.order,s=self.smoothing)
            ddx_f = f.derivative()
        else: 
            print "ProfileFromData : error: unsupported interpolation scheme!"
            raise ValueError('Unsupported interpolation scheme!')
        
        if self.extrapolation == "derivativeFraction":

            ddx_g = lambda x :  ddx_f(x)/(self.fraction*(x - self.x0) + self.offset)
            temp = numpy.vectorize(lambda x: quad(ddx_g,self.x0,x)[0])
            g = lambda x : f(self.x0) + temp(x)
            def h(x):
                if x > self.x0:
                    return g(x)
                else:
                    return f(x)

            def ddx_h(x):
                if x > self.x0:
                    return ddx_g(x)
                else:
                    return ddx_f(x)
            p = Profile(numpy.vectorize(h))
            ddx_p = Profile(numpy.vectorize(ddx_h))
        
            
        elif self.extrapolation == "istvan":
            #h=Log{ Exp[C f(psi)]+Exp[C (a - b psi) ] } /C
            h = lambda x: numpy.log(numpy.exp(self.C * f(x)) + numpy.exp(self.C * (self.a + self.b*x)))/self.C
            ddx_h =  lambda x: (ddx_f(x) * numpy.exp(self.C*f(x)) + self.b * numpy.exp(self.C*(self.a + self.b*x)) )/(numpy.exp(self.C * f(x)) + numpy.exp(self.C * (self.a + self.b*x)))
            p = Profile(h)
            ddx_p = Profile(ddx_h)
        else: 
            print "ProfileFromData : error: unsupported extension scheme!"
            raise ValueError('Unsupported extension scheme!')

        
    
        p.profile = self.profile
        p.species = self.species
        p.generator = 'extendedFromData'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'extendedFromData'
        ddx_p.generator_dict = vars(self)
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    
    # data
    xp = numpy.linspace(0,1,20)
    yp= 3 -2 * xp
    ddx_yp = -2 + 0.0*xp
    
    plt.plot(xp,yp,'o')
    plt.plot(xp,ddx_yp,'s')
    

    gen1 = ProfileExtendFromData(xp,yp,[0,1],kind=3,extrapolation = 2.0,offset=1.0)
    (f1, ddx_f1) = gen1.generate()
    x = numpy.linspace(0,5,100)
    plt.plot(x,f1(x))
    plt.plot(x,ddx_f1(x))
    plt.show()

    gen2 = ProfileExtendFromData(xp,yp,[0,1],kind=3,extrapolation = "istvan")
    (f2, ddx_f2) = gen2.generate()
    x = numpy.linspace(0,5,300)
    plt.plot(x,f2(x))
    plt.plot(x,ddx_f2(x))

    # crude comparison with numerical derivative
    dx = x[1]-x[0]
    ddx_f2_num = (f2(x[1:]) - f2(x[:-1]))/dx
    plt.plot(x[:-1],ddx_f2_num)
    plt.show()
