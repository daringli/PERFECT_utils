from profile import Profile
from scipy.interpolate import UnivariateSpline


class ProfileSplineFromLinear(object):
    """Profile generator that takes linear splines and implements a smooth transition between them.
    transXs: list of X at which unsmoothened spline transitions
    ddx_Ys: list of slope in each interval (transXs[i],transXs[i+1]), should have length len(transXs)-1
    X0,Y0: value of spline at one point X0"""
    
    
    def __init__(self,transXs,ddx_Ys,X0,Y0,kind=3,w=None,profile='',species=''):
        self.profile=profile
        self.species=species
        Nintervals = len(transXs) - 1 
        if len(ddx_Ys) != Nintervals:
            print "ProfileSplineFromLinear : error: number of specified gradients do not match the number of intervals!"
            raise ValueError('Incorrect number of intervals!')
        self.ddx_Ys = []
        self.intervals = []
        
        startInterval = None
        for i in range(len(transXs[1:])):
            self.intervals.append([transXs[i],transXs[i+1]])
            self.ddx_Ys.append(ddx_Ys[i])
            if  transXs[i] > transXs[i+1]:
                print "ProfileSplineFromLinear : error: transXs not increasing."
                raise ValueError('transXs should be increasing')
            if X0 >= transXs[i] and X0 < transXs[i+1]:
                startInterval = i
        if startInterval == None:
            print "ProfileSplineFromLinear : error: X0 was not found in any interval."
            raise ValueError('X0 not found in intervals')

        # create list of linear functions that will be the basis of the smooth
        # spline
        self.flist = [None]*Nintervals
        self.ddx_flist = [None]*Nintervals
        
        def setFunctionOnIntervals(i,X0,Y0,increment):
            if i < 0 or i >= Nintervals:
                print "ProfileSplineFromLinear : error: Invalid start interval"
                raise ValueError('Invalid start interval')
            self.flist[i] = lambda x,i=i: Y0  + self.ddx_Ys[i]*(x-X0)
            self.ddx_flist[i] = lambda x,i=i: self.ddx_Ys[i] + 0 *x
            
            if increment == -1:
                X1 = self.intervals[i][0]
                Y1 = self.flist[i](X1)
                if i>0:
                    setFunctionOnIntervals(i+increment,X1,Y1,increment)
            if increment == 1:
                X1 = self.intervals[i][1]
                Y1 = self.flist[i](X1)
                if i<(Nintervals-1):
                    setFunctionOnIntervals(i+increment,X1,Y1,increment)
        setFunctionOnIntervals(startInterval,X0,Y0,1)
        setFunctionOnIntervals(startInterval,X0,Y0,-1)
        
        self.y = [None]*len(transXs)
        for i,transX in enumerate(transXs):
            if i<Nintervals:
                self.y[i] = self.flist[i](transX)
            else:
                self.y[i] = self.flist[i-1](transX)
        self.x = transXs
        
        if type(kind) is int:
            self.type = "smoothingSpline"
            self.order = kind
            self.smoothing = 0
            self.w=w

        elif type(kind) is tuple:
            self.type = "smoothingSpline"
            self.order = kind[0]
            self.smoothing = kind[1]
            self.w=w

    def generate(self,profileList=[]):
        if self.type is "smoothingSpline":
            f = UnivariateSpline(self.x, self.y, k=self.order,s=self.smoothing,w=self.w)
            ddx_f = f.derivative()
        else: 
            print "ProfileSplineFromLinear : error: unsupported interpolation scheme!"
            raise ValueError('Unsupported interpolation scheme!')
        p = Profile(f)
        ddx_p = Profile(ddx_f)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'splineFromLinear'
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'splineFromLinear'
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    
    # pedestal profile
    xwidth=0.03
    Tped_e=0.9
    dTCoredx_e=-0.1*Tped_e/xwidth
    dTpeddx_e=-Tped_e/xwidth
    #print dTpeddx_e
    dTSOLdx_e=-0.05*Tped_e/xwidth

    xPed = 0.95
    xtrans=[xPed-2*xwidth,xPed-xwidth,xPed,xPed + xwidth,xPed + 2*xwidth,xPed + 3*xwidth]
    grads = [dTCoredx_e,dTCoredx_e,dTpeddx_e,dTSOLdx_e,dTSOLdx_e]

    x=numpy.linspace(xtrans[0],xtrans[-1],50)
    # linear spline
    gen1 = ProfileSplineFromLinear(xtrans,grads,xPed,Tped_e,kind=1)
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    #plt.plot(x,ddx_f1(x))
    
    # smoothened spline
    gen2 =  ProfileSplineFromLinear(xtrans,grads,xPed,Tped_e,kind=(1,None))
    (f2,ddx_f2) = gen2.generate()
    plt.plot(x,f2(x))
    #plt.plot(x,ddx_f2(x))
    # crazy interpolating spline
    gen3 =  ProfileSplineFromLinear(xtrans,grads,xPed,Tped_e,profile="test",species="dog",w=[1,1,0.0,0.0,1,1],kind=(2,0.001))
    (f3,ddx_f3) = gen3.generate()
    plt.plot(x,f3(x))
    #plt.plot(x,ddx_f3(x))
    print f3.profile
    print f3.species
    plt.show()
