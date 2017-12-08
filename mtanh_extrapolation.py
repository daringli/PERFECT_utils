from __future__ import division

from read_IDL_output import read_IDL_output,read_IDL_output_psiN
from mtanh import mtanh_profile
import numpy
from get_index_range import get_index_range
from buffer_extend_uniform_grid import buffer_extend_uniform_grid

# the workhorse of this script is:
from scipy.optimize import curve_fit

verbose = True

def linear_extrapolation(x1,y1,interval=None,endpoint=1.3):
    assert(x1.shape == y1.shape)
    
    #i1 = numpy.argmin(numpy.fabs(x1-0.95))
    if interval is None:
        # linearly extrapolate from last 2 points
        i1 = len(x1) - 2
        i2 = len(x1)
    else:
        [i1,i2] = get_index_range(x1,interval)
        i2 = i2 + 1 #want i2 = len(x) to get last index in x[i1:i2]

    p = numpy.polyfit(x1[i1:i2],y1[i1:i2],1)

    #extend the uniform x-grid to buffer zone
    dx = x1[-1] - x1[-2]
    x2 = buffer_extend_uniform_grid(x1[-1],endpoint,dx,first=False)

    y2 = numpy.polyval(p,x2)
    y = numpy.concatenate((y1[:i2],y2))
    x = numpy.concatenate((x1[:i2],x2))

def mtanh_extrapolation(x1,y1,interval=None,Vp=None,endpoint=1.3):
    """Extrapolates x1 and y1 data by fitting an mtanh function on the interval interval. Vp can be a guess for the value of y1 at the pedestal"""
    
    if interval is None:
        i1 = 0
        i2 = len(x1)
        interval = [x1[i1],x1[i2-1]]
    else:
        [i1,i2] = get_index_range(x1,interval)
        i2 = i2 + 1 #want i2 = len(x) to get last index in x[i1:i2]

    if Vp is None:
        Vp = y1[i1]
    #normalize data to make initial guesses for nonlinear solver standard
    y1 = y1/Vp 

    #fit an mtanh on the interval
    # mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)
    def f(r,a_ped,a_sol,a_etb,a_delta,a_slope):
        return mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)(r)

    # initial guesses
    a_ped0 = 1
    a_sol0 = 1e-10
    a_delta0 = (interval[1] - interval[0])/4 # if interval is the pedestal
    a_etb0 = interval[0]/2 + interval[1]/2
    a_slope0 = -(y1[i1+1] - y1[i1])/(x1[i1+1] - x1[i1]) #assymptotic core slope

    initial_guess = (a_ped0,a_sol0,a_etb0,a_delta0,a_slope0)
    bounds = ([0,0,interval[0],0,0],[numpy.inf,numpy.inf,interval[1],a_etb0+1e-10,numpy.inf]) # want positive values only

    if verbose:
        print "initial parameter guesses:"
        print (a_ped0,a_sol0,a_etb0,a_delta0,a_slope0)
    popt, pcov = curve_fit(f,x1[i1:i2], y1[i1:i2],initial_guess,bounds=bounds)
    [a_ped,a_sol,a_etb,a_delta,a_slope] = popt
    mtanh = mtanh_profile(a_ped,a_sol,a_etb,a_delta,a_slope)

    if verbose:
        print "final parameter values:"
        print (a_ped,a_sol,a_etb,a_delta,a_slope)

    #extend the uniform x-grid to buffer zone
    dx = x1[-1] - x1[-2]
    x2 = buffer_extend_uniform_grid(x1[-1],endpoint,dx,first=False)

    x = numpy.concatenate((x1[:i2],x2))
    y = Vp*mtanh(x) #denormalize y
    return (x,y)

if __name__=="__main__":
    #to visualize
    import matplotlib.pyplot as plt 
    #to get real data
    from read_IDL_output import read_IDL_output,read_IDL_output_psiN 

    
    
    x1 = read_IDL_output_psiN("Ti.dat")
    y1 = read_IDL_output("Ti.dat")
    print x1[-1]
    
    # interval must include pedestal midpoint
    # or the constrained optimization will not find right values
    interval = [0.8,1,0]

    x,y = mtanh_extrapolation(x1,y1,interval)
    plt.plot(x1,y1)
    plt.plot(x,y)

    plt.show()
