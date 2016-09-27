from __future__ import division
import numpy 
pi = numpy.pi
arctan=numpy.arctan

from perfectGlobalMultiplierProfileFile import create_globalMultiplier_of_Npsi
#for main
import matplotlib.pyplot as plt



#(1-c)+c\tilde{\chi}_{[a+\delta_a,b-\delta_b]}^{(\Delta)}.
def generate_smooth_heaviside(Delta):
    """Returns smooth approximation of Heaviside step function"""
    return lambda x: 0.5*(1+(2.0/pi)*arctan(x/Delta))
    

def generate_smooth_characteristic(Delta,a,b):
    """Returns smooth approximation of characteristic function of [a,b]"""
    h_Delta=generate_smooth_heaviside(Delta)
    return lambda x: h_Delta(x-a) - h_Delta(x-b)

def generate_smooth_upshifted_characteristic(Delta,a,b,delta_a,delta_b,c):
    """Smooth approximation of: (1-c) + c * chi_[a+delta_a,b-delta_b](x) """
    chi = generate_smooth_characteristic(Delta,a+delta_a,b-delta_b)
    middle = (b-delta_b + a + delta_a)/2.0
    print 
    d2 = (1-c)/(chi(middle)-chi(a))
    d1 = 1-d2*chi(middle)
    return lambda x: d1 + d2*chi(x)

def generate_global_multiplier(a,b,Delta=1/500,delta_a=0.1,delta_b=0.1,c=0.1):
    #psiN: the actual psiN (or psi, normalizations do not matter), as an array.
    #globalTermMultiplierFilename,psiN
    #a=psiN[0]
    #b=psiN[-1]
    width = b - a
    Delta = Delta * width
    delta_a = delta_a * width
    delta_b = delta_b * width
    C=generate_smooth_upshifted_characteristic(Delta,a,b,delta_a,delta_b,c)
    #create_globalMultiplier_of_Npsi(globalTermMultiplierFilename,Npsi,C(psiN))
    return C
    
if __name__=="__main__":
    a=1.
    b=2.
    Delta = (b-a)/500
    delta_a=50*Delta
    delta_b=50*Delta
    c=0.1
    C=generate_smooth_upshifted_characteristic(Delta,a,b,delta_a,delta_b,c)
    x=numpy.linspace(1,2,250)
    plt.plot(x,C(x))
    plt.show()
