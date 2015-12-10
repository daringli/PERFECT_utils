from coulomb_logarithms import lambda_ee
import scipy.constants
import sys
import numpy

def nu_r(RBar,nBar,TBar,logLambda):
    e=scipy.constants.e
    ep0=scipy.constants.epsilon_0
    pi=numpy.pi
    prefactor=(2**(1.0/2.0)/(12*pi**(3.0/2.0))) * (e**4/(ep0**2))
    return prefactor*logLambda*(RBar*nBar/TBar**2)

if __name__=='__main__':
    scipy.constants.e
    RBar=float(sys.argv[1])
    nBar=float(sys.argv[2])
    TBar=float(sys.argv[3])
    logLambda=lambda_ee(nBar,TBar)
    print logLambda
    print nu_r(RBar,nBar,TBar,logLambda)
