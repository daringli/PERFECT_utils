import scipy.constants
import sys
import numpy

def nu_r(RBar,nBar,TBar,logLambda):
    e=scipy.constants.e
    ep0=scipy.constants.epsilon_0
    pi=numpy.pi
    #2016-01-12: had a factor sqrt(2) too much due to vBar containing sqrt(2)
    # thus divided old result by sqrt(2).
    prefactor=(1/(12*pi**(3.0/2.0))) * (e**4/(ep0**2))
    return prefactor*logLambda*(RBar*nBar/TBar**2)

if __name__=='__main__':
    from coulomb_logarithms import logLambda_ee
    RBar=float(sys.argv[1])
    nBar=float(sys.argv[2])
    TBar=float(sys.argv[3])
    logLambda=logLambda_ee(nBar,TBar)
    print logLambda
    print nu_r(RBar,nBar,TBar,logLambda)
