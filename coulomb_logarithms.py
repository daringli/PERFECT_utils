import numpy
import sys
import scipy.constants


def logLambda_ee(ne,Te):
    """electron-electron logarithm from NRL page 34 2013
    ne in m^(-3)
    Te in J
    Note: these are coverted into cm^{-3} and eV, and then inserted into the formula:
    logLambda_ee = 23.5-log(ne^(1.0/2.0)*Te^(5.0/4.0))-(10^(-5)+(log(Te)-2)^2/16.0)^(1.0/2.0)

"""
    #converts to cm^(-3),eV from m^(-3),J
    e=scipy.constants.e
    Te=Te/e
    ne=ne/1e6
    return 23.5-numpy.log(ne**(1.0/2.0)*Te**(-5.0/4.0))-(10**(-5)+((numpy.log(Te)-2)**2)/16.0)**(1.0/2.0)

def logLambda_ie(ne,Te):
    """ion-electron logarithm from NRL page 34 2013
    ne in m^(-3)
    Te in J
    Note: these are coverted into cm^{-3} and eV, and then inserted into the formula:
    logLambda_ie = 24.0-log(ne^(1.0/2.0)/Te)
"""
    #converts to cm^(-3),eV from m^(-3),J
    e=scipy.constants.e
    Te=Te/e
    ne=ne/1e6
    return 24.0-numpy.log(ne**0.5/Te)


if __name__=='__main__':
    test_arrays=False
    
    n=float(sys.argv[1])
    T=float(sys.argv[2])

    if test_arrays:
        n=numpy.linspace(0.2e14,0.4e14)
        T=numpy.linspace(450/2.0,450.0)
    print "logLambda_ee:"
    print logLambda_ee(n,T)

    print "logLambda_ie:"
    print logLambda_ie(n,T)
