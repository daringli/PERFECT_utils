import numpy
import sys
import scipy.constants


def lambda_ee(ne,Te):
    """electron-electron logarithm from NRL page 34 2013
    ne in m^(-3) = 10^(-6) m^(-3)
    Te in keV
    lambda_ee=log(Lambda) = 23.5-log(ne^(1.0/2.0)*Te^(5.0/4.0))-(10^(-5)+(log(Te)-2)^2/16.0)^(1.0/2.0)

"""
    m_keV=True
    #converts to cm^(-3),eV from m^(-3),J
    if m_keV:
        e=scipy.constants.e
        Te=Te/e
        ne=ne/1e6
        #print "ne:"
        #print ne
        #print "Te:"
        #print Te
    return 23.5-numpy.log(ne**(1.0/2.0)*Te**(-5.0/4.0))-(10**(-5)+((numpy.log(Te)-2)**2)/16.0)**(1.0/2.0)

if __name__=='__main__':
    test_arrays=False
    
    n=float(sys.argv[1])
    T=float(sys.argv[2])

    if test_arrays:
        n=numpy.linspace(0.2e14,0.4e14)
        T=numpy.linspace(450/2.0,450.0)
    print "lambda_ee:"
    print lambda_ee(n,T)
