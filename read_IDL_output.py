import numpy

def read_IDL_output(filename):
    data = numpy.genfromtxt(filename,comments=';')
    return data[:,1]


def read_IDL_output_psiN(filename):
    data = numpy.genfromtxt(filename,comments=';')
    return data[:,0]


if __name__=="__main__":
    Ti = read_IDL_output("Ti.dat")
    psiN = read_IDL_output_psiN("Ti.dat")

    print Ti
    print psiN
