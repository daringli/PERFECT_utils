import h5py
import sys

def input_psiN_params_from_geometry_file(filename):
    f=h5py.File(filename,'r')
    groupname = f["/"].keys()[0]
    psiMin  = f["/" + groupname + "/psiMin"][()][0]
    psiMax  = f["/" + groupname + "/psiMax"][()][0]
    psiA = f["/" + groupname + "/psi0"][()][0]
    psiMid = (psiMin + psiMax)/2.0
    psiDiameter = psiMax - psiMin
    
    return (psiMid,psiDiameter,psiA)

if __name__=="__main__":
    argv = sys.argv
    if len(argv) > 1:
        filename = argv[1]
    else:
        # use a default filename
        filename = "input.geometry.h5" 
    (psiMid,psiDiameter,psiA) = input_psiN_params_from_geometry_file(filename)
    print "psiMid: " + "{0:0.16f}".format(psiMid)
    print "psiDiameter: " + "{0:0.16f}".format(psiDiameter)
    print "psiA: " + "{0:0.16f}".format(psiA)
    
