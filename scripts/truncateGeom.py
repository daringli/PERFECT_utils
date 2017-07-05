#!/usr/bin/python

import h5py
import sys


def find_indices(x,xmin,xmax):
    iMin = 0
    while x[iMin] < xmin:
        iMin = iMin + 1
    iMin = iMin - 1 #bottom is in range

    iMax = 0
    while x[iMax] < xmax:
        iMax = iMax + 1
    iMax = iMax # top is out of range

    return (iMin,iMax)


def truncateGeometry(filename,psiMin,psiMax,new_filename):
    f=h5py.File(filename,'r')
    groupname = f["/"].keys()[0]
    datanames = f["/" + groupname].keys()

    psi = f["/" + groupname + "/psiArray"][()]
    
    iMin,iMax = find_indices(psi,psiMin,psiMax)

    new_psiMin = psi[iMin]
    new_psiMax = psi[iMax-1] # since iMax will not be included


    g=h5py.File(new_filename,'w')
    new_Npsi = iMax - iMin # no +1 since top is out of range
    new_groupname = "Npsi" + str(new_Npsi) + "Ntheta" + groupname.split("Ntheta")[1] 
    new_group = g.create_group(new_groupname)
    
    for dataname in datanames:
        if dataname not in ["EFITFileName", "Npsi", "R0", "Z0", "psi0", "psiMax", "psiMin"]:
            data = f["/" + groupname + "/" + dataname][()][iMin:iMax] # no +1 since top is out of range
        elif dataname in ["EFITFileName", "R0", "Z0", "psi0"]:
            data = f["/" + groupname + "/" + dataname]
        elif dataname == "Npsi":
            data = new_Npsi
        elif dataname == "psiMin":
            data = new_psiMin
        elif dataname == "psiMax":
            data = new_psiMax
        else:
            raise ValueError("Dataset " + dataname + " not supposed to be found in PERFECT geometry file.")
        
        new_group.create_dataset(dataname,data=data)

if __name__=="__main__":
    argv = sys.argv
    filename = argv[1]
    psiMin = float(argv[2])
    psiMax = float(argv[3])
    new_filename = argv[4]
    
    truncateGeometry(filename,psiMin,psiMax,new_filename)
    
