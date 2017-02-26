from __future__ import division

from perfectProfilesFile import perfectProfiles
from perfectPsiAHatProfileFile import create_psiAHat_of_Npsi

from findInProfileList import extractProfilesFromList
import numpy

def sampleFunclist(f_i,x):
    #samples a list of 1d functions f_i
    #returns a 2D array A_ij = f_i(x_j)
    if isinstance(f_i, (list, tuple, numpy.ndarray)):
        f_ij = [f(x) for f in f_i]
    else:
        #if not a list, we just sample
        f_ij = f_i(x)
    return numpy.array(f_ij)


def writeProfileListToOutput(profileList,simul):
    outputname = simul.input_profiles_filename
    Npsi = simul.Npsi
    Nspecies = simul.Nspecies

    # create the nonuniform grid
    grid = extractProfilesFromList(profileList,"grid","none","any")[0]
    ddx_grid = extractProfilesFromList(profileList,"ddx_grid","none","any")[0]
    psiN = simul.inputs.psi
    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,ddx_grid(psiN),grid(psiN))
    psiAHat = simul.psiAHat
    
    # with the grid in place, we sample the profiles on it
    n = [None]*Nspecies
    dds_n = [None]*Nspecies
    eta = [None]*Nspecies
    dds_eta = [None]*Nspecies
    T = [None]*Nspecies
    dds_T = [None]*Nspecies
    
    for i,s in enumerate(simul.species_list):
        n[i] = extractProfilesFromList(profileList,"n",s,"any")[0]
        dds_n[i] = extractProfilesFromList(profileList,"ddx_n",s,"any")[0]*ddx_grid/psiAHat
        eta[i] = extractProfilesFromList(profileList,"eta",s,"any")[0]
        dds_eta[i] = extractProfilesFromList(profileList,"ddx_eta",s,"any")[0]*ddx_grid/psiAHat
        T[i] = extractProfilesFromList(profileList,"T",s,"any")[0]
        dds_T[i] = extractProfilesFromList(profileList,"ddx_T",s,"any")[0]*ddx_grid/psiAHat

    Phi = extractProfilesFromList(profileList,"Phi","none","any")[0]
    dds_Phi = extractProfilesFromList(profileList,"ddx_Phi","none","any")[0]*ddx_grid/psiAHat

    T = numpy.transpose(sampleFunclist(T,psiN))
    dds_T = numpy.transpose(sampleFunclist(dds_T,psiN))
    n = numpy.transpose(sampleFunclist(n,psiN))
    dds_n = numpy.transpose(sampleFunclist(dds_n,psiN))
    eta = numpy.transpose(sampleFunclist(eta,psiN))
    dds_eta = numpy.transpose(sampleFunclist(dds_eta,psiN))
    Phi = sampleFunclist(Phi,psiN)
    dds_Phi = sampleFunclist(dds_Phi,psiN)

    # write profiles to file
    pP = perfectProfiles(outputname)
    pP.create_profiles_for_Npsi(Npsi,Nspecies, Phi, dds_Phi, T, dds_T, n, dds_n, eta, dds_eta)
    
