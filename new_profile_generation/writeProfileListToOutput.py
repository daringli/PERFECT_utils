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
    s = simul.inputs.psi
    psiAHat = simul.psiAHat
    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,ddx_grid(s),grid(s),grid(s)/psiAHat)
    
    
    # with the grid in place, we sample the profiles on it
    n = [None]*Nspecies
    ddpsiN_n = [None]*Nspecies
    eta = [None]*Nspecies
    ddpsiN_eta = [None]*Nspecies
    T = [None]*Nspecies
    ddpsiN_T = [None]*Nspecies

    for i,specy in enumerate(simul.species_list):
        n[i] = extractProfilesFromList(profileList,"n",specy,"any")[0]
        ddpsiN_n[i] = extractProfilesFromList(profileList,"ddx_n",specy,"any")[0]
        eta[i] = extractProfilesFromList(profileList,"eta",specy,"any")[0]
        ddpsiN_eta[i] = extractProfilesFromList(profileList,"ddx_eta",specy,"any")[0]
        T[i] = extractProfilesFromList(profileList,"T",specy,"any")[0]
        ddpsiN_T[i] = extractProfilesFromList(profileList,"ddx_T",specy,"any")[0]

    Phi = extractProfilesFromList(profileList,"Phi","none","any")[0]
    ddpsiN_Phi = extractProfilesFromList(profileList,"ddx_Phi","none","any")[0]
    psiN = grid(s)/psiAHat

    scale=(ddx_grid/psiAHat)(s)
    Phi = sampleFunclist(Phi,psiN)
    dds_Phi = sampleFunclist(ddpsiN_Phi,psiN)*scale
    
    scale = scale[:,numpy.newaxis]
    T = numpy.transpose(sampleFunclist(T,psiN))
    dds_T = numpy.transpose(sampleFunclist(ddpsiN_T,psiN))*scale
    n = numpy.transpose(sampleFunclist(n,psiN))
    dds_n = numpy.transpose(sampleFunclist(ddpsiN_n,psiN))*scale
    eta = numpy.transpose(sampleFunclist(eta,psiN))
    dds_eta = numpy.transpose(sampleFunclist(ddpsiN_eta,psiN))*scale
    

    # write profiles to file
    pP = perfectProfiles(outputname)
    pP.create_profiles_for_Npsi(Npsi,Nspecies, Phi, dds_Phi, T, dds_T, n, dds_n, eta, dds_eta)

    #write profile parameters to profile specific output file for documentation
    for profile in profileList:
        profile.print_to_file()
    
