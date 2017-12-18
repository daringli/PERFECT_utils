#!/usr/bin/env python

from __future__ import division
from perfect_simulation import perfect_simulation
from divide_interval_remainder import divide_interval_remainder
from split_array import split_array
from perfectProfilesFile import perfectProfiles
from perfectGeometryFile import perfectGeometry
from perfectPsiAHatProfileFile import create_psiAHat_of_Npsi

def perfect_splitter(d,N):
    """splits the domain of a perfect simulation in the directory d into N chunks. If Npsi of the original simulation is not divisible by N, the last chunk (at highest psi) will be smaller."""
    
    input_filename = d+"/input.namelist"
    output_filename = d + "/perfectOutput.h5"
    species_filename = d + "/species"
    psiN_to_psi_filename = d  + "/psiAHat.h5"
    simul = perfect_simulation(input_filename,output_filename=output_filename,species_filename=species_filename,psiN_to_psi_filename=psiN_to_psi_filename,global_term_multiplier_filename=None,group_name=None,pedestal_start_stop=(None,None),pedestal_point = None,core_point=None)
    Npsi = simul.Npsi
    Ntheta = simul.Ntheta
    Nspecies = simul.Nspecies

    Npsis = divide_interval_remainder(Npsi,N)
    split_points = [sum(Npsis[0:i]) for i in range(0,N+1)]

    print split_points
    
    # we will split the inputs of our simulation with this array
    splitter = lambda x: split_array(x,split_points)
    
    psiNs = splitter(simul.psi)

    inputs = ["PhiHat","dPhiHatdpsiN","THat","dTHatdpsiN","nHat","dnHatdpsiN","etaHat","detaHatdpsiN","BHat","ddpsiN_BHat","ddtheta_BHat","RHat","JHat","FSA_IHat","ddpsiN_IHat","psiAHatArray","actual_psi","actual_psiN"]

    Ninputs = len(inputs)
    split_inputs = [None] * Ninputs
    for i in range(0,Ninputs):
        split_inputs[i] = splitter(getattr(simul,inputs[i]))
    
    print len(split_inputs)
    print len(split_inputs[i])
    print len(Npsis)

    for i in range(0,N):
        pp = perfectProfiles("p" + str(i) + ".h5")
        pp.create_profiles_for_Npsi(Npsis[i], Nspecies, split_inputs[0][i], split_inputs[1][i], split_inputs[2][i], split_inputs[3][i], split_inputs[4][i], split_inputs[5][i], split_inputs[6][i], split_inputs[7][i])

        pg = perfectGeometry("g" + str(i) + ".h5")
        pg.create_geometry_for_Npsi_Ntheta(Npsis[i], Ntheta, psiNs[i][0], psiNs[i][-1], split_inputs[8][i], split_inputs[9][i], split_inputs[10][i], split_inputs[11][i], split_inputs[12][i], split_inputs[13][i], split_inputs[14][i])

        create_psiAHat_of_Npsi("psiA" + str(i) + ".h5", Npsis[i], split_inputs[15][i],split_inputs[16][i],split_inputs[17][i])
        

if __name__=="__main__":
    from sys import argv
    d = argv[1]
    N = int(argv[2])

    perfect_splitter(d,N)
