#!/usr/bin/env python

from __future__ import division
from perfect_simulation import perfect_simulation
from divide_interval_remainder import divide_interval_remainder
from split_array import split_array
from perfectProfilesFile import perfectProfiles
from perfectGeometryFile import perfectGeometry
from perfectPsiAHatProfileFile import create_psiAHat_of_Npsi
from perfectInputFile import perfectInput

from subdir_creator import subdir_creator

import subprocess

def perfect_splitter_distributor(d,N,jobfile=None):
    profile_files, geometry_files,grid_files,input_files = perfect_splitter(d,N)
    subdirs = subdir_creator(d,N)

    for i in range(0,N):
        subprocess.call(["mv",profile_files[i],subdirs[i] +"/input.profiles.h5"])
        subprocess.call(["mv",geometry_files[i],subdirs[i] +"/input.geometry.h5"])
        subprocess.call(["mv",grid_files[i],subdirs[i] +"/psiAHat.h5"])
        subprocess.call(["mv",input_files[i],subdirs[i] +"/input.namelist"])
        subprocess.call(["cp","norms.namelist","species",subdirs[i]])
        if jobfile is not None:
            subprocess.call(["cp",jobfile,subdirs[i]])
    return subdirs

def perfect_splitter(d,N,profile_filename = "profile.input", geometry_filename = "geometry.input", psiA_filename = "psiA", namelist_filename="input", filename_index_offset=0):
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

    profile_filenames = [None]*N
    geometry_filenames = [None]*N
    psiA_filenames = [None]*N
    input_filenames = [None]*N
    for i in range(0,N):
        filename = profile_filename + str(filename_index_offset + i) + ".h5"
        profile_filenames[i] = filename
        pp = perfectProfiles(filename)
        pp.create_profiles_for_Npsi(Npsis[i], Nspecies, split_inputs[0][i], split_inputs[1][i], split_inputs[2][i], split_inputs[3][i], split_inputs[4][i], split_inputs[5][i], split_inputs[6][i], split_inputs[7][i])

        filename = geometry_filename + str(filename_index_offset + i) + ".h5"
        geometry_filenames[i] = filename
        pg = perfectGeometry(filename)
        pg.create_geometry_for_Npsi_Ntheta(Npsis[i], Ntheta, psiNs[i][0], psiNs[i][-1], split_inputs[8][i], split_inputs[9][i], split_inputs[10][i], split_inputs[11][i], split_inputs[12][i], split_inputs[13][i], split_inputs[14][i])

        filename = psiA_filename + str(filename_index_offset + i) + ".h5"
        psiA_filenames[i] = filename
        create_psiAHat_of_Npsi(filename, Npsis[i], split_inputs[15][i],split_inputs[16][i],split_inputs[17][i])

        filename = namelist_filename + str(filename_index_offset + i) + ".namelist"
        input_filenames[i] = filename
        subprocess.call(["cp",d +"/input.namelist", filename])
        #modify input.namelist with new limits
        pI = perfectInput(filename)
        psiMin=psiNs[i][0]
        psiMax=psiNs[i][-1]
        pI.psiDiameter = psiMax - psiMin
        pI.psiMid = (psiMin + psiMax)/2
        pI.widthExtender = 0
        pI.leftBoundaryShift = 0
        pI.rightBoundaryShift = 0
        
    return (profile_filenames,geometry_filenames,psiA_filenames,input_filenames)

if __name__=="__main__":
    from sys import argv
    if len(argv) == 3:
        d = argv[1]
        N = int(argv[2])
        jobfile = None
    elif len(argv) == 4:
        d = argv[1]
        N = int(argv[2])
        jobfile = argv[3]
    else:
        raise ValueError("perfect_splitter supports following format: <d N> or <d N jobname>")

    perfect_splitter_distributor(d,N,jobfile)
