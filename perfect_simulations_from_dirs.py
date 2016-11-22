import os,sys
from perfect_simulation import perfect_simulation,normalized_perfect_simulation

def perfect_simulations_from_dirs(dirlist,normlist,species,psiN_to_psi=None,global_term_multiplier=None,pedestal_start_stop_list=None,pedestal_point_list=None, core_point_list = None):
    #for now, assumes that the simulations have been previously finished
    if type(normlist) is not list:
        normlist=[normlist]*len(dirlist)
    if type(species) is not list:
        species=[species]*len(dirlist)
    if psiN_to_psi != None:
        if type(psiN_to_psi) is not list:
            psiN_to_psi=[psiN_to_psi]*len(dirlist)
    if  global_term_multiplier != None:
        if type(global_term_multiplier) is not list:
            global_term_multiplier=[global_term_multiplier]*len(dirlist)
            

    if pedestal_start_stop_list != None:
        if type(pedestal_start_stop_list) is not list:
            pedestal_start_stop_list = [pedestal_start_stop_list]*len(dirlist)
    else:
        pedestal_start_stop_list=[(None,None)]*len(dirlist)
        
    if  pedestal_point_list != None:
        if type(pedestal_point_list) is not list:
            pedestal_point_list=[pedestal_point_list]*len(dirlist)
    else:
        pedestal_point_list=[None]*len(dirlist)

    if  core_point_list != None:
        if type(core_point_list) is not list:
            core_point_list=[core_point_list]*len(dirlist)
    else:
        core_point_list=[None]*len(dirlist)

    input="input.namelist"
    output="perfectOutput.h5"
    
    i=0
    simuls=[]
    for dir in dirlist:
        input_filename=dir+"/"+input
        output_filename=dir+"/"+output
        norm_filename=normlist[i]
        species_filename=species[i]
        pedestal_start_stop = pedestal_start_stop_list[i]
        pedestal_point = pedestal_point_list[i]
        core_point = core_point_list[i]
        if psiN_to_psi != None:
            psiN_to_psi_filename = psiN_to_psi[i]
        else:
            psiN_to_psi_filename = None
        if global_term_multiplier != None:
            global_term_multiplier_filename = global_term_multiplier[i]
        else:
            global_term_multiplier_filename = None
        
        print species_filename
        print psiN_to_psi_filename
        print global_term_multiplier_filename
        
        simuls.append(normalized_perfect_simulation(input_filename,norm_filename,species_filename,output_filename,psiN_to_psi_filename,global_term_multiplier_filename,pedestal_start_stop = pedestal_start_stop,pedestal_point = pedestal_point,core_point = core_point))
        simuls[i].description=dir
        #simuls[i].species=species
        i=i+1
    return simuls
        
if __name__=='__main__':
    dirlist=sys.argv[1:]
    normlist="norms.namelist"
    species="species"
    simuls=perfect_simulations_from_dirs(dirlist,normlist,species)
