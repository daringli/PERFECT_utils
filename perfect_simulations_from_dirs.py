import os,sys
from perfect_simulation import perfect_simulation,normalized_perfect_simulation


sentinel=object()

def perfect_simulations_from_dirs(dirlist,normlist,species,psiN_to_psi=sentinel,global_term_multiplier=sentinel,pedestal_points_list=[]):
    #for now, assumes that the simulations have been previously finished
    if type(normlist) is not list:
        normlist=[normlist]*len(dirlist)
    if type(species) is not list:
        species=[species]*len(dirlist)
    if psiN_to_psi != sentinel:
        if type(psiN_to_psi) is not list:
            psiN_to_psi=[psiN_to_psi]*len(dirlist)
    if  global_term_multiplier != sentinel:
        if type(global_term_multiplier) is not list:
            global_term_multiplier=[global_term_multiplier]*len(dirlist)
            

    if pedestal_points_list != []:
        if type(pedestal_points_list) is not list:
            print "perfect_simulations_from_dir: ERROR: pedestal points should be a list or list of lists"
            sys.exit(1)
        elif type(pedestal_points_list[0]) is not list:
            # only simple list, make list of lists
            pedestal_points_list = [pedestal_points_list]*len(dirlist)
    else:
        pedestal_points_list=[[]]*len(dirlist)
    input="input.namelist"
    output="perfectOutput.h5"
    
    i=0
    simuls=[]
    for dir in dirlist:
        input_filename=dir+"/"+input
        output_filename=dir+"/"+output
        norm_filename=normlist[i]
        species_filename=species[i]
        pedestal_points = pedestal_points_list[i]
        if psiN_to_psi != sentinel:
            psiN_to_psi_filename = psiN_to_psi[i]
        else:
            psiN_to_psi_filename = None
        if global_term_multiplier != sentinel:
            global_term_multiplier_filename = global_term_multiplier[i]
        else:
            global_term_multiplier_filename = None
        
        print species_filename
        print psiN_to_psi_filename
        print global_term_multiplier_filename
        
        simuls.append(normalized_perfect_simulation(input_filename,norm_filename,species_filename,output_filename,psiN_to_psi_filename,global_term_multiplier_filename,pedestal_points=pedestal_points))
        simuls[i].description=dir
        #simuls[i].species=species
        i=i+1
    return simuls
        
if __name__=='__main__':
    dirlist=sys.argv[1:]
    normlist="norms.namelist"
    species="species"
    simuls=perfect_simulations_from_dirs(dirlist,normlist,species)
