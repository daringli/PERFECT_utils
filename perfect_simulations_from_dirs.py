import os,sys
from perfect_simulation import perfect_simulation,normalized_perfect_simulation


sentinel=object()

def perfect_simulations_from_dirs(dirlist,normlist,species,psiN_to_psi=sentinel):
    #for now, assumes that the simulations have been previously finished
    if type(normlist) is not list:
        normlist=[normlist]*len(dirlist)
    if type(species) is not list:
        species=[species]*len(dirlist)
    if psiN_to_psi != sentinel:
        if type(psiN_to_psi) is not list:
            psiN_to_psi=[psiN_to_psi]*len(dirlist)
            

    input="input.namelist"
    output="perfectOutput.h5"
    
    i=0
    simuls=[]
    for dir in dirlist:
        input_filename=dir+"/"+input
        output_filename=dir+"/"+output
        norm_filename=normlist[i]
        species_filename=species[i]
        if psiN_to_psi != sentinel:
            psiN_to_psi_filename = psiN_to_psi[i]
        else:
            psiN_to_psi_filename = None
        print psiN_to_psi_filename
        print species_filename
        
        simuls.append(normalized_perfect_simulation(input_filename,norm_filename,species_filename,output_filename,psiN_to_psi_filename))
        simuls[i].description=dir
        #simuls[i].species=species
        i=i+1
    return simuls
        
if __name__=='__main__':
    dirlist=sys.argv[1:]
    normlist="norms.namelist"
    species="species"
    simuls=perfect_simulations_from_dirs(dirlist,normlist,species)
