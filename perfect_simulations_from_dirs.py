import os,sys
from perfect_simulation import perfect_simulation,normalized_perfect_simulation


sentinel=object()

def perfect_simulations_from_dirs(dirlist,normlist,species):
    #for now, assumes that the simulations have been previously finished
    if type(normlist) is not list:
        #assume that we only have one file describing the normalization
        #rather than one in each subdirectory
        #normlist=[x+"/"+normlist for x in dirlist]
        normlist=[normlist]*len(dirlist)

    input="input.namelist"
    output="perfectOutput.h5"
    
    i=0
    simuls=[]
    for dir in dirlist:
        input_filename=dir+"/"+input
        output_filename=dir+"/"+output
        #print "output_filename:"
        #print output_filename
        norm_filename=normlist[i]
        #print "norm_filename:"
        #print norm_filename
        
        simuls.append(normalized_perfect_simulation(input_filename,norm_filename,output_filename))
        simuls[i].description=dir
        simuls[i].species=species
        i=i+1
    return simuls
        
if __name__=='__main__':
    dirlist=sys.argv[1:]
    normlist="norms.namelist"
    species="d,e-"
    species=species.split(',')
    simuls=perfect_simulations_from_dirs(dirlist,normlist,species)
