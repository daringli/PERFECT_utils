from p_subplot import perfect_subplot
from p_subplot_group import perfect_subplot_group

from p_visualizer import perfect_visualizer
from perfect_simulations_from_dirs import perfect_simulations_from_dirs
from array_rank import arraylist_rank
from is_attribute_species_dependent import is_attribute_species_dependent
from sort_species_list import sort_species_list
from generic_species_labels import generic_species_labels

from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy
import sys

import scipy.integrate



def perfect_1d_plot(dirlist,attribs,xattr="psi",normname="norms.namelist",speciesname="species",cm=cm.rainbow,lg=True,xlims=[0.9,1.0],same_plot=False,outputname="default",ylabels=None,label_all=False,global_ylabel="",sort_species=True,first=["D","He"],last=["e"],generic_labels=True,label_dict={"D":"i","He":"z","N":"z","e":"e"},vlines=None,hlines=None):
    #dirlist: list of simulation directories
    #attribs: list of fields to plot from simulation
    #speciesname: species filename in the simuldir
    #normname: norm filename in the simuldir
    #lg: controls whether to interpret nearby simulations as being paired
    
    if type(attribs) is not list:
        attribs=[attribs]
    if ylabels is not None:
        if type(ylabels) is not list:
            ylabels=[ylabels]*len(attribs)
        else:
            if len(ylabels) != len(attribs):
                print "p_1d_plot: error: ylabels not the same size as attribs"
                exit(1)
    else:
        ylabels=['']*len(attribs)
                
    
    normlist=[x + "/" + normname for x in dirlist]
    specieslist=[x + "/" + speciesname for x in dirlist]
    simulList=perfect_simulations_from_dirs(dirlist,normlist,specieslist)

    #some logic to differentiate between species independent
    #and species quantities further down will fail if there are simulations
    #one species (which you really shouldn't do due to quasi-neutrality!!)
    num_species_array=numpy.array([len(simul.species) for simul in simulList])
    one_species_list=numpy.where(num_species_array<=1)
    if len(num_species_array[one_species_list])>0:
        print "p_1d_plot: warning: there are simulations with one (or) less species! Logic to determine whether attribute is a species property will not work."
        
    
    if generic_labels:
        for simul in simulList:
            simul.species_list=generic_species_labels(simul.species_list,label_dict)
            first=generic_species_labels(first,label_dict)
            last=generic_species_labels(last,label_dict)
    species_set=set([])
    for simul in simulList:
        species_set = species_set | set(simul.species)
    
    #print gridspec
    

    i=-1
    psp_lists=[] #list of lists of perfect subplot object
    gridspec_list=[]
    local=[simul.local for simul in simulList]
    colors=cm(numpy.linspace(0,1,len(simulList)-sum(local)))
    #assign colors to simulations
    all_linecolors=[]
    color_index=0
    if lg:
        for loc in local:
            all_linecolors.append(colors[color_index])
            if loc == True:
                color_index=color_index+1
    else:
        all_linecolors=colors
        
    linestyles=[] #linestyles generated from local or global
    for l in local:
        if l:
            linestyles=linestyles+["dashed"]
        else:
            linestyles=linestyles+["solid"]
    
    for i_a,attrib in enumerate(attribs):
        psp_list=[]
        #check whether attribute is species dependent
        #it will be assumed to be if the lenght along the 1 axis
        #is equal to the number of species in a simulation for all simulations
        #at least as long as there is more than one species in the simulation.
        attrib_sp_dep = is_attribute_species_dependent(simulList,attrib)
        #we will assign data to the following attribute related groups
        attrib_groupname=attrib
        if attrib_sp_dep:
            species_attrib_groupname="species_dependent"
        else:
            species_attrib_groupname="species_independent"
        if same_plot:
            if i_a == len(attribs)-1:
                perhaps_last=True
            else:
                perhaps_last=False
        else:
            perhaps_last=True
        if sort_species:
            species_set=sort_species_list(list(species_set),first,last)
        if attrib_sp_dep:
            for i_sp,s in enumerate(species_set):
                i=i+1
                #data is taken for a given species for all simulations
                #index of species in simulation given by index to index
                index=[[ind for ind,spec in enumerate(simul.species) if spec==s] for simul in simulList]
                
                if all(len(ind)<=1 for ind in index):
                    data=[getattr(simul,attrib)[:,index[i_si][0]] for i_si,simul in enumerate(simulList) if s in simul.species]
                else:
                    print "p_1d_plot: warning: more than one of the same species in the simulation. Will add contributions."
                    data=[numpy.sum(getattr(simul,attrib)[:,index[i_si]],axis=1) for i_si,simul in enumerate(simulList) if s in simul.species]
                
                x=[getattr(simul,xattr) for simul in simulList if s in simul.species]
                linecolors=[all_linecolors[i_si] for i_si,simul in enumerate(simulList) if s in simul.species]
                coordinates=(i,0)
                if perhaps_last and (i_sp == len(species_set) - 1):
                    last_groupname="last"
                    gridspec_list.append([i+1,1])
                else:
                    last_groupname="not_last"

                psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,groups=[s,attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,colors=linecolors)) #yaxis_label=ylabels[i_a]
                

        else:
            i=i+1
            if perhaps_last:
                last_groupname="last"
                gridspec_list.append([i+1,1])
            else:
                last_groupname="not_last"
            #species independent plot
            data=[getattr(simul,attrib) for simul in simulList]
            x=[getattr(simul,xattr) for simul in simulList]
            linecolors=all_linecolors
            coordinates=(i,0)
            psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,groups=[attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,colors=linecolors)) #yaxis_label=ylabels[i_a]
            
        
        psp_lists.append(psp_list)
        if not same_plot:
            i=-1

    #merge the psp_lists if everything is supposed to go in the same plot
    if same_plot:
        final_psp_lists=[]
        for psp_list in psp_lists:
            final_psp_lists = final_psp_lists + psp_list
        psp_lists=[final_psp_lists]
        
            
    for i_li,psp_list in enumerate(psp_lists):
        for psp in psp_list:
            print psp.groups
            psp.xlims=xlims
            psp.data=psp.data_inrange()
            psp.x=psp.x_inrange()
        if same_plot:
            attrib_groups=[perfect_subplot_group(psp_list,groups=[a]) for a in attribs]
            for ylabel,attrib_group in zip(ylabels,attrib_groups):
                if label_all:
                    attrib_group.setattrs("yaxis_label",ylabel)
                else:
                    attrib_group.set_middle_attr("yaxis_label",ylabel)
        else:
            attrib_groups=[perfect_subplot_group(psp_list,groups=[a]) for a in [attribs[i_li]]]
            for ylabel,attrib_group in zip([ylabels[i_li]],attrib_groups):
                if label_all:
                    attrib_group.setattrs("yaxis_label",ylabel)
                else:
                    attrib_group.set_middle_attr("yaxis_label",ylabel)
        
        species_groups=[perfect_subplot_group(psp_list,groups=[s]) for s in species_set]
        species_indep_groups=perfect_subplot_group(psp_list,groups=["species_independent"])
        local_group = perfect_subplot_group(psp_list,groups=["local"])
        global_group = perfect_subplot_group(psp_list,groups=["global"])
        last_group = perfect_subplot_group(psp_list,groups=["last"])
        all_group=perfect_subplot_group(psp_list,groups='',get_all=True)
        
        for species_group,s in zip(species_groups,species_set):
            species_group.setattrs("title",s)
            
        all_group.setattrs("show_yaxis_ticklabel",True)
        all_group.setattrs("vlines",vlines)
        all_group.setattrs("hlines",hlines)
        last_group.setattrs("show_xaxis_ticklabel",True)
        #print gridspec_list[i_li]
        if xattr=="psi":
            global_xlabel=r"$\psi_N$"
        elif xattr=="theta":
            global_xlabel=r"$\theta$"
            
        perfect_visualizer(psp_list,gridspec_list[i_li],global_xlabel=global_xlabel,dimensions=1,global_ylabel=global_ylabel)
        if same_plot:
            plt.savefig(outputname+".pdf")
        else:
            plt.savefig(attribs[i_li]+".pdf")
        

    
        
if __name__=="__main__":
    dirlist=sys.argv[1:]
    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
    attribs=["particle_flux","heat_flux","momentum_flux","Phi"]
    perfect_1d_plot(dirlist,attribs,same_plot="True")
