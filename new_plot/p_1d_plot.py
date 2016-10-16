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

from mpldatacursor import datacursor




def perfect_1d_plot(dirlist,attribs,xattr="psi",normname="norms.namelist",speciesname="species",psiN_to_psiname="psiAHat.h5",global_term_multiplier_name="globalTermMultiplier.h5",cm=cm.rainbow,lg=True,markers=None,linestyles=None,xlims=None,same_plot=False,outputname="default",ylabels=None,label_all=False,global_ylabel="",sort_species=True,first=["D","He"],last=["e"],generic_labels=True,label_dict={"D":"i","He":"z","N":"z","e":"e"},vlines=None,hlines=None,share_scale=[],interactive=False,colors=None):
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
    psiN_to_psiList=[x + "/" + psiN_to_psiname for x in dirlist]
    global_term_multiplierList=[x + "/" + global_term_multiplier_name for x in dirlist]
    if vlines== None:
        here_vlines=[]
    else:
        here_vlines=vlines
    simulList=perfect_simulations_from_dirs(dirlist,normlist,specieslist,psiN_to_psiList,global_term_multiplierList,pedestal_points_list=here_vlines)
    if markers == None:
        markers=['']*len(simulList)

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

    #assign colors to simulations
    color_index=0

    local=[simul.local for simul in simulList]
    noddpsi=[simul.no_ddpsi for simul in simulList]

    
    if colors == None:
        colors=cm(numpy.linspace(0,1,len(simulList)))
        all_linecolors=[]
    
        if lg:
            colors=cm(numpy.linspace(0,1,len(simulList)-sum(local)-sum(noddpsi)))
            if sum(local)>0:
                #if we have local simulations, increment color after them
                for loc in local:
                    all_linecolors.append(colors[color_index])
                    if loc == True:
                        color_index=color_index+1
            elif sum(noddpsi)>0:
                #otherwise, increment color after noddpsi simulation
                for nd in noddpsi:
                    all_linecolors.append(colors[color_index])
                    if nd == True:
                        color_index=color_index+1
            else:
                #else, always increment
                all_linecolors=cm(numpy.linspace(0,1,len(simulList)))
        else:
            # force always increment behavior
            all_linecolors=cm(numpy.linspace(0,1,len(simulList)))
    else:
        all_linecolors=colors

    if linestyles == None:
        linestyles=[] #linestyles generated from local or global
        for n,l in zip(noddpsi,local):
            # local overrides noddpsi in code as well
            if l:
                linestyles=linestyles+["dashed"]
            elif n:
                linestyles=linestyles+["dashdot"]
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

                if xattr=="theta":
                    x_scale=1/numpy.pi
                else:
                    x_scale=1
                if xattr != None:
                    x=[getattr(simul,xattr)*x_scale for simul in simulList if s in simul.species]
                else:
                    # If xattrib is None, we plot against the index of the data
                    # This probably will not work if we are not plotting against
                    # the first index of data
                    x=[numpy.array(range(len(getattr(simul,attrib)))) for simul in simulList if s in simul.species]
                if xlims == None:
                    # min to max among all the simulations
                    xlims = [numpy.min(x),numpy.max(x)]
                    
                linecolors=[all_linecolors[i_si] for i_si,simul in enumerate(simulList) if s in simul.species]
                coordinates=(i,0)
                if perhaps_last and (i_sp == len(species_set) - 1):
                    last_groupname="last"
                    gridspec_list.append([i+1,1])
                else:
                    last_groupname="not_last"

                psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,groups=[s,attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,colors=linecolors,markers=markers)) #yaxis_label=ylabels[i_a]
                

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
            psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,groups=[attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,colors=linecolors,markers=markers)) #yaxis_label=ylabels[i_a]
            
        
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
                    attrib_group.set_middle_ylabel(ylabel)
        else:
            attrib_groups=[perfect_subplot_group(psp_list,groups=[a]) for a in [attribs[i_li]]]
            for ylabel,attrib_group in zip([ylabels[i_li]],attrib_groups):
                if label_all:
                    attrib_group.setattrs("yaxis_label",ylabel)
                else:
                    attrib_group.set_middle_ylabel(ylabel)
                if len(share_scale)>0:
                    share_scale_group=perfect_subplot_group(psp_list,groups=share_scale,logic="or")
                    share_scale_group.setattrs("ylims",[share_scale_group.get_min("data",margin=0.1),share_scale_group.get_max("data",margin=0.1)])
                    

                    
        
        species_groups=[perfect_subplot_group(psp_list,groups=[s]) for s in species_set]
        species_indep_groups=perfect_subplot_group(psp_list,groups=["species_independent"])
        local_group = perfect_subplot_group(psp_list,groups=["local"])
        global_group = perfect_subplot_group(psp_list,groups=["global"])
        last_group = perfect_subplot_group(psp_list,groups=["last"])
        all_group=perfect_subplot_group(psp_list,groups='',get_all=True)
        
        for species_group,s in zip(species_groups,species_set):
            species_group.setattrs("title",s)

        

        for attrib_group in attrib_groups:
            this_species_groups=[perfect_subplot_group(attrib_group.p_subplot_list,groups=[s]) for s in species_set]
            for this_species_group in this_species_groups:
                if len(this_species_group.p_subplot_list)>0:
                    this_species_group.setattrs("ylims",[this_species_group.get_min("data",margin=0.1),this_species_group.get_max("data",margin=0.1)])
                    #print [this_species_group.get_min("data",margin=0.1),this_species_group.get_max("data",margin=0.1)]
            this_species_indep_group=perfect_subplot_group(attrib_group.p_subplot_list,groups=["species_independent"])
            if len(this_species_indep_group.p_subplot_list)>0:
                this_species_indep_group.setattrs("ylims",[this_species_indep_group.get_min("data",margin=0.1),this_species_indep_group.get_max("data",margin=0.1)])
            this_share_scale_group=perfect_subplot_group(attrib_group.p_subplot_list,groups=share_scale,logic="or")
            if len(this_share_scale_group.p_subplot_list)>0:
                this_share_scale_group.setattrs("ylims",[this_share_scale_group.get_min("data",margin=0.1),this_share_scale_group.get_max("data",margin=0.1)])
            
            
        
            
        #print gridspec_list[i_li]
        if xattr=="psi":
            global_xlabel=r"$\psi_N$"
        elif xattr=="theta":
            global_xlabel=r"$\theta/\pi$"
        elif xattr=="sqrtpsi":
            global_xlabel=r"$\sqrt{\psi_N}$"
        elif xattr=="psiOverOrbitWidth":
            global_xlabel=r"$\sqrt{\psi_N}$"
        elif xattr=="actual_psi":
            global_xlabel=r"$\hat{\psi}$"
        elif xattr=="actual_psiN":
            global_xlabel=r"$\psi_N$"
        elif xattr=="psi_index":
            global_xlabel=r"$i_\psi$"
        elif xattr=="psi_o":
            global_xlabel=r"$\psi^o$"
        elif xattr==None:
            global_xlabel=r"$i$"
            
        if xattr == "psi_o":
            this_vlines = simulList[0].pedestal_points_psi_o
        else:
            this_vlines = vlines
        all_group.setattrs("show_yaxis_ticklabel",True)
        all_group.setattrs("vlines",this_vlines)
        all_group.setattrs("hlines",hlines)
        last_group.setattrs("show_xaxis_ticklabel",True)
        
        perfect_visualizer(psp_list,gridspec_list[i_li],global_xlabel=global_xlabel,dimensions=1,global_ylabel=global_ylabel)
        if same_plot:
            plt.savefig(outputname+".pdf")
        else:
            plt.savefig(attribs[i_li]+".pdf")
        if interactive:
            #dangerous, since it will (for some reason) be executed after all 1d_plot calls and show everything plotted in the given script.
            datacursor(display='multiple', draggable=True)
            plt.show()
        

    
        
if __name__=="__main__":
    dirlist=sys.argv[1:]
    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
    attribs=["particle_flux","heat_flux","momentum_flux","Phi"]
    perfect_1d_plot(dirlist,attribs,same_plot="True")
