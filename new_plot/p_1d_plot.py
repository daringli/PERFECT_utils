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
import collections

import scipy.integrate

from mpldatacursor import datacursor

def perfect_1d_plot(dirlist,attribs,xattr="psi",normname="norms.namelist",speciesname="species",psiN_to_psiname="psiAHat.h5",global_term_multiplier_name="globalTermMultiplier.h5",cm=cm.rainbow,lg=True,markers=None,markeverys=None,markersizes=None,fillstyles=None,linestyles=None,linewidths=None,xlims=None,same_plot=False,outputname="default",ylabels=None,label_all=False,global_ylabel="",sort_species=True,first=["D","He"],last=["e"],generic_labels=True,label_dict={"D":"i","H":"i","T":"i","He":"z","N":"z","e":"e"},vlines=[],hlines=[],share_scale=[],interactive=False,colors=None,skip_species = [],yaxis_powerlimits=(0,0),hidden_xticklabels=[],yaxis_label_x=-0.15,yaxis_label_y=0.5,ylims=None,simulList=None,pedestal_start_stop=None,pedestal_point=None,core_point=None,putboxes=[],return_psps=False,set_species_title=True,group_mode="species",visualizer_rows=None,label_size=None,filetype=".pdf",legend=None):
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

    if label_size is not None:
        yaxis_label_size = label_size
        xaxis_label_size = label_size
    else:
        yaxis_label_size = 12
        xaxis_label_size = 12
        
    all_linestyles = linestyles
                
    if simulList is None:
        normlist=[x + "/" + normname for x in dirlist]
        specieslist=[x + "/" + speciesname for x in dirlist]
        psiN_to_psiList=[x + "/" + psiN_to_psiname for x in dirlist]
        global_term_multiplierList=[x + "/" + global_term_multiplier_name for x in dirlist]

        # If no pedestal start stop is given:
        #deduce pedestal start stop from vlines
        #(this is for backwards compatiblity)
        if pedestal_start_stop is None:
            if vlines is None:
                pedestal_start_stop=[]
            else:
                pedestal_start_stop=(vlines[0],vlines[-1])
            
                
        # if no pedestal point is given, we get this from
        # the middle of the pedestal
        if pedestal_point is None:
            if pedestal_start_stop is not None:
                if isinstance(pedestal_start_stop[0], collections.Iterable):
                    pedestal_point = [(pss[0] + pss[1])/2.0 for pss in pedestal_start_stop]
                else:
                    pedestal_point = (pedestal_start_stop[0] + pedestal_start_stop[1])/2.0
        # if no core point is given, we get this from
        # three pedestal widths into the core
        if core_point is None:
            if pedestal_start_stop is not None:
                if isinstance(pedestal_start_stop[0], collections.Iterable):
                    w = [pss[1] - pss[0] for pss in pedestal_start_stop]
                    core_point = [pss[0] - 3*(pss[1] - pss[0]) for pss in pedestal_start_stop] 
                else:
                    w = pedestal_start_stop[1] - pedestal_start_stop[0] 
                    core_point = pedestal_start_stop[0] - 3*w    

        simulList=perfect_simulations_from_dirs(dirlist, normlist, specieslist, psiN_to_psiList, global_term_multiplierList, pedestal_start_stop_list=pedestal_start_stop, pedestal_point_list=pedestal_point, core_point_list=core_point)      
    else:
        "p_1d_plot: simulList specified externally, ignoring dirlist, etc."

    for simul in simulList:
        pass
        #print "max deltaT: " + str(simul.max_deltaT_insideLCFS)
        #print "max deltaN: " +str(simul.max_deltaN_insideLCFS)
        #print "max deltaEta: " +str(simul.max_deltaEta_insideLCFS)
        #print "max U: " +str(simul.max_U_insideLCFS)
        #print "nuHat_97: " + str(simul.collisionality_at_point)
        #print "epsilon32_97: " + str(simul.epsilon32_at_point)
        
    if markers is not None:
        if markeverys is None:
            markeverys = [None]*len(simulList)
        if markersizes is None:
            markersizes = [1]*len(simulList)



        if fillstyles is None:
            if lg:
                fillstyles=[]
                for n,l in zip(noddpsi,local):
                    # local overrides noddpsi in code as well
                    if l:
                        fillstyles=all_linestyles+["none"]
                    elif n:
                        fillstyles=all_linestyles+["left"]
                    else:
                        fillstyles=all_linestyles+["full"]
            else:
                fillstyles = ['full']*len(simulList)
    else:
        markers = [None]*len(simulList)
        fillstyles = ["none"]*len(simulList)
        markeverys = [None]*len(simulList)
        markersizes = [None]*len(simulList)


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
        nonexcluded = [s for s in simul.species if s not in skip_species]
        #simul.species_list = nonexcluded
        # add nonexcluded species to total species set
        species_set = species_set | set(nonexcluded)
    if sort_species:
        species_set=sort_species_list(list(species_set),first,last)
    
    i=-1
    psp_lists=[] #list of lists of perfect subplot object
    gridspec_list=[]

    #assign colors to simulations
    color_index=0

    local=[simul.local for simul in simulList]
    noddpsi=[simul.no_ddpsi for simul in simulList]

    
    if colors is None:
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

    if all_linestyles is None:
        if lg:
            all_linestyles=[] #linestyles generated from local or global
            for n,l in zip(noddpsi,local):
                # local overrides noddpsi in code as well
                if l:
                    all_linestyles=all_linestyles+["dashed"]
                elif n:
                    all_linestyles=all_linestyles+["dashdot"]
                else:
                    all_linestyles=all_linestyles+["solid"]
        else:
            all_linestyles = ['solid']*len(simulList)
                
    
                    

    if linewidths is None:
        linewidths = [1]*len(simulList)
    
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
        if attrib_sp_dep:
            this_iterable = None
            if group_mode == "species":
                this_iterable = species_set
            elif group_mode == "simulations":
                this_iterable = range(len(simulList))
                set_species_title = False
       
            for i_sp,s in enumerate(this_iterable):
                i=i+1
                data = None
                
                if group_mode == "species":
                    #data is taken for a given species for all simulations
                    #index of species in simulation given by index to index
                    index=[[ind for ind,spec in enumerate(simul.species) if spec==s] for simul in simulList]
                    if all(len(ind)<=1 for ind in index):
                        data=[getattr(simul,attrib)[:,index[i_si][0]] for i_si,simul in enumerate(simulList) if s in simul.species]
                    else:
                        print "p_1d_plot: warning: more than one of the same species in the simulation. Will add contributions."
                        data=[numpy.sum(getattr(simul,attrib)[:,index[i_si]],axis=1) for i_si,simul in enumerate(simulList) if s in simul.species]
                        
                elif group_mode == "simulations":
                    # i_sp indices simulation and i_si species
                    data=[getattr(simulList[i_sp],attrib)[:,i_si] for i_si in range(len(simul.species))]

                if xattr=="theta":
                    x_scale=1/numpy.pi
                    x_period = 2.0
                else:
                    x_scale=1
                    x_period = None
                    
                if xattr != None:
                    if group_mode == "species":
                        x=[getattr(simul,xattr)*x_scale for simul in simulList if s in simul.species]
                    elif group_mode == "simulations":
                        x=[getattr(simulList[i_sp],xattr)*x_scale for s in simulList[i_sp].species]
                else:
                    # If xattrib is None, we plot against the index of the data
                    # This probably will not work if we are not plotting against
                    # the first index of data
                    if group_mode == "species":
                        x=[numpy.array(range(len(getattr(simul,attrib)))) for simul in simulList if s in simul.species]
                    elif group_mode == "simulations":
                        x=[numpy.array(range(len(getattr(simulList[i_sp],attrib)))) for s in simulList[i_sp].species]
                if xlims is None:
                    # min to max among all the simulations
                    xlims = [numpy.min(x),numpy.max(x)]

                if group_mode == "species":
                    linecolors=[all_linecolors[i_si] for i_si,simul in enumerate(simulList) if s in simul.species]
                    linestyles=[all_linestyles[i_si] for i_si,simul in enumerate(simulList) if s in simul.species]
                elif group_mode == "simulations":
                    linecolors=[all_linecolors[i_si] for i_si in range(len(simulList[i_sp].species))]
                    linestyles=[all_linestyles[i_si] for i_si in range(len(simulList[i_sp].species))]
                coordinates=(i,0)

                if perhaps_last and (i_sp == len(this_iterable) - 1):
                    last_groupname="last"
                    gridspec_list.append([i+1,1])
                else:
                    last_groupname="not_last"

                if ylims is None:
                    ylim = None
                else:
                    if isinstance(ylims[0], collections.Iterable):
                        ylim = ylims[i]
                    else:
                        ylim=ylims

                psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,dimensions=1,groups=[s,attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,linewidths=linewidths,colors=linecolors,markers=markers,fillstyles=fillstyles,markeverys=markeverys,markersizes=markersizes,yaxis_powerlimits=yaxis_powerlimits,hidden_xticklabels=hidden_xticklabels,yaxis_label_x=yaxis_label_x,yaxis_label_y=yaxis_label_y,ylims=ylim,x_period=x_period,yaxis_label_size=yaxis_label_size,xaxis_label_size=xaxis_label_size)) #yaxis_label=ylabels[i_a]
                

        else:
            i=i+1
            if perhaps_last:
                last_groupname="last"
                gridspec_list.append([i+1,1])
            else:
                last_groupname="not_last"
            #species independent plot
            data=[getattr(simul,attrib) for simul in simulList]

            if xattr=="theta":
                x_scale=1/numpy.pi
                x_period = 2.0
            else:
                x_scale=1
                x_period = None

            if xattr != None:
                x=[x_scale*getattr(simul,xattr) for simul in simulList]

            else:
                # If xattrib is None, we plot against the index of the data
                # This probably will not work if we are not plotting against
                # the first index of data
                try:
                    L=len(getattr(simul,attrib))
                except TypeError:
                    # does not have a length, so probably a single element
                    L=1
                
                x=[numpy.array(range(L)) for simul in simulList]
            if xlims is None:
                # min to max among all the simulations
                xlims = [numpy.min(x),numpy.max(x)]
            
            

            linecolors=all_linecolors
            linestyles=all_linestyles
            coordinates=(i,0)
            
            if ylims is None:
                ylim = None
            else:
                if isinstance(ylims[0], collections.Iterable):
                    ylim = ylims[i]
                else:
                    ylim=ylims
            
            print ylim
            psp_list.append(perfect_subplot(data,x,subplot_coordinates=coordinates,dimensions=1,groups=[attrib_groupname,species_attrib_groupname,last_groupname],linestyles=linestyles,colors=linecolors,markers=markers,fillstyles=fillstyles,markeverys=markeverys,yaxis_powerlimits=yaxis_powerlimits,hidden_xticklabels=hidden_xticklabels,yaxis_label_x=yaxis_label_x,ylims=ylim,x_period=x_period,yaxis_label_size=yaxis_label_size,xaxis_label_size=xaxis_label_size)) #yaxis_label=ylabels[i_a]
            
        
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
                    if ylims is None:
                        share_scale_group.setattrs("ylims",[share_scale_group.get_min("data",margin=0.1),share_scale_group.get_max("data",margin=0.1)])

                    
        
        species_groups=[perfect_subplot_group(psp_list,groups=[s]) for s in species_set]
        species_indep_groups=perfect_subplot_group(psp_list,groups=["species_independent"])
        local_group = perfect_subplot_group(psp_list,groups=["local"])
        global_group = perfect_subplot_group(psp_list,groups=["global"])
        last_group = perfect_subplot_group(psp_list,groups=["last"])
        all_group=perfect_subplot_group(psp_list,groups='',get_all=True)

        if set_species_title:
            for species_group,s in zip(species_groups,species_set):
                species_group.setattrs("title",s)

        

        for attrib_group in attrib_groups:
            this_species_groups=[perfect_subplot_group(attrib_group.p_subplot_list,groups=[s]) for s in species_set]
            for this_species_group in this_species_groups:
                if len(this_species_group.p_subplot_list)>0:
                    if ylims is None:
                        this_species_group.setattrs("ylims",[this_species_group.get_min("data",margin=0.1),this_species_group.get_max("data",margin=0.1)])
                    #print [this_species_group.get_min("data",margin=0.1),this_species_group.get_max("data",margin=0.1)]
                    
            this_species_indep_group=perfect_subplot_group(attrib_group.p_subplot_list,groups=["species_independent"])
            if len(this_species_indep_group.p_subplot_list)>0:
                if ylims is None:
                    this_species_indep_group.setattrs("ylims",[this_species_indep_group.get_min("data",margin=0.1),this_species_indep_group.get_max("data",margin=0.1)])
            this_share_scale_group=perfect_subplot_group(attrib_group.p_subplot_list,groups=share_scale,logic="or")
            if len(this_share_scale_group.p_subplot_list)>0:
                if ylims is None:
                    this_share_scale_group.setattrs("ylims",[this_share_scale_group.get_min("data",margin=0.1),this_share_scale_group.get_max("data",margin=0.1)])
            
        if xattr=="psi":
            global_xlabel=r"$\psi_N$"
        elif xattr=="psiN2":
            global_xlabel=r"$\psi_N$"
        elif xattr=="psiN3":
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
            global_xlabel=r"$\psi^{\circ}$"
        elif xattr==None:
            global_xlabel=r"$i$"

        if xattr == "psi_o":
            this_vlines = simulList[0].pedestal_start_stop_psi_o
        else:
            this_vlines = vlines

        all_group.setattrs("show_yaxis_ticklabel",True)
        all_group.setattrs("vlines",this_vlines)
        all_group.setattrs("hlines",hlines)
        all_group.setattrs("show_xaxis_ticklabel",False)
        last_group.setattrs("show_xaxis_ticklabel",True)

        if return_psps:
            return psp_list

        print gridspec_list[i_li]
        if visualizer_rows is None:
            visualizer_rows = gridspec_list[i_li][0]
        
        perfect_visualizer(psp_list,gridspec_list[i_li],global_xlabel=global_xlabel,dimensions=1,global_ylabel=global_ylabel,interactive=interactive,putboxes=putboxes,rows=visualizer_rows,global_xlabel_size=xaxis_label_size,legend=legend)
        if same_plot:
            plt.savefig(outputname+filetype)
        else:
            plt.savefig(attribs[i_li]+filetype)

    
        
if __name__=="__main__":
    dirlist=sys.argv[1:]
    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
    attribs=["particle_flux","heat_flux","momentum_flux","Phi"]
    perfect_1d_plot(dirlist,attribs,same_plot="True")
