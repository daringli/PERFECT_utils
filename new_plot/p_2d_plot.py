from p_subplot import perfect_subplot
from p_subplot_group import perfect_subplot_group

from p_visualizer import perfect_visualizer
from perfect_simulations_from_dirs import perfect_simulations_from_dirs
from array_rank import arraylist_rank
from sort_species_list import sort_species_list
from generic_species_labels import generic_species_labels
from is_attribute_species_dependent import is_attribute_species_dependent

from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy
import sys


def perfect_2d_plot(dirlist,attribs,xattr="psi",yattr="theta",normname="norms.namelist",speciesname="species",psiN_to_psiname="psiAHat.h5",global_term_multiplier_name="globalTermMultiplier.h5",cm=cm.rainbow,lg=True,xlims=[90,100],ylims=[0,2],zlims=None,species=True,sort_species=True,first=["D","He"],last=["e"],generic_labels=True,label_dict={"D":"i","H":"i","T":"i","He":"z","N":"z","e":"e"},vlines=[],hlines=[],share_scale=[],skip_species = [],simulList=None,outputname=None,colors=None,linestyles=None,zlabels=None,zaxis_powerlimits=(0,0),pedestal_start_stop=None,pedestal_point=None,core_point=None,putboxes=[]):
    #dirlist: list of simulation directories
    #attribs: list of fields to plot from simulation
    #speciesname: species filename in the simuldir
    #normname: norm filename in the simuldir
    #lg: controls whether to interpret nearby simulations as being paired
    
    if type(attribs) is not list:
        attribs=[attribs]

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
                pedestal_start_stop=(None,None)
            else:
                pedestal_start_stop=(vlines[0],vlines[-1])

        # if no pedestal point is given, we get this from
        # the middle of the pedestal
        if pedestal_point is None:
            if pedestal_start_stop is not None:
                pedestal_point = (pedestal_start_stop[0] + pedestal_start_stop[1])/2.0

        # if no core point is given, we get this from
        # three pedestal widths into the core
        if core_point is None:
            if pedestal_start_stop is not None:
                w = pedestal_start_stop[1] - pedestal_start_stop[0] 
                core_point = pedestal_start_stop[0] - 3*w    

        simulList=perfect_simulations_from_dirs(dirlist, normlist, specieslist, psiN_to_psiList, global_term_multiplierList, pedestal_start_stop_list=pedestal_start_stop, pedestal_point_list=pedestal_point, core_point_list=core_point)
    else:
        "p_2d_plot: simulList specified externally, ignoring dirlist, etc."
    

    #if we translate species labels to generic ion and impurity labels
    #we still need to use original species for sorting preferences
    if generic_labels:
        for simul in simulList:
            simul.species_list=generic_species_labels(simul.species_list,label_dict)
            first=generic_species_labels(first,label_dict)
            last=generic_species_labels(last,label_dict)
            
        
    #get the total list of species found in each simulation
    species_set=set([])
    for simul in simulList:
        # exclude species that are in the skip list
        #nonexcluded = set(simul.species).difference(skip_species)
        nonexcluded = [s for s in simul.species if s not in skip_species]
        
        simul.species_list = nonexcluded
        # add nonexcluded species to total species set
        species_set = species_set | set(simul.species)
        
    # sort species according to preferences
    if sort_species:
        species_set=sort_species_list(list(species_set),first,last)
    print species_set
        
    
    for attrib in attribs:
        psp_list=[] #list of perfect subplot objects
        attrib_sp_dep = is_attribute_species_dependent(simulList,attrib)
        #we will assign data to the following attribute related groups
        attrib_groupname=attrib
        if attrib_sp_dep and species:
            species_attrib_groupname="species_dependent"
            this_species_set=species_set
            gridspec=[len(simulList),len(species_set)]
        else:
            species_attrib_groupname="species_independent"
            this_species_set=set([''])
            gridspec=[len(simulList),1]
        for i,simul in enumerate(simulList):
            data0=getattr(simul,attrib) #need to split this into species data
            #i_s will be coordinate of plot
            for i_s,s in enumerate(this_species_set):
                if attrib_sp_dep:
                    index=[i2 for i2,s2 in enumerate(simul.species) if s2 == s]
                    if len(index)==0:
                        #no data for this species in this simulation. Set data to zeros
                        data=numpy.zeros(data0[:,:,0].shape)
                    elif len(index)==1:
                        data=data0[:,:,index[0]]
                    else:
                        print "perfect_2d_plot: warning: more than one of the same species in the simulation. Will add contributions, but this is untested."
                        data=numpy.sum(data0[:,:,index],axis=2)
                else:
                    data=data0
                subplot_coordinates=(i,i_s)
                #print subplot_coordinates

                
                if i == 0:
                    show_zaxis_ticklabel=True
                    if zlabels is None:
                        title = s
                    else:
                        title = zlabels[i_s]
                else:
                    show_zaxis_ticklabel=False
                    title=''
                    
                if i == len(simulList)-1:
                    show_xaxis_ticklabel=True
                else:
                    show_xaxis_ticklabel=False


                if i_s == 0:
                    show_yaxis_ticklabel=True
                else:
                    show_yaxis_ticklabel=False
                #print simul.local
                if simul.local:
                    gl_grp="local"
                else:
                    gl_grp="global"

                if (yattr == "theta") or (yattr == "theta_shifted"):
                    y_scale = 1/numpy.pi
                    y_period = 2.0
                else:
                    y_scale=1
                    y_period=None
                    
                if (xattr == "theta") or (xattr == "theta_shifted"):
                    x_scale = 1/numpy.pi
                    x_period = 2.0
                elif (xattr == "psi") or (xattr == "actual_psiN"):
                    x_scale=100
                    x_period = None
                else:
                    x_scale=1
                    x_period=None

                if xattr=="psi_o":
                    vlines=simul.pedestal_start_stop_psi_o
                x = getattr(simul,xattr)*x_scale
                y = getattr(simul,yattr)*y_scale
                psp_list.append(perfect_subplot(data,x=x,y=y,subplot_coordinates=subplot_coordinates,show_zaxis_ticklabel=show_zaxis_ticklabel,show_yaxis_ticklabel=show_yaxis_ticklabel,show_xaxis_ticklabel=show_xaxis_ticklabel,title=title,groups=[s,"sim"+str(i),gl_grp,"pair"+str(i/2),species_attrib_groupname],dimensions=2,zaxis_powerlimits=zaxis_powerlimits,x_period=x_period,y_period=y_period))
        #end simulList loop
        for psp in psp_list:
            print psp.groups
        species_groups=[perfect_subplot_group(psp_list,groups=[s,"species_dependent"],logic="and") for s in species_set if len(perfect_subplot_group(psp_list,groups=[s,"species_dependent"],logic="and").p_subplot_list)>0]
        nospecies_group=perfect_subplot_group(psp_list,groups=["species_independent"])

        sim_groups = [perfect_subplot_group(psp_list,groups=["sim"+str(i)]) for i in range(len(simulList))]
        pair_groups = [perfect_subplot_group(psp_list,groups=["pair"+str(i)]) for i in range(len(simulList)/2)]
        local_group = perfect_subplot_group(psp_list,groups=["local"])
        global_group = perfect_subplot_group(psp_list,groups=["global"])
        all_group=perfect_subplot_group(psp_list,groups='',get_all=True)
        share_scale_group=perfect_subplot_group(psp_list,groups=share_scale,logic="or")
        
        if lg==False:
            if colors is None:
                color=iter(cm(numpy.linspace(0,1,len(simulList))))
            else:
                color=iter(colors)
            for sim_group in sim_groups:
                c=next(color)
                sim_group.setattrs("border_color",c)
        else:
            if colors is None:
                color=iter(cm(numpy.linspace(0,1,len(pair_groups))))
            else:
                color=iter(colors)
            for pair_group in pair_groups:
                c=next(color)
                pair_group.setattrs("border_color",c)

        if linestyles is None:
            local_group.setattrs("border_linestyle","dashed")
        else:
            for isg,grp in enumerate(sim_groups):
                grp.setattrs("border_linestyle",linestyles[isg])
            
        all_group.setattrs("xlims",xlims)
        #all_group.setattrs("ylims",[all_group.get_min("y",xlim=False,ylim=False),all_group.get_max("y",xlim=False,ylim=False)])
        all_group.setattrs("ylims",ylims)
        all_group.setattrs("vlines",vlines)
        all_group.setattrs("hlines",hlines)

        if zlims is None:
            for species_group in species_groups:
                species_group.setattrs("zlims",[species_group.get_min("data"),species_group.get_max("data")])
                print [species_group.get_min("data"),species_group.get_max("data")]

            if len(share_scale_group.p_subplot_list)>0:
                share_scale_group.setattrs("zlims",[share_scale_group.get_min("data"),share_scale_group.get_max("data")])

            if len(nospecies_group.p_subplot_list)>0:
                nospecies_group.setattrs("zlims",[nospecies_group.get_min("data"),nospecies_group.get_max("data")])
        else:
            all_group.setattrs("zlims",zlims)
            
        for i,sim_group in enumerate(sim_groups):
            if i != 0:
                sim_group.setattrs("show_zaxis",False)
        
                
        if xattr is not "psi_o":
            global_xlabel=r"$100\psi_N$"
        else:
            global_xlabel=r"$\psi^{\circ}$"
        perfect_visualizer(psp_list,gridspec,global_xlabel=global_xlabel,dimensions=2,global_ylabel=r"$\theta/\pi$",putboxes=putboxes)

        if outputname is None:
            plt.savefig(attrib+'.pdf')
        else:
            plt.savefig(outputname+'.pdf')

if __name__=="__main__":
    dirlist=sys.argv[1:]
    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
    attribs=["density_perturbation","flow"]
    perfect_2d_plot(dirlist,attribs)
