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


def perfect_2d_plot(dirlist,attribs,normname="norms.namelist",speciesname="species",cm=cm.rainbow,lg=True,xlims=[90,100],ylims=[0,2],species=True,sort_species=True,first=["D","He"],last=["e"],generic_labels=True,label_dict={"D":"i","He":"z","N":"z","e":"e"},vlines=None,hlines=None):
    #dirlist: list of simulation directories
    #attribs: list of fields to plot from simulation
    #speciesname: species filename in the simuldir
    #normname: norm filename in the simuldir
    #lg: controls whether to interpret nearby simulations as being paired
    
    if type(attribs) is not list:
        attribs=[attribs]
    
    normlist=[x + "/" + normname for x in dirlist]
    specieslist=[x + "/" + speciesname for x in dirlist]
    simulList=perfect_simulations_from_dirs(dirlist,normlist,specieslist)

    if generic_labels:
        for simul in simulList:
            simul.species_list=generic_species_labels(simul.species_list,label_dict)
            first=generic_species_labels(first,label_dict)
            last=generic_species_labels(last,label_dict)
            
        
    #get number of columns needed to fit all species
    species_set=set([])
    for simul in simulList:
        species_set = species_set | set(simul.species)
    if sort_species:
        species_set=sort_species_list(list(species_set),first,last)
    
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
            if len(simul.species) == 1:
                #adding a species dimension to single-species data
                data0=data0[:,:,numpy.newaxis]
            for i_s,s in enumerate(this_species_set):
                if attrib_sp_dep:
                    index=[i2 for i2,s2 in enumerate(simul.species) if s2 == s]
                    if len(index)==0:
                        #no data for this species in this simulation. Set data to zeros
                        data=numpy.zeros(data0[:,:,0].shape)
                    elif len(index)==1:
                        data=data0[:,:,index[0]]
                    else:
                        print "perfect_2d_plot: warning: more than one of the same species in the simulation. Will add contributions."
                        data=numpy.sum(data0[:,:,index],axis=2)
                else:
                    data=data0
                #print data
                subplot_coordinates=(i,i_s)
                #print subplot_coordinates 
                if i == 0:
                    show_zaxis_ticklabel=True
                    show_xaxis_ticklabel=False
                    title=s
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
                #print data
                psp_list.append(perfect_subplot(data,x=100*simul.psi,y=simul.theta/numpy.pi,subplot_coordinates=subplot_coordinates,show_zaxis_ticklabel=show_zaxis_ticklabel,show_yaxis_ticklabel=show_yaxis_ticklabel,show_xaxis_ticklabel=show_xaxis_ticklabel,title=title,groups=[s,"sim"+str(i),gl_grp,"pair"+str(i/2),species_attrib_groupname],dimensions=2))
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

        
        color=iter(cm(numpy.linspace(0,1,len(pair_groups))))
        
        if lg==False:
            for sim_group in sim_groups:
                c=next(color)
                sim_group.setattrs("border_color",c)
        else:
            for pair_group in pair_groups:
                c=next(color)
                pair_group.setattrs("border_color",c)

        local_group.setattrs("border_linestyle","dashed")
        
        all_group.setattrs("xlims",xlims)
        #all_group.setattrs("ylims",[all_group.get_min("y",xlim=False,ylim=False),all_group.get_max("y",xlim=False,ylim=False)])
        all_group.setattrs("ylims",ylims)
        all_group.setattrs("vlines",vlines)
        all_group.setattrs("hlines",hlines)
        
        for species_group in species_groups:
            species_group.setattrs("zlims",[species_group.get_min("data"),species_group.get_max("data")])

        if len(nospecies_group.p_subplot_list)>0:
            nospecies_group.setattrs("zlims",[nospecies_group.get_min("data"),nospecies_group.get_max("data")])
            
        for i,sim_group in enumerate(sim_groups):
            if i != 0:
                sim_group.setattrs("show_zaxis",False)
        

            
        perfect_visualizer(psp_list,gridspec,global_xlabel=r"$100\psi_N$",dimensions=2,global_ylabel=r"$\theta/\pi$")

        plt.savefig(attrib+'.pdf')

if __name__=="__main__":
    dirlist=sys.argv[1:]
    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
    attribs=["density_perturbation","flow"]
    perfect_2d_plot(dirlist,attribs)
