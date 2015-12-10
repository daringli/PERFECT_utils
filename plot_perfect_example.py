from plot_object import plot_object
from perfect_simulation import perfect_simulation
from perfect_simulations_from_dirs import perfect_simulations_from_dirs
import sys
import matplotlib.pyplot as plt



if __name__=='__main__':
    dirlist=sys.argv[1:]
    normlist="norms.namelist"
    species="d,N,e"
    species=species.split(',')
    #print species
    simulList=perfect_simulations_from_dirs(dirlist,normlist,species)
    #print simulList

    plotObjList=[]
    plotObjList.append(plot_object("normed_particle_flux",len(species),1,"particle flux"))
    plotObjList.append(plot_object("normed_conductive_heat_flux",len(species),1,"heat flux"))
    plotObjList.append(plot_object("T",len(species),1,"T"))
    plotObjList.append(plot_object("n",len(species),1,"n"))

    for plotObj in plotObjList:
        for simul in simulList:
            #print "simul:"
            #print simul
            plotObj.plot(simul)

    #plotObjList[0].show_figure()
    #plotObjList[1].show_figure()
    #fig0=plotObjList[0].fig
    #fig1=plotObjList[1].fig
    plt.show(block=True)
    
        
            
