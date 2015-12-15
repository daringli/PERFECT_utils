from plot_object import plot_object
from perfect_simulation import perfect_simulation
from perfect_simulations_from_dirs import perfect_simulations_from_dirs
import sys
import matplotlib.pyplot as plt



if __name__=='__main__':
    dirlist=sys.argv[1:]

    if len(dirlist)==0:
        print "No directories specified!"
        sys.exit("Nothing to plot")

    for i in range(len(dirlist)):
        dir=dirlist[2*i]
        try:
            localdir=dir.rsplit('/',1)[0]+"-local"+"/"+dir.rsplit('/',1)[1]
            print localdir
        except IndexError:
            "local dir does not exist. Please specify what to plot as dirname/subdir where dirname-local/subdir exists!"
        else:
            dirlist.insert(2*i+1,localdir)


    print dirlist
    normlist="norms.namelist"
    species="d,N,e"
    species=species.split(',')
    #print species
    simulList=perfect_simulations_from_dirs(dirlist,normlist,species)
    #print simulList

    plotObjList=[]
    ncolors=len(simulList)/2
    plotObjList.append(plot_object("normed_particle_flux",len(species),1,"particle flux",ncolors))
    plotObjList.append(plot_object("normed_conductive_heat_flux",len(species),1,"heat flux",ncolors))
    plotObjList.append(plot_object("T",len(species),1,"T",ncolors))
    plotObjList.append(plot_object("n",len(species),1,"n",ncolors))
    plotObjList.append(plot_object("normed_heat_source",len(species),1,"normed heat source",ncolors))
    plotObjList.append(plot_object("normed_particle_source",len(species),1,"normed particle source",ncolors))

    tf=[False,True]
    i=0

    desc=["$Z_{eff}-1=0.02$","$Z_{eff}-1=0.02$ local","$Z_{eff}-1=0.1$","$Z_{eff}-1=0.1$ local"]

    for plotObj in plotObjList:
        for simul in simulList:
            print i
            simul.description=desc[i%4]
            #print "simul:"
            #print simul
            #print "to call plotObj.plot()"
            plotObj.plot(simul,tf[i%2])
            i=i+1 #not python3 safe
        plotObj.save_figure()

    #plotObjList[0].show_figure()
    #plotObjList[1].show_figure()
    #fig0=plotObjList[0].fig
    #fig1=plotObjList[1].fig
    #plt.show(block=True)
    
        
            
