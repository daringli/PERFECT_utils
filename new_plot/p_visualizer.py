from p_subplot import perfect_subplot
from kwarg_default import kwarg_default
from p_subplot_group import perfect_subplot_group
from latexify import latexify
from monkey_patch_axes import monkey_patch
from symmetrize_colormap import symmetrize_colormap

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker


# for the __name__="__main__" part only:
from perfect_simulations_from_dirs import perfect_simulations_from_dirs


def perfect_visualizer(p_subplot_list,gridspec_params,**kwargs):
    #p_subplot_list: list of perfect_subplot
    #gridspec_params = [rows,cols] of the subplot grid. Should be compatible with the coordinates of the perfect_subplots
    #kwargs: contains instructions to visualize everything

    
    dimensions=kwarg_default("dimensions",1,**kwargs)
    height=kwarg_default("height",0.15,**kwargs)
    textcols=kwarg_default("textcols",1,**kwargs)
    rows=kwarg_default("rows",gridspec_params[0],**kwargs)
    cols=kwarg_default("cols",gridspec_params[1],**kwargs)

    
    #set the height differently depending on if most plots are colormaps (2D), and thus need colorbars
    if dimensions == 1:
        adjust_bottom=kwarg_default("adjust_bottom",0.15,**kwargs)
        adjust_left=kwarg_default("adjust_left",0.17,**kwargs)
        adjust_top=kwarg_default("adjust_top",0.97,**kwargs)
        adjust_right=kwarg_default("adjust_right",0.95,**kwargs)
        base_rows=3.0
    elif dimensions == 2:
        adjust_bottom=kwarg_default("adjust_bottom",0.15,**kwargs)
        adjust_left=kwarg_default("adjust_left",0.12,**kwargs)
        adjust_top=kwarg_default("adjust_top",0.85,**kwargs)
        adjust_right=kwarg_default("adjust_right",0.99,**kwargs)
        base_rows=2.0
    base_topbot_adjust=adjust_bottom+(1-adjust_top)

    def row_to_height(rows):
            return (rows/base_rows)*3.39*(1-base_topbot_adjust)*(numpy.sqrt(5)-1.0)/2.0 +3.39*base_topbot_adjust*(numpy.sqrt(5)-1.0)/2.0


    height=row_to_height(rows)
    base_height=row_to_height(base_rows)
    adjust_bottom=adjust_bottom*base_height/height
    adjust_top=1-(1-adjust_top)*base_height/height
    
    latexify(fig_height=height,columns=textcols)
    fig=plt.figure()


    fig.subplots_adjust(bottom=adjust_bottom)
    fig.subplots_adjust(left=adjust_left)
    fig.subplots_adjust(top=adjust_top)
    fig.subplots_adjust(right=adjust_right)

    adjust_hspace=kwarg_default("adjust_hspace",0.01,**kwargs)
    adjust_wspace=kwarg_default("adjust_wspace",0.01,**kwargs)
    fig.subplots_adjust(hspace=adjust_hspace)
    fig.subplots_adjust(wspace=adjust_wspace)


    global_xlabel=kwarg_default("global_xlabel",None,**kwargs)
    if global_xlabel is not None:
        if dimensions == 1:
            fig.text(0.5+adjust_left/4, 0.01, global_xlabel, ha='center')
        elif dimensions == 2:
            fig.text(0.5+adjust_left/4, 0.01*base_height/height, global_xlabel, ha='center')
    
    global_ylabel=kwarg_default("global_ylabel",None,**kwargs)
    if global_ylabel is not None:
        if dimensions == 1:
            fig.text(0.01, 0.5+adjust_bottom/4, global_ylabel, va='center', rotation='vertical')
        elif dimensions == 2:
            fig.text(0.01, 0.5+adjust_bottom/4.0, global_ylabel, va='center', rotation='vertical')

    
    gs=gridspec.GridSpec(gridspec_params[0],gridspec_params[1])

    ax_list=[]
    
    for p_subplot in p_subplot_list:
        ax_list.append(plt.subplot(gs[p_subplot.subplot_coordinates]))

    for ax,p_subplot in zip(ax_list,p_subplot_list):
        monkey_patch(ax.xaxis,'x')
        monkey_patch(ax.yaxis,'y')

        #keywords in the p_subplot to control looks:
        ax.spines['bottom'].set_color(p_subplot.border_color)
        ax.spines['top'].set_color(p_subplot.border_color)
        ax.spines['left'].set_color(p_subplot.border_color)
        ax.spines['right'].set_color(p_subplot.border_color)

        ax.spines['bottom'].set_linestyle(p_subplot.border_linestyle)
        ax.spines['top'].set_linestyle(p_subplot.border_linestyle)
        ax.spines['left'].set_linestyle(p_subplot.border_linestyle)
        ax.spines['right'].set_linestyle(p_subplot.border_linestyle)
        if p_subplot.vlines is not None:
            for p in p_subplot.vlines:
                ax.axvline(x=p,color='k',linestyle=':')
        if p_subplot.hlines is not None:
            for p in p_subplot.hlines:
                ax.axhline(y=p,color='silver',linestyle=':')
        
        if p_subplot.show_xaxis_ticklabel == False:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.xaxis.get_major_formatter().set_powerlimits((float("-inf"), float("inf")))
        else:
            plt.setp(ax.get_xticklabels(), visible=True)
            
        if p_subplot.show_yaxis_ticklabel == False:
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.yaxis.get_major_formatter().set_powerlimits((float("-inf"), float("inf")))

        if p_subplot.show_xaxis == False:
            plt.setp(ax.xaxis, visible=False)
        if p_subplot.show_yaxis == False:
            plt.setp(ax.yaxis, visible=False)

                
        if p_subplot.xscale == 'linear':
            x_locator=ticker.MaxNLocator(p_subplot.xticks)
        if p_subplot.xscale == 'log':
            x_locator=ticker.LogLocator(p_subplot.xticks)
        ax.xaxis.set_major_locator(x_locator)
            
        if p_subplot.yscale == 'linear':
            y_locator=ticker.MaxNLocator(p_subplot.yticks)
        if p_subplot.yscale == 'log':
            y_locator=ticker.LogLocator(p_subplot.yticks)
        ax.yaxis.set_major_locator(y_locator)

       
        ############## 1D #########
        if p_subplot.dimensions == 1:
            x=p_subplot.x
            y=p_subplot.data            
            ax.yaxis.offset_text_position="there"
            

            
            if type(p_subplot.linestyles) is not list:
                p_subplot.linestyles=[p_subplot.linestyles]*len(y)
                p_subplot.linewidths=[p_subplot.linewidths]*len(y)
            for i in range(len(y)):
                try:
                    ax.plot(x[i],y[i],linestyle=p_subplot.linestyles[i],color=p_subplot.colors[i],linewidth=p_subplot.linewidths[i],marker=p_subplot.markers[i],fillstyle=p_subplot.fillstyles[i])
                except IndexError:
                    print "Index error in ax.plot(). Most likely, linestyles, linewidths and colors have the wrong lengths."
                    print "len(colors): " + str(len(p_subplot.colors))
                    print "len(linestyles): " + str(len(p_subplot.linestyles))
                    print "len(linewidths): " + str(len(p_subplot.linewidths))
                    print "len(markers): " + str(len(p_subplot.markers))
                    print "len(fillstyles): " + str(len(p_subplot.fillstyles))
                    print "len(x): " + str(len(x))
                    print "len(y): " + str(len(y))
            ax.yaxis.set_label_coords(-0.15, 0.5)            

        ############# 2D ###########
        if p_subplot.dimensions == 2:
            X,Y=numpy.meshgrid(p_subplot.x,p_subplot.y)
            z=numpy.transpose(p_subplot.data)
  
            ax.pcolor(X, Y, z,rasterized=True,linewidth=0)
            
            if p_subplot.show_zaxis == True:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("top", "-10%", pad="10%")
                cb=fig.colorbar(ax.collections[0],cax=cax,orientation="horizontal")
                #cb.solids.set_edgecolor("face")
                monkey_patch(cb.ax.xaxis, 'x')
                if p_subplot.yscale == 'linear':
                    tick_locator = ticker.MaxNLocator(nbins=p_subplot.zticks)
                if p_subplot.yscale == 'log':
                    tick_locator = ticker.LogLocator(nbins=p_subplot.zticks)
                cb.locator=tick_locator
                cb.ax.xaxis.set_ticks_position('top')
                cb.ax.xaxis.set_label_position('top')
                cb.ax.set_xscale(p_subplot.zscale)
                cb.ax.tick_params(direction='out', pad=1)
                cb.formatter.set_powerlimits((0, 0))
                cb.update_ticks()
                if p_subplot.show_zaxis_ticklabel == False:
                    plt.setp(cb.ax.get_xticklabels(), visible=False)
                    #makes it impossible to have an offset without a ticklabel
                    cb.formatter.set_powerlimits((float("-inf"), float("inf")))
                    cb.update_ticks()
                cb.ax.xaxis.offset_text_position="there2"

            if p_subplot.zaxis_label is not None:
                cb.ax.set_xlabel(p_subplot.zaxis_label)
            if p_subplot.zlims is not None:
                zmin=p_subplot.zlims[0]
                zmax=p_subplot.zlims[1]
                ax.collections[0].set_clim(zmin,zmax)
                cm=p_subplot.cm
                if p_subplot.symmetrize_cm:
                    cm=symmetrize_colormap(cm,zmin,zmax)
                ax.collections[0].set_cmap(cm)
                for i in p_subplot.hidden_zticklabels:
                    plt.setp(cb.ax.get_xticklabels()[i],visible=False)
            cb.solids.set_edgecolor("face")

        ### QUIVER PLOT #####
        if p_subplot.dimensions == -2:
            X,Y=numpy.meshgrid(p_subplot.x,p_subplot.y)
            U=numpy.transpose(p_subplot.x_data)
            V=numpy.transpose(p_subplot.y_data)
            ax.quiver(X,Y,U,V)


        #we need to sort out ticks after plotting, since they are generated
        #during the plot operation
        if p_subplot.xaxis_label is not None:
            ax.set_xlabel(p_subplot.xaxis_label)
        if p_subplot.yaxis_label is not None:
            #print str(p_subplot.yaxis_label)+": " +"("+str(p_subplot.yaxis_label_x)+","+str(p_subplot.yaxis_label_y)+")"
            ax.set_ylabel(p_subplot.yaxis_label)
            ax.yaxis.set_label_coords(p_subplot.yaxis_label_x,p_subplot.yaxis_label_y) 
        
        if p_subplot.xlims is not None:
            ax.set_xlim(p_subplot.xlims)
        if p_subplot.xlims is not None:
            ax.set_ylim(p_subplot.ylims)

        if p_subplot.xaxis_powerlimits is not None:
            print "Warning: xlabel power limits not fully implemented: offset position unknown, offset not shown."
            ax.ticklabel_format(axis='x', style='sci', scilimits=p_subplot.xaxis_powerlimits)
            #ax.xaxis.get_major_formatter().set_powerlimits(p_subplot.xaxis_powerlimits)
        
        
        if p_subplot.yscale == 'log':
            ax.set_yscale(p_subplot.yscale, nonposy='clip',subsy=[10])
            for ticklabel in ax.get_yticklabels()[0:2:-1]:
                ticklabel.set_visible(False)
            ax.get_yticklabels()[0].set_visible(False)
            ax.get_yticklabels()[-1].set_visible(False)
            ax.get_yticklabels()[1].set_visible(False)
            ax.get_yticklabels()[-2].set_visible(False)
        elif p_subplot.yscale=='linear':
            if p_subplot.yaxis_powerlimits is None:
                ax.ticklabel_format(axis='y', style='sci', scilimits=(-0,1))
            else:
                ax.ticklabel_format(axis='y', style='sci', scilimits=p_subplot.yaxis_powerlimits)
            ax.get_yticklabels()[0].set_visible(False)
            ax.get_yticklabels()[-1].set_visible(False)
        elif p_subplot.yscale=='symlog':
            ax.set_yscale(p_subplot.yscale)
        for i in p_subplot.hidden_xticklabels:
            plt.setp(ax.get_xticklabels()[i],visible=False)
        for i in p_subplot.hidden_yticklabels:  
            plt.setp(ax.get_yticklabels()[i],visible=False)
        

        
        #print str(p_subplot.title)+": " +"("+str(p_subplot.title_x)+","+str(p_subplot.title_y)+")"
        ax.set_title(p_subplot.title,x=p_subplot.title_x,y=p_subplot.title_y,fontsize=8)



        
if __name__=="__main__":


    base_dir='../../../../../../../perfectSimuls/he_paper/nonuniform/const_Phi-m2tanh/'
    sub_dirs=['0.01/','0.01-local']
    dirlist=[base_dir+x for x in sub_dirs]
    normlist=[dir+'/norms.namelist' for dir in dirlist]
    species=[dir+'/species' for dir in dirlist]
    psiAHats=[dir+'/psiAHat.h5' for dir in dirlist]
    simulList=perfect_simulations_from_dirs(dirlist,normlist,species,psiAHats)
    
    psp1=perfect_subplot(simulList[0].flow[:,:,0],x=100*simulList[0].psi,y=simulList[0].theta/numpy.pi,subplot_coordinates=(0,0),show_zaxis_ticklabel=True,show_yaxis_ticklabel=True,zaxis_label="test z",yaxis_label="test y",xaxis_label="test bad x",zscale='linear',xscale='linear',groups=["i","sim1"])
    psp2=perfect_subplot(simulList[0].flow[:,:,1],x=100*simulList[0].psi,y=simulList[0].theta/numpy.pi,subplot_coordinates=(0,1),show_zaxis_ticklabel=True,groups=["z","sim1"])
    psp3=perfect_subplot(simulList[0].flow[:,:,2],x=100*simulList[0].psi,y=simulList[0].theta/numpy.pi,subplot_coordinates=(0,2),show_zaxis_ticklabel=True,groups=["e","sim1"])
    psp4=perfect_subplot(simulList[1].flow[:,:,0],x=100*simulList[1].psi,y=simulList[1].theta/numpy.pi,subplot_coordinates=(1,0),show_xaxis_ticklabel=True,show_yaxis_ticklabel=True,xaxis_label="test x",show_zaxis=False,groups=["i","sim2"])
    psp5=perfect_subplot(simulList[1].flow[:,:,1],x=100*simulList[1].psi,y=simulList[1].theta/numpy.pi,subplot_coordinates=(1,1),show_xaxis_ticklabel=True,show_zaxis=False,groups=["z","sim2"])
    psp6=perfect_subplot(simulList[1].flow[:,:,2],x=100*simulList[1].psi,y=simulList[1].theta/numpy.pi,subplot_coordinates=(1,2),show_xaxis_ticklabel=True,show_zaxis=False,groups=["e","sim2"])

    gridspec_params=[2,3]
    psp_list=[psp1,psp2,psp3,psp4,psp5,psp6]

    igroup=perfect_subplot_group(psp_list,groups=["i"])
    zgroup=perfect_subplot_group(psp_list,groups=["z"])
    egroup=perfect_subplot_group(psp_list,groups=["e"])

    sim1group=perfect_subplot_group(psp_list,groups=["sim1"])
    sim2group=perfect_subplot_group(psp_list,groups=["sim2"])

    #igroup.setattrs("xlims",[igroup.get_min("x"),igroup.get_max("x")])
    #zgroup.setattrs("xlims",[zgroup.get_min("x"),zgroup.get_max("x")])
    #egroup.setattrs("xlims",[egroup.get_min("x"),egroup.get_max("x")])
    igroup.setattrs("xlims",[90,100])
    zgroup.setattrs("xlims",[90,100])
    egroup.setattrs("xlims",[90,100])

    igroup.setattrs("ylims",[igroup.get_min("y"),igroup.get_max("y")])
    zgroup.setattrs("ylims",[zgroup.get_min("y"),zgroup.get_max("y")])
    egroup.setattrs("ylims",[egroup.get_min("y"),egroup.get_max("y")])

    igroup.setattrs("zlims",[igroup.get_min("data"),igroup.get_max("data")])
    zgroup.setattrs("zlims",[zgroup.get_min("data"),zgroup.get_max("data")])
    egroup.setattrs("zlims",[egroup.get_min("data"),egroup.get_max("data")])

    sim1group.setattrs("border_color",'magenta')
    sim1group.setattrs("border_linestyle",'dashed')

    perfect_visualizer(psp_list,gridspec_params,global_xlabel="X",dimensions=2,global_ylabel="Y")

    plt.savefig('test2D.pdf')

    data1=[simulList[i].particle_flux[:,0] for i in range(len(simulList))]
    data2=[simulList[i].particle_flux[:,1] for i in range(len(simulList))]
    data3=[simulList[i].particle_flux[:,2] for i in range(len(simulList))]
    
    x=[simulList[i].psi for i in range(len(simulList))]

    
    psp1=perfect_subplot(data1,x=100*simulList[0].psi,subplot_coordinates=(0,0),show_yaxis_ticklabel=True,yaxis_label="test y",groups=["i"],linestyles=['solid','dashed'])
    psp2=perfect_subplot(data2,x=100*simulList[0].psi,subplot_coordinates=(1,0),groups=["z"],linestyles=['solid','dashed'])
    psp3=perfect_subplot(data3,x=100*simulList[0].psi,subplot_coordinates=(2,0),groups=["e"],linestyles=['solid','dashed'],show_xaxis_ticklabel=True)

    gridspec_params=[3,1]
    psp_list=[psp1,psp2,psp3]

    igroup=perfect_subplot_group(psp_list,groups=["i"])
    zgroup=perfect_subplot_group(psp_list,groups=["z"])
    egroup=perfect_subplot_group(psp_list,groups=["e"])

    igroup.setattrs("xlims",[90,100])
    zgroup.setattrs("xlims",[90,100])
    egroup.setattrs("xlims",[90,100])

    igroup.setattrs("ylims",[igroup.get_min("data"),igroup.get_max("data")])
    zgroup.setattrs("ylims",[zgroup.get_min("data"),zgroup.get_max("data")])
    egroup.setattrs("ylims",[egroup.get_min("data"),egroup.get_max("data")])

    perfect_visualizer(psp_list,gridspec_params,global_xlabel="X",dimensions=1,global_ylabel="Y")

    plt.savefig('test1D.pdf')
