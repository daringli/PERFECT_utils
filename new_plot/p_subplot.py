import numpy
from array_rank import arraylist_rank
from kwarg_default import kwarg_default
from array2D_to_array_list import array2D_to_array_list
from sys import exit
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from symmetrize_colormap import invert_cm
from colormap_remove_middle import cm_remove_middle
from diverging_red_blue_colormap import diverging_rb_cm
from get_index_range import get_index_range

from mpl_toolkits.axes_grid1 import make_axes_locatable
from monkey_patch_axes import monkey_patch
from matplotlib import ticker
from symmetrize_colormap import symmetrize_colormap
from streamplot import streamplot

import cmocean # has some nice perceptually uniform diverging colormaps
from radar_marker_plot import radar_marker_plot

import matplotlib


class perfect_subplot:
    def __init__(self,data,x,subplot_coordinates,dimensions=None,y=None,x_data=None, y_data=None,radarvalues=None,**kwargs):
        #this object contains instructions for what to plot in a numpy subplot with gridspec
        
        #dimensions: whether to visualize data as a 2D colormap (2), or several 1D plots in the same subplot (1)
        #            this can be partially infered from the dimensions of the data
        #data: the data to go into this plot. Should be a 1D, list of 1D, or 2D numpy array
        # x_data,y_data: same as data, but x and y component for stream plots
        #x,y: x/y grid points corresponding to the indices in the data array. y value needed for 2D arrays
        #subplot_coordinates: where to put this plot on the gridspec
        #**kwargs: extra visualization options that will be different depending on "dimensions".
        #data_<X>axis: how to intepret the axes of the data array. <X>=x,(line),(y)
        #x_<X>axis: how to intepret the axes of the x array. <X>=x,(line),(y)
        #show_<X>axis: whether to show axis when plotting. <X>=x,y,(z)
        #show_<X>axis_label: whether to show axis label when plotting. <X>=x,y,(z)
        #period_x: None if x coordinate is not periodic, otherwise period
        #period_y: None if y coordinate is not periodic, otherwise period
        

        
        #kwargs that may change how we intepret the data
        data_xaxis=kwarg_default("data_xaxis",0,**kwargs)
        data_otheraxis=(data_xaxis+1)%2

        #use rank of input data to try to determine dimensions
        #and reshape appropriately
        rank=arraylist_rank(data)
        if rank == 1:
            #put it in a list for consistency with multi-array code
            data=[data]
            if dimensions==None:
                dimensions=1
        elif rank == 2:
            #can be a 2D function or several 1D.
            #please manually specify to avoid ambiguity
            #default assumes second index reperesents species/simulation
            #if it is under or equal to 4.#
            if dimensions==None:
                print "perfect_subplot: warning: ambigious array rank."
                if y == None:
                    dimensions=1
                else:
                    dimensions=2
        elif type(rank) is list:
            if all(x == 1 for x in rank):
                if dimensions==None:
                    dimensions=1
            else:
                print rank
                print "perfect_subplot: error: input list has array(s) of rank above 1."
                raise ValueError
            

        #Kwargs related to tagging this subplot
        self.groups=kwarg_default("groups",[],**kwargs)
                
        #kwargs related to flagging how the subplot should be visualized
        self.border_color=kwarg_default("border_color",'k',**kwargs)
        self.border_linestyle=kwarg_default("border_linestyle",'solid',**kwargs)
        self.title=kwarg_default("title","",**kwargs)
        
        #axis propert kwargs
        self.show_xaxis=kwarg_default("show_xaxis",True,**kwargs)
        self.show_yaxis=kwarg_default("show_yaxis",True,**kwargs)

        #powerlimit should be a tuple: (-a,b)
        self.xaxis_powerlimits=kwarg_default("xaxis_powerlimits",None,**kwargs)
        self.yaxis_powerlimits=kwarg_default("yaxis_powerlimits",None,**kwargs)

        
        self.xaxis_label=kwarg_default("xaxis_label",None,**kwargs)
        self.yaxis_label=kwarg_default("yaxis_label",None,**kwargs)

        self.xaxis_label_color=kwarg_default("xaxis_label_color",'Black',**kwargs)
        self.yaxis_label_color=kwarg_default("yaxis_label_color",'Black',**kwargs)

        self.xaxis_label_size=kwarg_default("xaxis_label_size",12,**kwargs)
        self.yaxis_label_size=kwarg_default("yaxis_label_size",12,**kwargs)
        print self.xaxis_label_size
        
        self.show_xaxis_label=kwarg_default("show_xaxis_label",True,**kwargs)
        self.show_yaxis_label=kwarg_default("show_yaxis_label",True,**kwargs)
        
        self.show_xaxis_ticklabel=kwarg_default("show_xaxis_ticklabel",True,**kwargs)
        self.show_yaxis_ticklabel=kwarg_default("show_yaxis_ticklabel",True,**kwargs)

        self.xscale=kwarg_default("xscale",'linear',**kwargs)
        self.yscale=kwarg_default("yscale",'linear',**kwargs)

        self.xlims=kwarg_default("xlims",None,**kwargs)
        self.ylims=kwarg_default("ylims",None,**kwargs)

        self.vlines=kwarg_default("vlines",[],**kwargs)
        self.hlines=kwarg_default("hlines",[],**kwargs)
        self.vlines_colors=kwarg_default("vlines_colors","k",**kwargs)
        self.hlines_colors=kwarg_default("hlines_colors","silver",**kwargs)

        self.yaxis_label_x=kwarg_default("yaxis_label_x",-0.15,**kwargs)
        self.yaxis_label_y=kwarg_default("yaxis_label_y",0.5,**kwargs)
        self.xaxis_label_x=kwarg_default("xaxis_label_x",0.5,**kwargs)
        self.xaxis_label_y=kwarg_default("xaxis_label_y",-0.35,**kwargs)

        self.datastyle=kwarg_default("datastyle","color",**kwargs)
        
        self.bg_alpha=kwarg_default("bg_alpha",1.0,**kwargs)

        self.subplot_coordinates=subplot_coordinates

                
        #############################################################################
        #
        #                               1D
        #
        #############################################################################
        
        if dimensions==1:
            self.radarvalues = radarvalues
            #going to visualize 1D data
            self.dimensions=1

            #default plot-type is simply 1D plot (not really used for 1D plot)
            self.plot_type=kwarg_default("plot_type","line",**kwargs) 
            
            #split 2D array into list of 1D arrays
            self.data=array2D_to_array_list(data,data_otheraxis)
            #print data
            #print  self.data

            #terminology suitable for 1D case
            data_lineaxis=data_otheraxis

            #kwarg needed to intepret x, which can be a list of arrays here
            x_xaxis=kwarg_default("x_xaxis",0,**kwargs)
            x_lineaxis=(x_xaxis+1)%2

            if type(arraylist_rank(x)) is not list:
                if arraylist_rank(x) == 1:
                    if all(len(x) != len(self.data[i]) for i in range(len(self.data))):
                        print "perfect_subplot: 1D: error: data and x have different dimensions"
                        raise ValueError
                    else:
                        #make the x list the same size as the data list
                        self.x=[x]*len(self.data)                       
                elif arraylist_rank(x) == 2:
                    self.x=array2D_to_array_list(x,x_lineaxis)
                else:
                    print "perfect_subplot: 1D: error: x has too high rank to be plotted"
                    raise ValueError
            else:
                if all(len(x[i]) != len(self.data[i]) for i in range(len(self.data))):
                    print "perfect_subplot: 1D: error: data and x have different dimensions"
                    raise ValueError
                else:
                    self.x=x
  
            #x is now a list of 1D arrays. Check compatiblity with data
            if len(self.x) != len(self.data):
                print "perfect_subplot: 1D: error: data and x list have different lengths"
                raise ValueError
            elif not all(len(self.x[i]) == len(self.data[i]) for i in range(len(self.data))):
                print [self.x[i].shape for i in range(len(self.x))]
                print [self.data[i].shape for i in range(len(self.data))]
                print "perfect_subplot: 1D: error: some data and x have different lengths"
                raise ValueError
                
            self.x_period=kwarg_default("x_period",None,**kwargs)
            #aesthetic kwargs
            self.xticks=kwarg_default("xticks",6,**kwargs)
            if self.yscale=='linear':
                self.yticks=kwarg_default("yticks",5,**kwargs)
            if self.yscale=='log':
                self.yticks=kwarg_default("yticks",4,**kwargs)
                
            self.hidden_xticklabels=kwarg_default("hidden_xticklabels",[],**kwargs)
            self.hidden_yticklabels=kwarg_default("hidden_yticklabels",[0,-1],**kwargs)
            self.linestyles=kwarg_default("linestyles",["solid"]*len(x),**kwargs)
            self.markers=kwarg_default("markers",[""]*len(x),**kwargs)
            self.markeverys=kwarg_default("markeverys",[1]*len(x),**kwargs)
            self.markersizes=kwarg_default("markersizes",[3]*len(x),**kwargs)
            self.fillstyles=kwarg_default("fillstyles",["full"]*len(x),**kwargs)
            self.linewidths=kwarg_default("linewidths",[1]*len(x),**kwargs)
            self.colors=kwarg_default("colors",["k"]*len(x),**kwargs)


            #title position
            self.title_y=kwarg_default("title_y",0.4,**kwargs)
            self.title_x=kwarg_default("title_x",1.04,**kwargs)
            

        #############################################################################
        #
        #                               2D
        #
        #############################################################################
        
        if dimensions==2:
            #going to visualize 2D data
            self.dimensions=2

            #z powerlimits
            self.zaxis_powerlimits=kwarg_default("zaxis_powerlimits",(0,0),**kwargs)
            
            self.linecolor=kwarg_default("linecolor","k",**kwargs)
            #terminology suitable for 2D case
            data_yaxis=data_otheraxis

            if (arraylist_rank(x) !=1) or (arraylist_rank(y) != 1):
                print "perfect_subplot: 2D: error: x and y should have rank 1"
            self.x=x
            self.y=y
            self.x_period=kwarg_default("x_period",None,**kwargs)
            self.y_period=kwarg_default("y_period",2*numpy.pi,**kwargs)
            
            
            if (data_xaxis == 0) and (data_yaxis == 1):
                self.data=data
            elif (data_xaxis == 1) and (data_yaxis == 0):
                #we always want our x data along the first axis
                self.data=numpy.transpose(data)
            else:
                print "perfect_subplot: 2D: error: data has too high rank to be plotted"
                raise ValueError

            
            if y is None:
                print "perfect_subplot: 2D: error: need y for 2D plot."
                raise ValueError
            if (x_data!=None) and (y_data!=None):
                self.xy_data_exists=True
                #default plot-type is streamplot
                self.plot_type=kwarg_default("plot_type","streamplot",**kwargs)
                self.arrow_density = kwarg_default("arrow_density",(2,3),**kwargs)


                if (data_xaxis == 0) and (data_yaxis == 1):
                    self.y_data = y_data
                    self.x_data = x_data
                elif (data_xaxis == 1) and (data_yaxis == 0):
                    #we always want our x data along the first axis
                    self.y_data = numpy.transpose(y_data)
                    self.x_data = numpy.transpose(x_data)
                    
            elif (x_data==None and y_data!=None) or (x_data!=None and y_data==None):
                print "perfect_subplot: 2D: error: need x_data and y_data for vector plot."
                raise ValueError
            else:
                self.xy_data_exists=False
                #default plot-type is pcolormesh 
                self.plot_type=kwarg_default("plot_type","pcolormesh",**kwargs) 
                

            
            #kwargs related to flagging how the subplot should be visualized
            #z-axis refers to the colorbar for the colormap
            self.show_zaxis=kwarg_default("show_zaxis",True,**kwargs)
            
            self.zaxis_label=kwarg_default("zaxis_label",None,**kwargs)
            
            self.show_zaxis_ticklabel=kwarg_default("show_zaxis_ticklabel",False,**kwargs)

            self.zscale=kwarg_default("zscale",'linear',**kwargs)

            self.xticks=kwarg_default("xticks",4,**kwargs)
            self.yticks=kwarg_default("yticks",4,**kwargs)
            self.zticks=kwarg_default("zticks",4,**kwargs)

            self.zlims=kwarg_default("zlims",None,**kwargs)

            #self.cm=kwarg_default("cm",invert_cm(cm_remove_middle(cm.RdBu,cut_radius=0.1)),**kwargs)
            self.cm=kwarg_default("cm",cm_remove_middle(cmocean.cm.delta,cut_radius=0.05),**kwargs) #GB
            #self.cm=kwarg_default("cm",cm_remove_middle(cmocean.cm.balance,cut_radius=0.05),**kwargs) #RB
            #self.cm=kwarg_default("cm",cmocean.cm.balance,**kwargs)
            #self.cm=kwarg_default("cm",invert_cm(cm.RdBu),**kwargs)
            #self.cm=kwarg_default("cm",diverging_rb_cm(),**kwargs)
            self.symmetrize_cm=kwarg_default("symmetrize_cm",True,**kwargs)

            self.hidden_xticklabels=kwarg_default("hidden_xticklabels",[0,-1],**kwargs)
            self.hidden_yticklabels=kwarg_default("hidden_yticklabels",[],**kwargs)
            self.hidden_zticklabels=kwarg_default("hidden_zticklabels",[],**kwargs)

            #title position
            self.title_y=kwarg_default("title_y",1.19,**kwargs)
            self.title_x=kwarg_default("title_x",0.5,**kwargs)

        

            

    #functions data in the range that will be plotted

    def x_inrange(self):
        if self.dimensions==1:
            xs=[]
            for x in self.x:
                ixrange=get_index_range(x,self.xlims,True,self.x_period)
                xs=xs + [x[ixrange]]
            return xs
        if self.dimensions==2:
            ixrange=get_index_range(self.x,self.xlims,True,self.x_period)
            return self.x[ixrange]

    def y_inrange(self):
        if self.dimensions==1:
            print "warning-error: Cannot take yrange of 1D data."
        if self.dimensions==2:
            iyrange=get_index_range(self.y,self.ylims,True,self.y_period)
            return self.y[iyrange]

    def data_inrange(self):
        #print self.dimensions
        if self.dimensions==1:
            #since a subplot can have 1D data of different length
            data=[]
            for x,d in zip(self.x,self.data):
                ixrange=get_index_range(x,self.xlims,True,self.x_period)
                data=data+[d[ixrange]]
            return data
        if self.dimensions==2:
            xlims=self.xlims
            ylims=self.ylims
            ixrange=numpy.array(get_index_range(self.x,xlims,True,self.x_period))
            iyrange=numpy.array(get_index_range(self.y,ylims,True,self.y_period))
            data=self.data[ixrange[:,numpy.newaxis],iyrange]
            return data

    def plot(self,fig,ax):
        if self.dimensions == 1:
            x=self.x
            y=self.data
            if self.plot_type == "csv_all":
                f = open(str(ax)+".perfect_subplot",'w')
                #export this subplot to a csv file
                l=[]
                for (_x,_y) in zip(x,y):
                    for (__x,__y) in zip(_x,_y):
                        l.append(str(__x)+","+str(__y))
                f.write("\n".join(l))
                f.close()
                return None
            
            ax.yaxis.offset_text_position="there"
            
            if type(self.linestyles) is not list:
                self.linestyles=[self.linestyles]*len(y)
                self.linewidths=[self.linewidths]*len(y)
            for i in range(len(y)):
                if self.plot_type == "line":
                    try:
                        #print str(i) + ": " + str(self.linestyles[i])
                        ax.plot(x[i],y[i],linestyle=self.linestyles[i],color=self.colors[i],linewidth=self.linewidths[i],marker=self.markers[i],markevery=self.markeverys[i],fillstyle=self.fillstyles[i],markersize=self.markersizes[i],)
                    except IndexError:
                        print "Index error in ax.plot(). Most likely, linestyles, linewidths and colors have the wrong lengths."
                        print "len(colors): " + str(len(self.colors))
                        print "len(linestyles): " + str(len(self.linestyles))
                        print "len(linewidths): " + str(len(self.linewidths))
                        print "len(markers): " + str(len(self.markers))
                        print "len(markersizes: " + str(len(self.markersizes))
                        print "len(markeverys): " + str(len(self.markeverys))
                        print "len(fillstyles): " + str(len(self.fillstyles))
                        print "len(x): " + str(len(x))
                        print "len(y): " + str(len(y))
                elif self.plot_type == "radarmarker":
                    try:
                        radar_marker_plot(ax,x[i],y[i],values=self.radarvalues[i],colors=self.colors[i],radius=self.markersizes[i],markeredgewidth=0.1)
                    except IndexError:
                        print "Index error radar_marker_plot(). Most likely, linestyles, linewidths and colors have the wrong lengths."
                        print "len(colors): " + str(len(self.colors))
                        print "len(markersizes): " + str(len(self.markersizes))
                        print "len(radarvalues): " + str(len(self.radarvalues))


                elif self.plot_type == "csv":
                    #export this dataset in this subplot to a csv file
                    f = open(str(i)+".perfect_subplot",'w')
                    for (_x,_y) in zip(x[i],y[i]):
                        f.write(str(_x)+","+str(_y)+"\n")
                    f.close()

                
                    
                ax.yaxis.set_label_coords(-0.15, 0.5)

        ############ 2D #################
            
        if self.dimensions == 2:
            if self.zlims is not None:
                #modify colormap
                zmin=self.zlims[0]
                zmax=self.zlims[1]
            else:
                zmin=numpy.min(self.data)
                zmax=numpy.max(self.data)

            if self.symmetrize_cm:
                if zmax*zmin<=0:
                    if abs(zmax)>abs(zmin):
                        vmin = -zmax
                        vmax = zmax
                    else:
                        vmin = zmin
                        vmax = -zmin
                elif zmin>0:
                    vmin=-zmax
                    vmax=zmax
                else:
                    vmin=zmin
                    vmax=-zmin
            else:
                vmin=zmin
                vmax=zmax

            if self.xy_data_exists:                
                #streamplot
                z=numpy.transpose(self.data)
                U=numpy.transpose(self.x_data)
                V=numpy.transpose(self.y_data)
                if self.datastyle == "color":
                    color=z
                    linecolor=None
                elif (self.datastyle == "") or (self.datastyle is None):
                    color=None
                    #linecolor=self.linecolor
                    linecolor=None
                if self.plot_type == "streamplot":
                    streamplot(ax,self.x,self.y,U,V,density=self.arrow_density,cmap=self.cm,color=color,vmin=vmin,vmax=vmax,linecolor=linecolor,rasterized=True)
                    #rasterizes collection of lines and arrows in streamplot as one
                else:
                    print "perfect_subplot: 2D: error: unrecognized plot-type for vector data"
                    raise ValueError
            else:
                #pcolor
                X,Y=numpy.meshgrid(self.x,self.y)
                z=numpy.transpose(self.data)
                if self.plot_type == "pcolormesh":
                    ax.pcolormesh(X, Y, z,rasterized=True,linewidth=0,cmap=self.cm,vmin=vmin,vmax=vmax)
                elif self.plot_type == "contour":
                    levels = [0.0]
                    ax.contour(X, Y, z,colors='k',rasterized=True,levels=levels)
                else:
                    print "perfect_subplot: 2D: error: unrecognized plot-type for 2D data"
                    raise ValueError
                
            if self.show_zaxis == True:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("top", "-10%", pad="10%")
                cb=fig.colorbar(ax.collections[0],cax=cax,orientation="horizontal")
                #cb.solids.set_rasterized(True)
                #cb.solids.set_edgecolor("face")
                monkey_patch(cb.ax.xaxis, 'x')
                if self.yscale == 'linear':
                    tick_locator = ticker.MaxNLocator(nbins=self.zticks)
                if self.yscale == 'log':
                    tick_locator = ticker.LogLocator(nbins=self.zticks)
                cb.locator=tick_locator
                cb.ax.xaxis.set_ticks_position('top')
                cb.ax.xaxis.set_label_position('top')
                cb.ax.set_xscale(self.zscale)
                cb.ax.tick_params(direction='out', pad=1)
                cb.formatter.set_powerlimits(self.zaxis_powerlimits)
                cb.update_ticks()
                if self.show_zaxis_ticklabel == False:
                    plt.setp(cb.ax.get_xticklabels(), visible=False)
                    #makes it impossible to have an offset without a ticklabel
                    cb.formatter.set_powerlimits((float("-inf"), float("inf")))
                    cb.update_ticks()
                cb.ax.xaxis.offset_text_position="there2"
                cb.solids.set_edgecolor("face")
                


            if self.zaxis_label is not None:
                cb.ax.set_xlabel(self.zaxis_label)
            if self.zlims is not None:
                #modify colorbar
                if self.symmetrize_cm:
                    cm=symmetrize_colormap(self.cm,zmin,zmax)
                else:
                    cm=self.cm
                ax.collections[0].set_cmap(cm)
                ax.collections[0].set_clim(zmin,zmax)
                for i in self.hidden_zticklabels:
                    plt.setp(cb.ax.get_xticklabels()[i],visible=False)
                if self.show_zaxis == True:
                    cb.solids.set_rasterized(True)
        ax.patch.set_alpha(self.bg_alpha)

                
                

        

if __name__=="__main__":
    data=numpy.array([[1,2,3],[4,5,6]])
    x=numpy.array([0,1])
    y=numpy.array([-1,-2,-3])
    psp=perfect_subplot(data,x,(0,0))
    print repr(psp.data)
    print repr(psp.x)

    psp=perfect_subplot(data,x,(0,1),y=y)
    print repr(psp.data)
    print repr(psp.x)
    print repr(psp.y)
