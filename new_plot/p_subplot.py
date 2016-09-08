import numpy
from array_rank import arraylist_rank
from kwarg_default import kwarg_default
from array2D_to_array_list import array2D_to_array_list
from sys import exit
from matplotlib.pyplot import cm
from symmetrize_colormap import invert_cm
from colormap_remove_middle import cm_remove_middle
from diverging_red_blue_colormap import diverging_rb_cm
from get_index_range import get_index_range

class perfect_subplot:
    def __init__(self,data,x,subplot_coordinates,dimensions=None,y=None,**kwargs):
        #this object contains instructions for what to plot in a numpy subplot with gridspec
        
        #dimensions: whether to visualize data as a 2D colormap, or several 1D plots in the same subplot
        #            this can be partially infered from the dimensions of the data
        #data: the data to go into this plot. Should be a 1D, list of 1D, or 2D numpy array
        #x,y: x/y grid points corresponding to the indices in the data array. y value needed for 2D arrays
        #subplot_coordinates: where to put this plot on the gridspec
        #**kwargs: extra visualization options that will be different depending on "dimensions".
        #data_<X>axis: how to intepret the axes of the data array. <X>=x,(line),(y)
        #x_<X>axis: how to intepret the axes of the x array. <X>=x,(line),(y)
        #show_<X>axis: whether to show axis when plotting. <X>=x,y,(z)
        #show_<X>axis_label: whether to show axis label when plotting. <X>=x,y,(z)

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
                exit(1)


        #kwargs related to tagging this subplot
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
            
        self.show_xaxis_ticklabel=kwarg_default("show_xaxis_ticklabel",False,**kwargs)
        self.show_yaxis_ticklabel=kwarg_default("show_yaxis_ticklabel",False,**kwargs)

        self.xscale=kwarg_default("xscale",'linear',**kwargs)
        self.yscale=kwarg_default("yscale",'linear',**kwargs)

        self.xlims=kwarg_default("xlims",None,**kwargs)
        self.ylims=kwarg_default("ylims",None,**kwargs)

        self.vlines=kwarg_default("vlines",None,**kwargs)
        self.hlines=kwarg_default("hlines",None,**kwargs)

        self.yaxis_label_x=kwarg_default("yaxis_label_x",-0.15,**kwargs)
        self.yaxis_label_y=kwarg_default("yaxis_label_y",0.5,**kwargs)

        
                
        #############################################################################
        #
        #                               1D
        #
        #############################################################################
        
        if dimensions==1:
            #going to visualize 1D data
            self.dimensions=1
            
            #split 2D array into list of 1D arrays
            self.data=array2D_to_array_list(data,data_otheraxis)
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
                        exit(1)
                    else:
                        #make the x list the same size as the data list
                        self.x=[x]*len(self.data)                       
                elif arraylist_rank(x) == 2:
                    self.x=array2D_to_array_list(x,x_lineaxis)
                else:
                    print "perfect_subplot: 1D: error: x has too high rank to be plotted"
                    exit(1)
            else:
                if all(len(x[i]) != len(self.data[i]) for i in range(len(self.data))):
                    print "perfect_subplot: 1D: error: data and x have different dimensions"
                    exit(1)
                else:
                    self.x=x
  
            #x is now a list of 1D arrays. Check compatiblity with data
            if len(self.x) != len(self.data):
                print "perfect_subplot: 1D: error: data and x list have different lengths"
                exit(1)
            elif not all(len(self.x[i]) == len(self.data[i]) for i in range(len(self.data))):
                print "perfect_subplot: 1D: error: some data and x have different lengths"
                exit(1)

            self.subplot_coordinates=subplot_coordinates

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
            
            if y is None:
                print "perfect_subplot: 2D: error: need y for 2D plot."
                exit(1)
                
            #terminology suitable for 2D case
            data_yaxis=data_otheraxis


            
            if (data_xaxis == 0) and (data_yaxis == 1):
                self.data=data
            elif (data_xaxis == 1) and (data_yaxis == 0):
                #we always want our x data along the first axis
                self.data=numpy.transpose(data)
            else:
                print "perfect_subplot: 2D: error: data has too high rank to be plotted"
                exit(1)

            if (arraylist_rank(x) !=1) or (arraylist_rank(y) != 1):
                print "perfect_subplot: 2D: error: x and y should have rank 1"
            self.x=x
            self.y=y
                
            self.subplot_coordinates=subplot_coordinates

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

            self.cm=kwarg_default("cm",invert_cm(cm_remove_middle(cm.RdBu,cut_radius=0.1)),**kwargs)
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
                ixrange=get_index_range(x,self.xlims,True)
                xs=xs + [x[ixrange]]
            return xs
        if self.dimensions==2:
            ixrange=get_index_range(self.x,self.xlims,True)
            return self.x[ixrange]

    def y_inrange(self):
        if self.dimensions==1:
            print "warning-error: Cannot take yrange of 1D data."
        if self.dimensions==2:
            iyrange=get_index_range(self.y,self.ylims,True)
            return self.y[iyrange]

    def data_inrange(self):
        #print self.dimensions
        if self.dimensions==1:
            #since a subplot can have 1D data of different length
            data=[]
            for x,d in zip(self.x,self.data):
                ixrange=get_index_range(x,self.xlims,True)
                data=data+[d[ixrange]]
        if self.dimensions==2:
            xlims=self.xlims
            ylims=self.ylims
            ixrange=numpy.array(get_index_range(self.x,xlims,True))
            iyrange=numpy.array(get_index_range(self.y,ylims,True))
            #print ixrange
            #print iyrange
            data=self.data[ixrange[:,numpy.newaxis],iyrange]
            #print numpy.shape(data)
            #print numpy.size(ixrange)
            #print numpy.size(iyrange)
        return data

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
