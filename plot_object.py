import matplotlib.pyplot as plt
import matplotlib
from matplotlib.pyplot import cm
from perfect_simulation import perfect_simulation
import scipy.constants
import numpy
from latexify import latexify
from get_index_range import get_index_range
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

from truncate_colormap import truncate_colormap
from shiftedColorMap import shiftedColorMap



#needed to place the offset text (the 10^x text by the axis) at arbitrary position: https://github.com/matplotlib/matplotlib/issues/4476
import types
import matplotlib.transforms as mtransforms

def monkey_patch(axis, func):
    axis._update_offset_text_position = types.MethodType(func, axis)

def x_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'bottom':
        if bboxes:
            bbox = mtransforms.Bbox.union(bboxes)
        else:
            bbox = self.axes.bbox
        y = bbox.ymin - self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='top', ha='right')

    if self.offset_text_position == 'there':
        # y in axes coords, x i0n display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top+15# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 1.03# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    if self.offset_text_position == 'there2':
        # y in axes coords, x i0n display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top+15# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 0.97# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    else:
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        y = bbox.ymax + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='bottom', ha='left')

    self.offsetText.set_position((x, y))

def y_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'left':
        # y in axes coords, x in display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    if self.offset_text_position == 'there':
        # y in axes coords, x in display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax
        right=self.axes.bbox.xmin
        #print x
        #print y
        y = top-10# - 2*self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        x = 0.03# + 0.1*self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    else:
        # x & y in display coords
        self.offsetText.set_transform(mtransforms.IdentityTransform())

        # Northwest of upper-right corner of right-hand extent of tick labels
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        top, right = bbox.ymax, bbox.xmax
        x = right + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    self.offsetText.set_position((x, y))

sentinel=object()

class plot_object(object):

    def __init__(self,plot_func_str=sentinel,numRows=sentinel,numCols=sentinel,title=sentinel,numColors=sentinel,xlim=sentinel,ylim=sentinel):

        #self.set_style_counter=0
        
        #defines what to plot with this object
        if plot_func_str==sentinel:
            plot_func_str="default"
        plot_funcs = {
            "default": self.default_plot_func,
            "normed_particle_flux": self.normed_particle_flux_plot_func,
            "particle_flux": self.particle_flux_plot_func,
            "normed_particle_flux_over_n": self.normed_particle_flux_over_n_plot_func,
            "particle_flux_over_n": self.particle_flux_over_n_plot_func,
            "normed_momentum_flux": self.normed_momentum_flux_plot_func,
            "momentum_flux": self.momentum_flux_plot_func,
            "normed_momentum_flux_over_n": self.normed_momentum_flux_over_n_plot_func,
            "momentum_flux_over_n": self.momentum_flux_over_n_plot_func,
            
            "normed_conductive_heat_flux": self.normed_conductive_heat_flux_plot_func,
            "conductive_heat_flux": self.conductive_heat_flux_plot_func,
            "normed_conductive_heat_flux_over_n": self.normed_conductive_heat_flux_over_n_plot_func,
            "conductive_heat_flux_over_n": self.conductive_heat_flux_over_n_plot_func,
            "normed_particle_source": self.normed_particle_source_plot_func,
            "particle_source": self.particle_source_plot_func,
            "normed_heat_source": self.normed_heat_source_plot_func,
            "heat_source": self.heat_source_plot_func,
            "normed_ambipolarity": self.normed_ambipolarity_plot_func,
            "ambipolarity": self.ambipolarity_plot_func,
            "normed_flow_outboard": self.normed_flow_outboard_plot_func,
            "flow_outboard": self.flow_outboard_plot_func,
            "normed_flow_inboard": self.normed_flow_inboard_plot_func,
            "flow_inboard": self.flow_inboard_plot_func,
            "kPar_outboard": self.kPar_outboard_plot_func,
            "kPar_inboard": self.kPar_inboard_plot_func,
            "main_main_collisionality": self.main_main_collisionality_plot_func,
            "normed_FSABJPar": self.normed_FSABJPar_plot_func,
            "FSABJPar": self.FSABJPar_plot_func,
            "normed_FSABFlow": self.normed_FSABFlow_plot_func,
            "FSABFlow": self.FSABFlow_plot_func,
            "T": self.T_plot_func,
            "n": self.n_plot_func,
            "eta": self.eta_plot_func,
            "Phi": self.Phi_plot_func,
            "U": self.U_plot_func,
            "deltaN": self.deltaN_plot_func,
            "deltaEta": self.deltaEta_plot_func,
            "deltaT": self.deltaT_plot_func,
            "density_perturbation": self.density_perturbation_plot_func,
            "potential_perturbation": self.potential_perturbation_plot_func,
            "total_density_perturbation": self.total_density_perturbation_plot_func,
            "baseline_inputs": self.baseline_inputs_plot_func,
            "baseline_outputs": self.baseline_outputs_plot_func,
            "baseline_sources": self.baseline_sources_plot_func,
            "baseline_sources_over_conc": self.baseline_sources_over_conc_plot_func,
            "baseline_kPar": self.baseline_kPar_plot_func,
            "baseline_kPar_inboard": self.baseline_kPar_inboard_plot_func,
            "baseline_flows": self.baseline_flows_plot_func,
            "baseline_fluxes": self.baseline_fluxes_plot_func,
            "baseline_all_fluxes": self.baseline_all_fluxes_plot_func,
            "baseline_momentum_flux": self.baseline_momentum_flux_plot_func,
            "n_i_z": self.n_i_z_plot_func,
            "particle_source_i_z": self.particle_source_i_z_plot_func,
            "heat_source_i_z": self.heat_source_i_z_plot_func,
            "baseline_flows_i_z_noJ": self.baseline_flows_i_z_noJ_plot_func,
        }

        self.xlim=xlim
        self.ylim=ylim


        #self.cm=truncate_colormap(cm.nipy_spectral,0,0.95)
        #self.cm=truncate_colormap(cm.gnuplot2,0.1,0.75)
        self.cm=cm.rainbow
        
        self.ls=['solid','dashed','dashdot','dotted']
        self.species_colors=[[r'#08519c',r'#a50f15',r'#006d2c'],[r'#6baed6',r'#fb6a4a',r'#74c476']]
        self.lw=[1.0,0.5]
        
        self.lsi=0 #picks a linestyle
        self.ci=0 #picks a color
        self.lwi=-1
        
        
        self.current_row=-1 #picks a column in xyz plots
        if numColors==sentinel:
            #uses a setter
            self.num_colors=10
        else:
            self.num_colors=numColors

        self.plot_func=plot_funcs[plot_func_str]
        #print self.plot_func
    
        #self.plotNum = 1
        if numRows==sentinel:
            self.numRows = 1
        else:
            self.numRows=numRows
        if numCols==sentinel:
            self.numCols = 1
        else:
            self.numCols=numCols
        if title==sentinel:
            self.title=''
        else:
            self.title=title
        
        self.species_plot_dict={}
        self.species_data_minmax={}
        if self.ylim==sentinel:
            self.manual_scale=True #uses above minmax to set plotted datarange
        else:
            self.manual_scale=False
        self.maxPlotNum=0


        if numRows != sentinel:
            self.numRows= numRows
        if numCols != sentinel:
            self.numCols= numCols

        self.postproc=self.no_post_proc
        self.ylimBottom0=False

        self.set_size()
        #print self.height
        latexify(fig_height=self.height) #call this before plt.figure??
        
        self.fig=plt.figure()
        #self.fig.suptitle(self.title)

        self.background="white"

    @property
    def background(self):
        return self.background_color
    @background.setter
    def background(self,color):
        self.background_color=color
        self.fig.patch.set_facecolor(self.background_color)

    @property
    def num_colors(self):
        return self.numColors
    
    @num_colors.setter
    def num_colors(self,nc):
        self.numColors=nc
        #print self.numColors
        self.color=iter(self.cm(numpy.linspace(0,1,nc)))

    def set_size(self):
        twoD_plots=[self.density_perturbation_plot_func,self.potential_perturbation_plot_func,self.total_density_perturbation_plot_func]
        #sets the margins of the figure to accomodate common labels

        if self.plot_func in twoD_plots:
            self.adjust_bottom=0.15
            self.adjust_left=0.12
            self.adjust_top=0.85
            self.adjust_right=0.99
            base_rows=2.0
        else:
            self.adjust_left=0.17
            self.adjust_bottom=0.15
            self.adjust_top=0.97
            self.adjust_right=0.95
            base_rows=3.0
            
        base_topbot_adjust=self.adjust_bottom+(1-self.adjust_top)
        def row_to_height(rows):
            return (rows/base_rows)*3.39*(1-base_topbot_adjust)*(numpy.sqrt(5)-1.0)/2.0 +3.39*base_topbot_adjust*(numpy.sqrt(5)-1.0)/2.0
        
        self.height=row_to_height(self.numRows)
        self.base_height=row_to_height(base_rows)
        self.adjust_bottom=self.adjust_bottom*self.base_height/self.height
        self.adjust_top=1-(1-self.adjust_top)*self.base_height/self.height

    def update_min_max(self,data,species_plot_index):
        if len(data)==0:
            print "update_min_max: warning, input data is empty."
        newmax=numpy.max(data)
        if newmax>self.species_data_minmax[species_plot_index][1]:
            self.species_data_minmax[species_plot_index][1]=newmax
        newmin=numpy.min(data)
        if newmin<self.species_data_minmax[species_plot_index][0]:
            self.species_data_minmax[species_plot_index][0]=newmin

    def set_y_scale(self):
            for specy in self.species_plot_dict.keys():
                if self.manual_scale==True:

                    if self.ylimBottom0==False:
                        lower=self.species_data_minmax[specy][0]-0.1*numpy.fabs(self.species_data_minmax[specy][0])
                        upper=self.species_data_minmax[specy][1]+0.1*numpy.fabs(self.species_data_minmax[specy][1])
                    else:
                        lower=0
                        upper=self.species_data_minmax[specy][1]+0.1*numpy.fabs(self.species_data_minmax[specy][1])                    
                    
                else:
                    for specy in self.species_plot_dict.keys():
                        if self.ylimBottom0==False:
                            lower=self.ylim[0]-0.1*numpy.fabs(self.ylim[0])
                            upper=self.ylim[1]-0.1*numpy.fabs(self.ylim[1])
                        else:
                            lower=0
                            upper=self.ylim[1]-0.1*numpy.fabs(self.ylim[1])               
                self.fig.axes[self.species_plot_dict[specy]-1].set_ylim(lower,upper)


    def default_plot_func(self,simul,same_color=False):
        print "Cannot plot '"+str(simul)+"', no specific plot function specified!"

    def plot_xy_data_sameplot_species_multiplot(self,x,ys,species,ylabels=None,ylimBottom0s=None,ylogscales=None,mark_zeros=None,share_x=True,same_color=False):
        #should have the same number of species for all the datasets!!

        if same_color==False:
            self.lwi+=1

        if ylabels==None:
            ylabels=['']*len(ys)

        if ylimBottom0s==None:
            ylimBottom0s=[False]*len(ys)

        if ylogscales==None:
            ylogscales=["lin"]*len(ys)

        if mark_zeros==None:
            mark_zeros=[False]*len(ys)
        
        if len(ys[0][0,:]) != self.numRows:
            print "Warning: number of rows in plot_object not equal to number of species to plot!"
        if 1 != self.numCols:
            print "Warning: number of columns in plot_object not equal to 1!"
        
        if len(ylabels) != len(ys):
            print "Warning: number of ylabels not equal to number of ys!"
            
        if len(ylimBottom0s) != len(ys):
            print "Warning: number of ylimBottom0s not equal to number of ys!"

        if len(ylogscales) != len(ys):
            print "Warning: number of ylogscales not equal to number of ys!"

        if len(mark_zeros) != len(ys):
            print "Warning: number of markzeros not equal to number of ys!"

        if share_x==True:
            self.fig.subplots_adjust(hspace=0.01)
            xticklabels=[]
        
        i=0 #data index
        for y in ys:
            for i_s in range(len(ys[i][0,:])):
                if i_s==0:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, 1);
                    first_ax=self.ax
                else:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, i_s+1,sharex=first_ax);




                #self.ax.autoscale(True)
                if self.xlim!=sentinel:
                    self.ax.set_xlim(self.xlim)
                
                monkey_patch(self.ax.yaxis, y_update_offset_text_position)
                self.ax.yaxis.offset_text_position="there"

            
                if same_color:
                    #local
                    linestyle="dashed"
                    linewidth=1
                    marker=''
                    markevery=10
                    markersize=3
                else:
                    #global
                    linestyle="solid"
                    linewidth=1
                    marker=''
                    markevery=10
                    markersize=3
                colors=self.species_colors[i]
                linewidth=self.lw[self.lwi]
                self.ax.plot(x, ys[i][:,i_s], ls=linestyle,c=colors[i_s],marker=marker,markevery=markevery,linewidth=linewidth,markersize=markersize)

                if ylogscales[i_s]=='log':
                    self.ax.set_yscale('log', nonposy='clip',subsy=[10])
                    self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                    for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                        ticklabel.set_visible(False)
                    self.ax.get_yticklabels()[0].set_visible(False)
                    self.ax.get_yticklabels()[-1].set_visible(False)
                    self.ax.get_yticklabels()[1].set_visible(False)
                    self.ax.get_yticklabels()[-2].set_visible(False)
                    #self.ax.get_yticklabels()[2].set_visible(False)
                    #self.ax.get_yticklabels()[-3].set_visible(False)
                elif ylogscales[i_s]=='lin':
                    self.ax.ticklabel_format(axis='y', style='sci', scilimits=(-0,1))
                    self.ax.yaxis.set_major_locator(MaxNLocator(5))
                    self.ax.get_yticklabels()[0].set_visible(False)
                    self.ax.get_yticklabels()[-1].set_visible(False)
                elif ylogscales[i_s]=='symlog':
                    #self.ax.set_yscale('symlog',subsy=[10])
                    self.ax.set_yscale('symlog')
                    #self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                    #for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                    #    ticklabel.set_visible(False)
                    #self.ax.get_yticklabels()[0].set_visible(False)
                    #self.ax.get_yticklabels()[-1].set_visible(False)
                    #self.ax.get_yticklabels()[1].set_visible(False)
                    #self.ax.get_yticklabels()[-2].set_visible(False)
                else:
                    print "Warning: unrecognized log-scale!"

                self.ax.set_ylabel(ylabels[i])
                labelx=-0.15
                self.ax.yaxis.set_label_coords(labelx, 0.5)

                if ylimBottom0s[i_s]:
                    if ylogscales[i_s]:
                        self.ax.set_ylim(bottom=0.00001) #works for baseline n_z
                    else:
                        self.ax.set_ylim(bottom=0)

                else:
                    if mark_zeros[i_s]:
                        self.ax.axhline(y=0,color='k',linestyle=':')

                if share_x==True:
                    if i_s !=len(ys[i][0,:])-1:
                        xticklabels=xticklabels+self.ax.get_xticklabels()
                self.ax.set_title(species[i_s],y=0.40,x=1.04,fontsize=8)
            i=i+1
        if share_x==True:
            plt.setp(xticklabels, visible=False)
        self.postproc=self.xy_species_postproc



    
    def plot_xy_species_multiplot_data_multiplot(self,x,ys,titles,ylabels=None,ylimBottom0s=None,ylogscales=None,mark_zeros=None,share_x=True,same_color=False):

        if same_color==False:
            self.lwi+=1
        
        length=0
        for y in ys:
            length+=len(y[0,:])
        
        if ylabels==None:
            ylabels=['']*length

        if ylimBottom0s==None:
            ylimBottom0s=[False]*length

        if ylogscales==None:
            ylogscales=["lin"]*length

        if mark_zeros==None:
            mark_zeros=[False]*length
        
        if length != self.numRows:
            print "Warning: number of rows in plot_object not equal to number of species + data to plot!"
        if 1 != self.numCols:
            print "Warning: number of columns in plot_object not equal to 1!"
        
        if len(ylabels) != length:
            print "Warning: number of ylabels not equal to number of species-data!"
            
        if len(ylimBottom0s) != length:
            print "Warning: number of ylimBottom0s not equal to number of species-data!"

        if len(ylogscales) != length:
            print "Warning: number of ylogscales not equal to number of species-data!"

        if len(mark_zeros) != length:
            print "Warning: number of markzeros not equal to number of species-data!"

        if share_x==True:
            self.fig.subplots_adjust(hspace=0.01)
            xticklabels=[]
        
        i=0 #data index
        i_t=0 #total index (data + species
        colors=self.species_colors[0]
        for y in ys:
            for i_s in range(len(ys[i][0,:])):
                if i_t==0:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, 1);
                    first_ax=self.ax
                else:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, i_t+1,sharex=first_ax);

                #self.ax.autoscale(True)
                if self.xlim!=sentinel:
                    self.ax.set_xlim(self.xlim)
                
                monkey_patch(self.ax.yaxis, y_update_offset_text_position)
                self.ax.yaxis.offset_text_position="there"
            
                if same_color:
                    #local
                    linestyle="dashed"
                    linewidth=1
                    marker=''
                    markevery=10
                    markersize=3
                else:
                    #global
                    linestyle="solid"
                    linewidth=1
                    marker=''
                    markevery=10
                    markersize=3
                linewidth=self.lw[self.lwi]
                self.ax.plot(x, ys[i][:,i_s], ls=linestyle,c=colors[i_s],marker=marker,markevery=markevery,linewidth=linewidth,markersize=markersize)

                if ylogscales[i_s]=='log':
                    self.ax.set_yscale('log', nonposy='clip',subsy=[10])
                    self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                    for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                        ticklabel.set_visible(False)
                    self.ax.get_yticklabels()[0].set_visible(False)
                    self.ax.get_yticklabels()[-1].set_visible(False)
                    self.ax.get_yticklabels()[1].set_visible(False)
                    self.ax.get_yticklabels()[-2].set_visible(False)
                    #self.ax.get_yticklabels()[2].set_visible(False)
                    #self.ax.get_yticklabels()[-3].set_visible(False)
                elif ylogscales[i_s]=='lin':
                    self.ax.ticklabel_format(axis='y', style='sci', scilimits=(-0,1))
                    self.ax.yaxis.set_major_locator(MaxNLocator(5))
                    self.ax.get_yticklabels()[0].set_visible(False)
                    self.ax.get_yticklabels()[-1].set_visible(False)
                elif ylogscales[i_s]=='symlog':
                    #self.ax.set_yscale('symlog',subsy=[10])
                    self.ax.set_yscale('symlog')
                    #self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                    #for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                    #    ticklabel.set_visible(False)
                    #self.ax.get_yticklabels()[0].set_visible(False)
                    #self.ax.get_yticklabels()[-1].set_visible(False)
                    #self.ax.get_yticklabels()[1].set_visible(False)
                    #self.ax.get_yticklabels()[-2].set_visible(False)
                else:
                    print "Warning: unrecognized log-scale!"

                self.ax.set_ylabel(ylabels[i])
                labelx=-0.15
                self.ax.yaxis.set_label_coords(labelx, 0.5)

                if ylimBottom0s[i_t]:
                    if ylogscales[i_t]:
                        self.ax.set_ylim(bottom=0.00001) #works for baseline n_z
                    else:
                        self.ax.set_ylim(bottom=0)

                else:
                    if mark_zeros[i_t]:
                        self.ax.axhline(y=0,color='k',linestyle=':')

                if share_x==True:
                    if i_t !=length-1:
                        xticklabels=xticklabels+self.ax.get_xticklabels()
                self.ax.set_title(titles[i_t],y=0.40,x=1.04,fontsize=8)
                i_t+=1
            i=i+1
        if share_x==True:
            plt.setp(xticklabels, visible=False)
        self.postproc=self.xy_species_postproc

        
    def plot_xy_species_sameplot_data_multiplot(self,x,ys,species,ylabels=None,ylimBottom0s=None,ylogscales=None,mark_zeros=None,share_x=True,same_color=False):


        if ylabels==None:
            ylabels=['']*len(ys)

        if ylimBottom0s==None:
            ylimBottom0s=[False]*len(ys)

        if ylogscales==None:
            ylogscales=["lin"]*len(ys)

        if mark_zeros==None:
            mark_zeros=[False]*len(ys)
        
        if len(ys) != self.numRows:
            print "Warning: number of rows in plot_object not equal to number of datasets to plot!"
        if 1 != self.numCols:
            print "Warning: number of columns in plot_object not equal to 1!"
        
        if len(ylabels) != len(ys):
            print "Warning: number of ylabels not equal to number of ys!"
            
        if len(ylimBottom0s) != len(ys):
            print "Warning: number of ylimBottom0s not equal to number of ys!"

        if len(ylogscales) != len(ys):
            print "Warning: number of ylogscales not equal to number of ys!"

        if len(mark_zeros) != len(ys):
            print "Warning: number of markzeros not equal to number of ys!"

        if share_x==True:
            self.fig.subplots_adjust(hspace=0.01)
            xticklabels=[]
        
        i=0
        for y in ys:
            if i==0:
                self.ax = self.fig.add_subplot(len(ys), self.numCols, 1);
                first_ax=self.ax
            else:
                self.ax = self.fig.add_subplot(self.numRows, self.numCols, i+1,sharex=first_ax);




            #self.ax.autoscale(True)
            if self.xlim!=sentinel:
                self.ax.set_xlim(self.xlim)
                
            monkey_patch(self.ax.yaxis, y_update_offset_text_position)
            self.ax.yaxis.offset_text_position="there"

            
            
            for i_s in range(len(ys[i][0])):
                if same_color:
                    linestyle=self.ls[i_s]
                    colors=self.species_colors[1]
                    marker=''
                else:
                    linestyle=self.ls[i_s]
                    colors=self.species_colors[0]
                    marker=''
                self.ax.plot(x, ys[i][:,i_s], ls=linestyle,c=colors[i_s],marker=marker)

            if ylogscales[i]=='log':
                self.ax.set_yscale('log', nonposy='clip',subsy=[10])
                self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                    ticklabel.set_visible(False)
                self.ax.get_yticklabels()[0].set_visible(False)
                self.ax.get_yticklabels()[-1].set_visible(False)
                self.ax.get_yticklabels()[1].set_visible(False)
                self.ax.get_yticklabels()[-2].set_visible(False)
                #self.ax.get_yticklabels()[2].set_visible(False)
                #self.ax.get_yticklabels()[-3].set_visible(False)
            elif ylogscales[i]=='lin':
                self.ax.ticklabel_format(axis='y', style='sci', scilimits=(-0,1))
                self.ax.yaxis.set_major_locator(MaxNLocator(5))
                self.ax.get_yticklabels()[0].set_visible(False)
                self.ax.get_yticklabels()[-1].set_visible(False)
            elif ylogscales[i]=='symlog':
                #self.ax.set_yscale('symlog',subsy=[10])
                self.ax.set_yscale('symlog')
                #self.ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
                #for ticklabel in self.ax.get_yticklabels()[0:2:-1]:
                #    ticklabel.set_visible(False)
                #self.ax.get_yticklabels()[0].set_visible(False)
                #self.ax.get_yticklabels()[-1].set_visible(False)
                #self.ax.get_yticklabels()[1].set_visible(False)
                #self.ax.get_yticklabels()[-2].set_visible(False)
            else:
                print "Warning: unrecognized log-scale!"

            self.ax.set_ylabel(ylabels[i])
            labelx=-0.15
            self.ax.yaxis.set_label_coords(labelx, 0.5)
            
            if ylimBottom0s[i] == True:
                if ylogscales[i]:
                    self.ax.set_ylim(bottom=0.00001) #works for baseline n_z
                else:
                    self.ax.set_ylim(bottom=0)
            elif ylimBottom0s[i] != False:
                self.ax.set_ylim(bottom=ylimBottom0s[i])
                if mark_zeros[i]:
                    self.ax.axhline(y=0,color='k',linestyle=':')
                
            if share_x==True:
                if i !=len(ys)-1:
                    xticklabels=xticklabels+self.ax.get_xticklabels()
                    
            i=i+1
        if share_x==True:
            plt.setp(xticklabels, visible=False)
        self.postproc=self.xy_species_postproc
        
        
    def plot_xy_legend_species_subplots(self,x,y,species,legend='',xlabel='',ylabel='',ylimBottom0=False,same_color=False,share_x=True):
        #see if species has a subplot, if not, create one for it
        #print x
        #print y
        
        if same_color==False:
            self.lsi=0
            #print self.color
            #print "color index:"+str(self.ci)
            self.ci=self.ci+1
            self.c=next(self.color)
        else:
            self.lsi=self.lsi+1
        linestyle=self.ls[self.lsi]

        

        
        if share_x==True:
            self.fig.subplots_adjust(hspace=0.01)
            xticklabels=[]

        i=0
        for specy in species:
            
            if specy not in self.species_plot_dict.keys():
                self.species_plot_dict[specy]=self.maxPlotNum+1
                self.maxPlotNum=self.maxPlotNum+1
                self.species_data_minmax[specy]=[float("inf"),float("-inf")]
            #print specy
            #print self.species_plot_dict[specy]

            if share_x==True:
                if i==0:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]);
                    first_ax=self.ax
                else:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy],sharex=first_ax);
            else:
                self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy])
            self.ax.ticklabel_format(axis='y', style='sci', scilimits=(-0,1))
            #self.ax.yaxis.set_offset_position('left') #can only do left or right.
            #self.ax.yaxis.offset_text_position="right" #does not work

            #dynamically modifies object to allow offset_text to be changed
            monkey_patch(self.ax.yaxis, y_update_offset_text_position)
            self.ax.yaxis.offset_text_position="there"
                
            self.ax.set_title(specy,y=0.40,x=1.04,fontsize=40)
            #self.ax.text(s=specy,y=0.1,x=0.1,fontsize=8)


            #self.ax.tick_params(axis='y', which='both', labelleft='off', labelright='on')

            self.ax.plot(x, y[:,i], '-',label=legend,ls=linestyle,c=self.c)


            
            self.ax.set_title(specy)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
            self.ax.yaxis.set_major_locator(MaxNLocator(5))
            self.ax.autoscale(True)
            if ylimBottom0:
                self.ax.set_ylim(bottom=0)
            else:
                self.ax.axhline(y=0,color='k',linestyle=':')
                
            if self.xlim!=sentinel:
                self.ax.set_xlim(self.xlim)
                xlim_indices= get_index_range(x,self.xlim)
            else:
                xlim_indices=[0,len(x)-1]
                
            if self.ylim!=sentinel:
                self.ax.set_ylim(self.ylim)
                ylim_indices= get_index_range(y,self.ylim)
            else:
                ylim_indices=[0,len(y)-1]

            #gets min/max values in xrange
            self.update_min_max(y[xlim_indices[0]:xlim_indices[1],i],specy)
            
            i=i+1
            
            if share_x==True:
                if specy!=species[-1]:
                    xticklabels=xticklabels+self.ax.get_xticklabels()
                #plt.setp(self.ax.get_yticklabels()[1::2], visible=False)
                plt.setp(self.ax.get_yticklabels()[0], visible=False)
                plt.setp(self.ax.get_yticklabels()[-1], visible=False)
        if share_x==True:
            plt.setp(xticklabels, visible=False)
        self.postproc=self.xy_species_postproc
        


    def plot_xy_legend(self,x,y,legend='',xlabel='',ylabel='',ylimBottom0=False,same_color=False):
        #here we do not keep track of species
        self.numRows=1
        self.plot_xy_legend_species_subplots(x,y[:,numpy.newaxis],species=[""],legend='',xlabel=xlabel,ylabel=ylabel,ylimBottom0=ylimBottom0,same_color=same_color)
        
        self.postproc=self.xy_postproc


    def plot_colormap_xyz_species_subplots(self,x,y,z,species,legend='',xlabel='',ylabel='',zlabel='',zlimBottom0=False,same_color=False,share_x=True,share_y=True):

        #print "in plot_colormap_xyz_species"
        
        X,Y=numpy.meshgrid(x,y)
        
        self.current_row=self.current_row+1
        
        if same_color==False:
            self.lsi=0
            #print self.color
            #print "color index:"+str(self.ci)
            self.ci=self.ci+1
            self.c=next(self.color)
        else:
            self.lsi=self.lsi+1
        linestyle=self.ls[self.lsi]


        

        
        if share_x==True:
            self.fig.subplots_adjust(hspace=0.01)
            yticklabels=[]

        if share_y==True:
            self.fig.subplots_adjust(wspace=0.01)
            yticklabels=[]
            
        i=0
        for specy in species:
            if specy not in self.species_plot_dict.keys():
                #here, PlotNum refers to the plots column number
                self.species_plot_dict[specy]=self.maxPlotNum+1
                self.maxPlotNum=self.maxPlotNum+1
                self.species_data_minmax[specy]=[float("inf"),float("-inf")]
            if share_x==True:
                if i==0:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]+self.current_row*self.numCols);
                    first_ax=self.ax
                else:
                    self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]+self.current_row*self.numCols,sharey=first_ax);
                self.ax.ticklabel_format(axis='y', style='sci', scilimits=(-1,2))
                #dynamically modifies object to allow offset_text to be changed
                monkey_patch(self.ax.yaxis, y_update_offset_text_position)
                self.ax.yaxis.offset_text_position="there"
                
                #self.ax.set_title(specy,y=1.14,x=0.5,fontsize=40)
                #self.ax.text(s=specy,y=0.1,x=0.1,fontsize=8)
            else:
                self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]+self.current_row*self.numCols)
            #print self.species_plot_dict[specy]+self.current_row*self.numCols
            self.ax.pcolormesh(X, Y, numpy.transpose(z[:,:,i]),label=legend)
            self.ax.spines['bottom'].set_color(self.c)
            self.ax.spines['top'].set_color(self.c)
            self.ax.spines['left'].set_color(self.c)
            self.ax.spines['right'].set_color(self.c)
            self.ax.spines['bottom'].set_linestyle(linestyle)
            self.ax.spines['top'].set_linestyle(linestyle)
            self.ax.spines['left'].set_linestyle(linestyle)
            self.ax.spines['right'].set_linestyle(linestyle)
 
            
            #,ls=linestyle,c=self.c
            #self.ax.set_title(specy)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
            #self.ax.yaxis.set_major_locator(MaxNLocator(5))
            #self.ax.autoscale(True)
            # unimplemented special treatment for positive data
            #if zlimBottom0:
            #    self.ax.set_zlim(bottom=0)
            #else:
            #    self.ax.axhline(y=0,color='k',linestyle=':')
            
            if self.xlim!=sentinel:
                self.ax.set_xlim(self.xlim)
                xlim_indices= get_index_range(x,self.xlim)
            else:
                xlim_indices=[0,len(x)-1]
                
            if self.ylim!=sentinel:
                self.ax.set_ylim(self.ylim)
                ylim_indices= get_index_range(y,self.ylim)
            else:
                ylim_indices=[0,len(y)-1]
                
            self.update_min_max(z[xlim_indices[0]:xlim_indices[1],ylim_indices[0]:ylim_indices[1],i],specy)
            
            i=i+1
            if share_x==True:
                if specy!=species[0]:
                    yticklabels=yticklabels+self.ax.get_yticklabels()
                #plt.setp(self.ax.get_yticklabels()[1::2], visible=False)
                plt.setp(self.ax.get_yticklabels()[0], visible=False)
                plt.setp(self.ax.get_yticklabels()[-1], visible=False)
        if share_x==True:
            plt.setp(yticklabels, visible=False)
        self.postproc=self.xyz_species_postproc
        self.old_data=z

    def plot_colormap_xyz(self,x,y,z,legend='',xlabel='',ylabel='',zlabel='',zlimBottom0=False,same_color=False,share_x=True,share_y=True):
        self.numCols=1
        self.plot_colormap_xyz_species_subplots(x,y,z[:,:,numpy.newaxis],species=[""],legend=legend,xlabel=xlabel,ylabel=ylabel,zlabel=zlabel,zlimBottom0=zlimBottom0,same_color=same_color,share_x=share_x,share_y=share_y)
        self.postproc=self.xyz_postproc
        
    def normed_particle_flux_plot_func(self,simul,same_color=False):
        
        x=simul.psi
        y=simul.normed_particle_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle \vec{\Gamma}\cdot \nabla \psi \rangle$/(T m^{-1} s^{-1})"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def particle_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.particle_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{\Gamma}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)


    def normed_particle_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_particle_flux/simul.n
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$(\langle \vec{\Gamma}\cdot \nabla \psi \rangle$/n)/(T m^{2} s^{-1})"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def particle_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.particle_flux/simul.nHat
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{\Gamma}/\hat{n}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_particle_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_particle_source
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$S_p/(s^2/m^6)$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def particle_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.particle_source
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{S}_p$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def particle_source_i_z_plot_func(self,simul,same_color=False):
        x=simul.psi
        species=[0,1]
        y=simul.particle_source[:,species]
        species_list=[simul.species[s] for s in species]
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{S}_p$"
        self.plot_xy_legend_species_subplots(x,y,species_list,legend,same_color=same_color)

    def normed_conductive_heat_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_conductive_heat_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle \vec{q}\cdot \nabla \psi \rangle$/(J T m^{-1} s^{-1})"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def conductive_heat_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.conductive_heat_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{q}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_conductive_heat_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_conductive_heat_flux/simul.n
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$(\langle \vec{q}\cdot \nabla\psi \rangle/n )/(J T m^{2} s^{-1})$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def conductive_heat_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.conductive_heat_flux/simul.nHat
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{q}/\hat{n}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_heat_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_heat_source
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$S_h/(s^2/m^6)$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def heat_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.heat_source
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{S}_h$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def heat_source_i_z_plot_func(self,simul,same_color=False):
        x=simul.psi
        species=[0,1]
        y=simul.heat_source[:,species]
        species_list=[simul.species[s] for s in species]
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{S}_h$"
        self.plot_xy_legend_species_subplots(x,y,species_list,legend,same_color=same_color)

    def normed_momentum_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_momentum_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle \vec{\Pi}\cdot \nabla \psi \rangle/(T m^{1} s^{-2})$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def momentum_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.momentum_flux
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{\Pi}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_momentum_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_momentum_flux/simul.n
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$n^{-1}\langle \vec{\Pi}\cdot \nabla \psi \rangle /(T m^{4} s^{-2})$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def momentum_flux_over_n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.momentum_flux/simul.nHat
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{\Pi}/\hat{n}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)


    def normed_ambipolarity_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_ambipolarity
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\sum_a Z_a \hat{\Gamma}_a/(CT m^{2} s^{-1})$"
        self.plot_xy_legend(x,y,legend,same_color=same_color)
        
    def ambipolarity_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.ambipolarity
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\sum_a Z_a \hat{\Gamma}_a$"
        self.plot_xy_legend(x,y,species,legend,same_color=same_color)




    def T_plot_func(self,simul,same_color=False):
        e=scipy.constants.e
        x=simul.psi
        y=simul.T/(1000*e)
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$T/keV$"
        self.ylimBottom0=True
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color,ylimBottom0=self.ylimBottom0)

    def n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.n*10**(-20)
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$n/(10^{20} m^{-3})$"
        self.ylimBottom0=True
        self.plot_xy_legend_species_subplots(x,y,species,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)

    def n_i_z_plot_func(self,simul,same_color=False):
        x=simul.psi
        species=[0,1]
        y=simul.n[:,species]*10**(-20)

        species_list=[simul.species[s] for s in species]
        
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$n/(10^{20} m^{-3})$"
        self.ylimBottom0=True
        self.plot_xy_legend_species_subplots(x,y,species_list,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)


    def eta_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.eta*10**(-20)
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\eta/(10^{20} m^{-3})$"
        self.ylimBottom0=True
        self.plot_xy_legend_species_subplots(x,y,species,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)

    def Phi_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.Phi/1000
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\Phi/kV$"
        self.ylimBottom0=True
        self.plot_xy_legend(x,y,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)

    def normed_FSABJPar_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_FSABJPar
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle B j_\parallel \rangle$/(AT m^{-2})"
        self.ylimBottom0=True
        self.plot_xy_legend(x,y,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)

    def FSABJPar_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.FSABJPar
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle \hat{B} \hat{j}_\parallel \rangle$"
        self.ylimBottom0=True
        self.plot_xy_legend(x,y,legend,ylimBottom0=self.ylimBottom0,same_color=same_color)
        
    def x_x_collisionality_plot_func(self,simul,species_index,same_color=False):
        x=simul.psi
        y=simul.collisionality[:,species_index]
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{\nu}$"
        ylabel2=r"$\nu^*$"
        #,two_axis=[simul.inputs.epsil**(-3.0/2.0),ylabel2]
        self.plot_xy_legend(x,y,legend,same_color=same_color)

        
    def main_main_collisionality_plot_func(self,simul,same_color=False):
        species_index=0
        self.x_x_collisionality_plot_func(simul,species_index,same_color)

    def U_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.U
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$U$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def deltaN_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.deltaN
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\delta_n$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def deltaT_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.deltaT
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\delta_T$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def deltaEta_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.deltaEta
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\delta_{\eta}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_flow_outboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_flow_outboard
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$V_\parallel/(m/s)$ Outboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def flow_outboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.flow_outboard
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{V}_\parallel$ Outboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color) 

    def normed_flow_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_flow_inboard
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$V_\parallel/(m/s)$ Inboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def flow_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.flow_inboard
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\hat{V}_\parallel$ Inboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def normed_FSABFlow_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_FSABFlow
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle B V_\parallel \rangle/(Tm/s)$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def FSABFlow_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.FSABFlow
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\langle \hat{B} \hat{V}_\parallel \rangle$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def kPar_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.kPar_inboard
        
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$k_\parallel$ Inboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def kPar_outboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.kPar_outboard
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$k_\parallel$ Outboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,same_color=same_color)

    def density_perturbation_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.theta/numpy.pi
        y=numpy.append(y,[2])
        z=simul.density_perturbation
        z=numpy.append(z,numpy.expand_dims(z[:,0,:],axis=1),axis=1)
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$100\,\psi_N$"
        self.ylabel=r"$\theta/\pi$"
        self.zlabel=r""
        self.plot_colormap_xyz_species_subplots(x,y,z,species,legend,same_color=same_color)

    def total_density_perturbation_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.theta/numpy.pi
        y=numpy.append(y,[2])
        z=simul.total_density_perturbation
        z=numpy.append(z,numpy.expand_dims(z[:,0,:],axis=1),axis=1)
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$100\,\psi_N$"
        self.ylabel=r"$\theta/\pi$"
        self.zlabel=r""
        self.plot_colormap_xyz_species_subplots(x,y,z,species,legend,same_color=same_color)

    def potential_perturbation_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.theta/numpy.pi
        y=numpy.append(y,[2])
        z=simul.potential_perturbation
        z=numpy.append(z,numpy.expand_dims(z[:,0],axis=1),axis=1)
        species=simul.species
        legend=simul.description
        self.share_y_label=True
        self.share_x_label=True
        self.xlabel=r"$100\,\psi_N$"
        self.ylabel=r"$\theta/\pi$"
        self.zlabel=r""
        self.plot_colormap_xyz(x,y,z,species,legend,same_color=same_color)

    def baseline_inputs_plot_func(self,simul,same_color=False):
        #NOTE: same color is ignored since local and global assumed to be the same.
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False

        main_index=0
        imp_n=numpy.append(0*simul.n[:,[1]],simul.n[:,[1]]*10**(-20),axis=1)
        ys=[simul.n*10**(-20),imp_n,simul.T/(1000*scipy.constants.e),numpy.expand_dims(simul.Phi/1000,axis=1),simul.U,simul.deltaN,simul.deltaT,numpy.expand_dims(simul.collisionality[:,main_index],axis=1)]
        ylabels=[r"$n/(10^{20} m^{-3})$",r"",r"$T/keV$",r"$\Phi/kV$",r"$U$",r"$\delta_n$",r"$\delta_T$",r"$\hat{\nu}$"]
        ylimBottom0s=[True,True,True,-0.1,True,True,True,False]
        ylogscales=['lin','lin','lin','lin','log','log','log','lin']
        mark_zeros=[False,False,False,False,False,False,False,False]
        self.plot_xy_species_sameplot_data_multiplot(x,ys,simul.species,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=False)

    def baseline_outputs_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False

        main_index=0
        ys=[simul.particle_flux,simul.conductive_heat_flux,simul.particle_source,simul.heat_source,simul.kPar_inboard,simul.kPar_outboard]
        ylabels=[r"$\hat{\Gamma}/\hat{n}$",r"$\hat{q}/\hat{n}$",r"$\hat{S}_p$",r"$\hat{S}_h$",r"$k_\parallel$ I",r"$k_\parallel$ O"]
        ylimBottom0s=[False,False,False,False,False,False]
        ylogscales=['lin','lin','lin','lin','lin','lin']
        #ylogscales=['symlog','symlog','symlog','symlog','lin','lin']
        mark_zeros=[True,True,True,True,True,True]
        self.plot_xy_species_sameplot_data_multiplot(x,ys,simul.species,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    
    def baseline_flows_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False
        ys=[numpy.expand_dims(simul.FSABJPar,axis=1),simul.FSABFlow]
        ylabels=[r"$\langle \hat{B} \hat{j}_\parallel \rangle$",r"$\langle \hat{B} \hat{V}_\parallel \rangle$"]
        ylimBottom0s=[False,False,False,False]
        titles=['']+simul.species
        ylogscales=['lin','lin','lin','lin']
        #ylogscales=['symlog','symlog','symlog','symlog','lin','lin']
        mark_zeros=[False,False,False,False]
        self.plot_xy_species_multiplot_data_multiplot(x,ys,titles,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_flows_i_z_noJ_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False
        species=[0,1]
        ys=[simul.FSABFlow[:,species]]
        ylabels=[r"$\langle \hat{B} \hat{V}_\parallel \rangle$"]
        ylimBottom0s=[False,False]
        titles=[simul.species[s] for s in species]
        ylogscales=['lin','lin']
        mark_zeros=[False,False]
        self.plot_xy_species_multiplot_data_multiplot(x,ys,titles,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

        
    def baseline_sources_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$S_p, S_h$"
        self.share_x_label=True
        self.share_y_label=True
        species=[0,1]
        species_list=[simul.species[s] for s in species]
        ys=[simul.particle_source[:,species],simul.heat_source[:,species]]
        ylabels=[r"",r""]
        ylimBottom0s=[False,False]
        ylogscales=['lin','lin']
        mark_zeros=[True,True]
        self.plot_xy_data_sameplot_species_multiplot(x,ys,species_list,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_sources_over_conc_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$\{S_p,\, S_h\} (n_i/n_a)_{\mathrm{core}}$"
        self.share_x_label=True
        self.share_y_label=True
        species=[0,1]
        species_list=[simul.species[s] for s in species]

        conc=simul.nHat[0,[species]]/simul.nHat[0,0]
        print conc
        
        ys=[simul.particle_source[:,species]/conc,simul.heat_source[:,species]/conc]
        
        ylabels=[r"",r""]
        ylimBottom0s=[False,False]
        ylogscales=['lin','lin']
        mark_zeros=[True,True]
        self.plot_xy_data_sameplot_species_multiplot(x,ys,species_list,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_kPar_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$k_\parallel$"
        self.share_x_label=True
        self.share_y_label=True
        ys=[simul.kPar_inboard,simul.kPar_outboard]
        ylabels=[r"",r"",r""]
        ylimBottom0s=[False,False,False]
        ylogscales=['lin','lin','lin']
        mark_zeros=[True,True,True]
        self.plot_xy_data_sameplot_species_multiplot(x,ys,simul.species,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_kPar_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r"$k_\parallel$"
        self.share_x_label=True
        self.share_y_label=True
        ys=[simul.kPar_inboard]
        ylabels=[r"",r"",r""]
        ylimBottom0s=[False,False,False]
        ylogscales=['lin','lin','lin']
        mark_zeros=[True,True,True]
        self.plot_xy_data_sameplot_species_multiplot(x,ys,simul.species,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_fluxes_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False
        titles=simul.species+simul.species
        ys=[simul.particle_flux/simul.nHat,simul.conductive_heat_flux/simul.nHat]
        ylabels=[r"$\hat{\Gamma}/\hat{n}$",r"$\hat{q}/\hat{n}$"]
        ylimBottom0s=[False,False,False,False,False,False]
        ylogscales=['lin','lin','lin','lin','lin','lin']
        mark_zeros=[True,True,True,True,True,True]
        self.plot_xy_species_multiplot_data_multiplot(x,ys,titles,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color)

    def baseline_momentum_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False
        titles=simul.species
        ys=[simul.momentum_flux/simul.nHat]
        ylabels=[r"$\hat{\Pi}/\hat{n}$"]
        ylimBottom0s=[False,False,False]
        ylogscales=['lin','lin','lin']
        mark_zeros=[True,True,True]
        self.plot_xy_species_multiplot_data_multiplot(x,ys,titles,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color) 

    def baseline_all_fluxes_plot_func(self,simul,same_color=False):
        x=simul.psi
        self.xlabel=r"$\psi_N$"
        self.ylabel=r""
        self.share_x_label=True
        self.share_y_label=False
        titles=simul.species+simul.species+simul.species
        ys=[simul.particle_flux/simul.nHat,simul.conductive_heat_flux/simul.nHat,simul.momentum_flux/simul.nHat]
        ylabels=[r"$\hat{\Gamma}/\hat{n}$",r"$\hat{q}/\hat{n}$",r"$\hat{\Pi}/\hat{n}$"]
        ylimBottom0s=[False,False,False,False,False,False,False,False,False]
        ylogscales=['lin','lin','lin','lin','lin','lin','lin','lin','lin']
        mark_zeros=[True,True,True,True,True,True,True,True,True]
        self.plot_xy_species_multiplot_data_multiplot(x,ys,titles,ylabels=ylabels,ylimBottom0s=ylimBottom0s,ylogscales=ylogscales,mark_zeros=mark_zeros,same_color=same_color) 
    
        
    def plot(self,simulation,same_color=False):
        #print "func to plot:"
        #print self.plot_func
        self.plot_func(simulation,same_color=same_color)

    def no_post_proc(self):
        #for when no post processing of the figure is needed
        pass


    def xy_postproc(self):
        if self.share_y_label==True:
            self.fig.text(0.01, 0.5+self.adjust_bottom/4, self.ylabel, va='center', rotation='vertical')

        if self.share_x_label==True:
            self.fig.text(0.5+self.adjust_left/4, 0.01, self.xlabel, ha='center')

        self.fig.subplots_adjust(left=self.adjust_left)
        self.fig.subplots_adjust(top=self.adjust_top)
        self.fig.subplots_adjust(bottom=self.adjust_bottom)
        self.fig.subplots_adjust(right=self.adjust_right)

        self.set_y_scale()

    def xy_species_postproc(self):
        self.xy_postproc()
        #if self.share_y_label==True:
        #    self.fig.text(0.01, 0.5, self.ylabel, va='center', rotation='vertical')
        #    adjust_left=0.17
        #    self.fig.subplots_adjust(left=0.17)
        #else:
        #    adjust_left=0

        #if self.share_x_label==True:
        #    self.fig.subplots_adjust(bottom=0.15)
        #    self.fig.subplots_adjust(top=0.97)
        #    self.fig.subplots_adjust(right=0.95)
        #    self.fig.text(0.5+adjust_left/4, 0.01, self.xlabel, ha='center')
        #self.set_y_scale()


    def xyz_species_postproc(self):
        if self.share_y_label==True:
            self.fig.text(0.01, 0.5+self.adjust_bottom/4.0, self.ylabel, va='center', rotation='vertical')

        if self.share_x_label==True:
            self.fig.text(0.5+self.adjust_left/4, 0.01*self.base_height/self.height, self.xlabel, ha='center')
        self.fig.subplots_adjust(bottom=self.adjust_bottom)
        self.fig.subplots_adjust(top=self.adjust_top)
        self.fig.subplots_adjust(right=self.adjust_right)
        self.fig.subplots_adjust(left=self.adjust_left)

        if self.manual_scale==True:
            for specy in self.species_plot_dict.keys():
                i_s=self.species_plot_dict[specy]
                for j in range(self.numRows):
                    if self.ylimBottom0:
                        vmin=0
                    else:
                        vmin=self.species_data_minmax[specy][0]
                        vmax=self.species_data_minmax[specy][1]
                    self.fig.axes[self.numCols*j+i_s-1].collections[0].set_clim(vmin,vmax)
                    #shift cmap to have 0 in the middle
                    shift=1 - vmax/(vmax + abs(vmin))
                    print "colormap shift: " + str(shift)
                    my_cm=shiftedColorMap(cm.coolwarm,midpoint=shift,name="mymap")
                    self.fig.axes[self.numCols*j+i_s-1].collections[0].set_clim(vmin,vmax)
                    self.fig.axes[self.numCols*j+i_s-1].collections[0].set_cmap(my_cm)
                    
                    if j == self.numRows-1:
                        #self.fig.axes[self.numCols*j+i_s-1].ticklabel_format(axis='x', style='sci', scilimits=(0,0))
                        scale_pow = 2
                        def my_formatter_fun(x, p):
                            return "%.1f"%(x * (10 ** scale_pow))

                        self.fig.axes[self.numCols*j+i_s-1].get_xaxis().set_major_formatter(ticker.FuncFormatter(my_formatter_fun))
                        self.fig.axes[self.numCols*j+i_s-1].xaxis.set_major_locator(MaxNLocator(4))
                        plt.setp(self.fig.axes[self.numCols*j+i_s-1].get_xticklabels()[-1], visible=False)
                        plt.setp(self.fig.axes[self.numCols*j+i_s-1].get_xticklabels()[0], visible=False)
                        #self.fig.axes[self.numCols*j+i_s-1].ticklabel_format(axis='x', style='sci', scilimits=(0,0))
                        #self.fig.axes[self.numCols*j+i_s-1].xaxis.set_powerlimits((0, 0))

                        #self.fig.axes[self.numCols*j+i_s-1].update_ticks()
                        
                        #cb=self.fig.colorbar(self.fig.axes[self.numCols*j+i_s-1].collections[0],ax=self.fig.axes[self.numCols*j+i_s-1],orientation="horizontal",use_gridspec=False,anchor=[0.5,1.0])#,anchor=[0.5,0.0]
                    elif j == 0:
                        self.fig.axes[self.numCols*j+i_s-1].xaxis.set_major_locator(MaxNLocator(4))

                        self.fig.axes[self.numCols*j+i_s-1].set_xticklabels([])
                        #self.fig.axes[self.numCols*j+i_s-1]

                        #to make coloraxis derive position from yaxis
                        divider = make_axes_locatable(self.fig.axes[self.numCols*j+i_s-1])
                        cax = divider.append_axes("top", "-10%", pad="10%")
                        cb=self.fig.colorbar(self.fig.axes[self.numCols*j+i_s-1].collections[0],cax=cax,orientation="horizontal")
                        tick_locator = ticker.MaxNLocator(nbins=4)
                        cb.locator=tick_locator
                        cb.ax.xaxis.set_ticks_position('top')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.tick_params(direction='out', pad=1)
                        #cb.format=ticker.FuncFormatter(fmt)
                        cb.formatter.set_powerlimits((0, 0))
                        cb.update_ticks()
                        #cb.ax.xaxis.get_offset_text().set_position((1,-0.5)) #this only works for x position! ;-;
                        self.fig.axes[self.numCols*j+i_s-1].set_title(specy,y=1.19,x=0.5,fontsize=8)
                        #cb.ax.xaxis.set_offset_position("top")
                        #cb.ax.text(-0.25, 1, r'$\times$10$^{-1}$', va='bottom', ha='left')
                        monkey_patch(cb.ax.xaxis, x_update_offset_text_position)
                        cb.ax.xaxis.offset_text_position="there2"
                        #ax=self.fig.axes[self.numCols*j+i_s-1]
                    else:
                        self.fig.axes[self.numCols*j+i_s-1].set_xticklabels([])
                        self.fig.axes[self.numCols*j+i_s-1].xaxis.set_major_locator(MaxNLocator(4))
                        self.fig.axes[self.numCols*j+i_s-1].set_title("",y=1.12,x=0.5,fontsize=8)


                            
                            
    def xyz_postproc(self):
        self.xyz_species_postproc()

    def save_figure(self):
        self.postproc()
        #print self.fig.axes
        self.fig.savefig(self.title+'.pdf')

