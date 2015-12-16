import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pyplot import cm
from perfect_simulation import perfect_simulation
import scipy.constants
import numpy


sentinel=object()
rc('text', usetex=True)

#object handles look of plots

#... could probably just be a function that takes a figure object and plotnumber


class plot_object(object):

    def __init__(self,plot_func_str=sentinel,numRows=sentinel,numCols=sentinel,title=sentinel,numColors=sentinel):
        #defines what to plot with this object
        if plot_func_str==sentinel:
            plot_func_str="default"
        plot_funcs = {
            "default": self.default_plot_func,
            "normed_particle_flux": self.particle_flux_plot_func,
            "normed_conductive_heat_flux": self.conductive_heat_flux_plot_func,
            "normed_heat_source": self.heat_source_plot_func,
            "normed_particle_source": self.particle_source_plot_func,
            "normed_flow_outboard": self.flow_outboard_plot_func,
            "normed_flow_inboard": self.flow_inboard_plot_func,
            "kPar_outboard": self.kPar_outboard_plot_func,
            "kPar_inboard": self.kPar_inboard_plot_func,
            "T": self.T_plot_func,
            "n": self.n_plot_func,
        }

        self.cm=cm.rainbow
        self.ls=['-',':','--','-.']
        self.lsi=0 #picks a linestyle
        self.ci=0
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
        self.maxPlotNum=0
        self.fig=plt.figure()
        #self.fig.suptitle(self.title)

        self.background="white"

        if numRows != sentinel:
            self.numRows= numRows
        if numCols != sentinel:
            self.numCols= numCols

    @property
    def num_colors(self):
        return self.numColors
    
    @num_colors.setter
    def num_colors(self,nc):
        self.numColors=nc
        print self.numColors
        self.color=iter(self.cm(numpy.linspace(0,1,nc)))

    def default_plot_func(self,simul,same_color=False):
        print "Cannot plot '"+str(simul)+"', no specific plot function specified!"
            
    def plot_xy_legend_species_subplots(self,x,y,species,legend='',xlabel='',ylabel='',ylimBottom0=False,same_color=False):
        #see if species has a subplot, if not, create one for it
        #print x
        #print y
        
        if same_color==False:
            self.lsi=0
            #print self.color
            print "color index:"+str(self.ci)
            self.ci=self.ci+1
            self.c=next(self.color)
        else:
            self.lsi=self.lsi+1
        linestyle=self.ls[self.lsi]

        i=0
        #print species
        self.plots=[]
        for specy in species:
            
            if specy not in self.species_plot_dict.keys():
                self.species_plot_dict[specy]=self.maxPlotNum+1
                self.maxPlotNum=self.maxPlotNum+1
            #print specy
            #print self.species_plot_dict[specy]
            
            self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]);
            self.ax.plot(x, y[:,i], '-',label=legend,ls=linestyle,c=self.c)
            self.ax.set_title(specy)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
            self.ax.autoscale(True)
            if ylimBottom0:
                self.ax.set_ylim(bottom=0)
            i=i+1
        self.ax.legend(bbox_to_anchor=(0.24,0.24,1,1), loc='upper right', borderaxespad=0.,bbox_transform = plt.gcf().transFigure,prop={'size':6})

        
    def particle_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_particle_flux
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$\langle \vec{\Gamma}\cdot \psi \rangle$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,same_color=same_color)

    def particle_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_particle_source
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$S_p/(s^2/m^6)$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,same_color=same_color)
        

    def conductive_heat_flux_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_conductive_heat_flux
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$\langle \vec{q}\cdot \psi \rangle$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,same_color=same_color)

    def heat_source_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_heat_source
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$S_h/(s^2/m^6)$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,same_color=same_color)

    def T_plot_func(self,simul,same_color=False):
        e=scipy.constants.e
        x=simul.psi
        y=simul.T/(1000*e)
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$T/keV$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    def n_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.n
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$n/m^{-3}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    def flow_outboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_flow_outboard
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$V_\parallel/(m/s)$ Outboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    def flow_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.normed_flow_inboard
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$V_\parallel/(m/s)$ Inboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    def kPar_inboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.kPar_inboard
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$k_\parallel$ Inboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    def kPar_outboard_plot_func(self,simul,same_color=False):
        x=simul.psi
        y=simul.kPar_outboard
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$k_\parallel$ Outboard"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True,same_color=same_color)

    @property
    def background(self):
        return self.background_color
    @background.setter
    def background(self,color):
        self.background_color=color
        self.fig.patch.set_facecolor(self.background_color)
    
    def plot(self,simulation,same_color=False):
        #print "func to plot:"
        #print self.plot_func
        self.plot_func(simulation,same_color)

    #def show_figure(self):
    #    self.fig.show()

    def save_figure(self):
        #could probably take some arguments for name or file extension.
        self.fig.savefig(self.title+'.png')

