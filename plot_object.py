import matplotlib.pyplot as plt
from matplotlib import rc
from perfect_simulation import perfect_simulation
import scipy.constants

sentinel=object()
rc('text', usetex=True)

#object handles look of plots

#... could probably just be a function that takes a figure object and plotnumber


class plot_object(object):

    def __init__(self,plot_func_str=sentinel,numRows=sentinel,numCols=sentinel,title=sentinel):
        #defines what to plot with this object
        if plot_func_str==sentinel:
            plot_func_str="default"
        plot_funcs = {
            "default": self.default_plot_func,
            "normed_particle_flux": self.particle_flux_plot_func,
            "normed_conductive_heat_flux": self.conductive_heat_flux_plot_func,
            "T": self.T_plot_func,
            "n": self.n_plot_func,
        }
        
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
        self.fig.suptitle(self.title)

        self.background="white"

        if numRows != sentinel:
            self.numRows= numRows
        if numCols != sentinel:
            self.numCols= numCols

    def default_plot_func(self,simul):
        print "Cannot plot '"+str(simul)+"', no specific plot function specified!"
            
    def plot_xy_legend_species_subplots(self,x,y,species,legend='',xlabel='',ylabel='',ylimBottom0=False):
        #see if species has a subplot, if not, create one for it
        #print x
        #print y
        i=0
        #print species
        for specy in species:
            
            if specy not in self.species_plot_dict.keys():
                self.species_plot_dict[specy]=self.maxPlotNum+1
                self.maxPlotNum=self.maxPlotNum+1
            #print specy
            #print self.species_plot_dict[specy]
            self.ax = self.fig.add_subplot(self.numRows, self.numCols, self.species_plot_dict[specy]);
            self.ax.plot(x, y[:,i], '-',label=legend)
            self.ax.set_title(specy)
            self.ax.legend(loc=1,prop={'size':6})
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
            if ylimBottom0:
                self.ax.set_ylim(bottom=0)
            i=i+1
        
    def particle_flux_plot_func(self,simul):
        x=simul.psi
        y=simul.normed_particle_flux
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$\langle \vec{\Gamma}\cdot \psi \rangle$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel)
        

    def conductive_heat_flux_plot_func(self,simul):
        x=simul.psi
        y=simul.normed_conductive_heat_flux
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$\langle \vec{q}\cdot \psi \rangle$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel)

    def T_plot_func(self,simul):
        e=scipy.constants.e
        x=simul.psi
        y=simul.T/(1000*e)
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$T/eV$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True)

    def n_plot_func(self,simul):
        x=simul.psi
        y=simul.n
        legend=simul.description
        species=simul.species
        xlabel=r"$\psi_N$"
        ylabel=r"$n/m^{-3}$"
        self.plot_xy_legend_species_subplots(x,y,species,legend,xlabel,ylabel,True)


    @property
    def background(self):
        return self.background_color
    @background.setter
    def background(self,color):
        self.background_color=color
        self.fig.patch.set_facecolor(self.background_color)
    
    def plot(self,simulation):
        #print "func to plot:"
        #print self.plot_func
        self.plot_func(simulation)

    def show_figure(self):
        #self.fig.show()
        self.fig.show()

