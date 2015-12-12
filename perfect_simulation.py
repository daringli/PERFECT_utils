import sys
from perfectInputFile import perfectInput
import numpy
import f90nml
import h5py
import shutil #copy to copy files
import os #to make directory
import os.path
import errno #to handle exception when making directory
import copy #to copy objects rather than get references to them
#for normed simulation
import scipy.constants
import subprocess
from subprocess import call




sentinel=object() #to detect default inputs

class perfect_simulation(object):

    def __init__(self,input_filename,output_filename=sentinel,group_name=sentinel):
        
        #description of simulation, for usage as legend in plots, etc.
        self.description=""
        self.species=[] #list of species identifiers

        #getting input
        self.input=input_filename
        #self.input_filename=input_filename
        #self.inputs=perfectInput(self.input_filename)
        #self.input_dir = "." + "/" + self.input_filename.rsplit('/',1)[:-1][0]

        #having no output is acceptable,
        #it just means that simulation has yet to be run
        self.output=output_filename

        #Output intepretation and group specification
        #Different group in the same output can in principle
        #be handled by different simulation objects
        self.version_group="/"
        self.program_mode_group="/"
        if group_name!=sentinel:
            self.group_name=group_name
        else:
            self.group_name="/run  1/"

        #hiding implementation details
        self.psi_name="psi"
        self.particle_flux_name="particleFlux"
        self.heat_flux_name="heatFlux"
        self.VPrimeHat_name="VPrimeHat"
        self.Delta_name="Delta"
        self.omega_name="omega"
        self.THat_name="THat"
        self.nHat_name="nHat"
        self.num_species_name="Nspecies"
        self.heat_source_name="heatSourceProfile"
        self.particle_source_name="particleSourceProfile"


    @property
    def input(self):
        return self.input_filename

    @input.setter
    def input(self,input_filename):
        self.input_filename="." + "/" + input_filename
        self.inputs=perfectInput(self.input_filename)
        self.input_dir = self.input_filename.rsplit('/',1)[:-1][0]
        self.input_filename_nodir = "." + "/" + self.input_filename.rsplit('/',1)[-1]
        if self.inputs.profilesScheme==7:
            self.input_profiles_filename=self.input_dir+"/"+self.inputs.profilesFilename
            self.input_profiles_groupname="/Npsi"+str(self.inputs.Npsi)+"/"
            #we only assign self.input_profiles if they exist
            if os.path.isfile(self.input_profiles_filename):
                self.input_profiles=h5py.File(self.input_profiles_filename,'r')

    @property
    def output(self):
        return self.output_filename

    @output.setter
    def output(self,output_filename):
        if (output_filename==sentinel) or (output_filename==None):
            self.output_filename=None
            self.outputs=None
            self.output_dir=self.input_dir
            self.output_filename_nodir =None
        else:
            self.output_filename="." + "/" + output_filename
            try:
                 self.outputs=h5py.File(self.output_filename)
            except:
                print "Could not read desired output file '"+self.output_filename+"', setting it to 'None'."
                self.output_filename=None
                self.outputs=None
                self.output_dir=self.input_dir
                self.output_filename_nodir =None
            else:
                self.output_dir = self.output_filename.rsplit('/',1)[:-1][0]
                self.output_filename_nodir = "." + "/" + self.output_filename.rsplit('/',1)[-1]

    def always_run_simulation(self,cmdline):
        print "Starting:\n>>"+cmdline+"\nin directory '"+self.input_dir+"'"      
        simulation = subprocess.Popen(
                cmdline, stdout=subprocess.PIPE,shell=True,stderr=subprocess.PIPE,cwd=self.input_dir
            )
        while True:
            out = simulation.stdout.readline()
            #errout = simulation.stderr.readline()
            if out == '' and simulation.poll() is not None:
                break
            if out != '':
                sys.stdout.write(out)
                sys.stdout.flush()
        self.output=self.input_dir+"/"+self.inputs.outputFilename
        

        def run_simulation(self,cmdline):
            if self.output_filename is not None:
                self.always_run_simulation(cmdline)

    
    

    @property
    def particle_flux(self):
        return self.outputs[self.group_name+self.particle_flux_name][()]

    @property
    def heat_flux(self):
        return self.outputs[self.group_name+self.heat_flux_name][()]

    @property
    def VPrimeHat(self):
        return self.outputs[self.group_name+self.VPrimeHat_name][()]

    @property
    def heat_source(self):
        return self.outputs[self.group_name+self.heat_source_name][()]

    @property
    def particle_source(self):
        return self.outputs[self.group_name+self.particle_source_name][()]

    @property
    def masses(self):
        return self.inputs.masses
    
    @property
    def THat(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"THats"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.THat_name][()]
        except KeyError:
            print "THat could not be obtained since no external profiles have been specified and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def nHat(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"nHats"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.nHat_name][()]
        except KeyError:
            print "nHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."
    
    @property
    def Delta(self):
        return self.inputs.Delta

    @property
    def omega(self):
        return self.inputs.omega

    @property
    def nu_r(self):
        return self.inputs.nu_r

    @property
    def psi(self):
        return self.outputs[self.group_name+self.psi_name][()]
        #self.inputs.psi

    @property
    def num_species(self):
        return self.outputs[self.group_name+self.num_species_name][()]

    
    def copy_simulation_to_dir(self,dir):
        #try to make the new directory or don't if it exists
        try:
            os.makedirs(dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        #For now, filename of the new will be the same as the old
        #but in the new directory.

        #the new output file
        if self.output_filename_nodir is not None:
            new_o_filename=dir+"/"+self.output_filename_nodir
            shutil.copy(self.output_filename,new_o_filename)
        else:
            new_o_filename=None
 
        #the new input file
        #assume that the input file always exists
        new_i_filename=dir+"/"+self.input_filename_nodir
        shutil.copy(self.input_filename,new_i_filename)
 
       
        #construct a new simulation with the copied inputs, and set some attributes
        sim_copy=perfect_simulation(new_i_filename,new_o_filename,self.group_name)
        sim_copy.description=self.description
        sim_copy.species=self.species

        return sim_copy
        
        





#####################################

class normalized_perfect_simulation(perfect_simulation):
    def __init__(self,input_filename,norm_filename,output_filename=sentinel,group_name=sentinel):
        perfect_simulation.__init__(self,input_filename,output_filename,group_name)

        self.norm=norm_filename

        tolerance=0.001
        #in the future, it might be convenient to change the quantities in the simulation to match those derived from the normfile
        if not self.verify_Delta(tolerance):
            print "Warning: Delta in simulation does not match Delta due to normalized quantities to relative precision " + str(tolerance)

        if not self.verify_omega(tolerance):
            print "Warning: omega in simulation does not match omega due to normalized quantities to relative precision " + str(tolerance)

    @property
    def norm(self):
        return self.norm_filename + "bajs"

    @norm.setter
    def norm(self,norm_filename):
        self.norm_filename="." + "/" + norm_filename
        self.normfile = f90nml.read(self.norm_filename)
        self.norm_dir = self.norm_filename.rsplit('/',1)[:-1][0]
        self.norm_filename_nodir = "." + "/" + self.norm_filename.rsplit('/',1)[-1]
        
        self.units=self.normfile["unitSpecification"]["units"]
        norm_group_name="normalizationParameters"
        self.BBar=self.normfile[norm_group_name]["BBar"]
        self.RBar=self.normfile[norm_group_name]["RBar"]
        self.TBar=self.normfile[norm_group_name]["TBar"]
        self.mBar=self.normfile[norm_group_name]["mBar"]
        self.nBar=self.normfile[norm_group_name]["nBar"]
        self.eBar=self.normfile[norm_group_name]["eBar"]
        self.ePhiBar=self.normfile[norm_group_name]["ePhiBar"]
        self.vBar=numpy.sqrt(2*self.TBar/float(self.mBar))

        
            
    def verify_Delta(self,tolerance):
        if self.units=="SI":
            Delta_norm=self.mBar*self.vBar/float(self.eBar*self.BBar*self.RBar)
        else:
            print "Units not supported, cannot calculate Delta."
        Delta_sim=self.Delta
        print "Delta in simulation: " + str(Delta_sim)
        print "Delta calculated from normalization constants: " + str(Delta_norm)
        if (Delta_sim <= (1+tolerance)*Delta_norm)and (Delta_sim >= (1-tolerance)*Delta_norm):

            return True
        else:
            return False

    def verify_omega(self,tolerance):
        if self.units=="SI":
            omega_norm=self.ePhiBar/float(self.vBar*self.eBar*self.BBar*self.RBar)
        else:
            print "Units not supported, cannot calculate omega."
        omega_sim=self.omega
        print "omega in simulation: " + str(omega_sim)
        print "omega calculated from normalization constants: " + str(omega_norm)
        if (omega_sim <= (1+tolerance)*omega_norm)and (omega_sim >= (1-tolerance)*omega_norm):
            return True
        else:
            return False

    #def verify_nu_r(self,tolerance,logLambda):
    #    if self.units=="SI":
    #        ep0=scipy.constants.epsilon_0
    #        prefactor=(2**(1.0/2.0)/(12*numpy.pi**(3.0/2.0))) * (self.eBar**4/(ep0**2))
    #        nu_r_norm=prefactor*logLambda*(self.RBar*self.nBar/self.TBar**2)
    #    else:
    #        print "Units not supported, cannot calculate omega."
    #    nu_r_sim=self.nu_r
    #    print "nu_r in simulation: " + str(nu_r_sim)
    #    print "nu_r calculated from normalization constants: " + str(omega_norm)
    #    if (nu_r_sim <= (1+tolerance)*nu_r_norm)and (nu_r_sim >= (1-tolerance)*nu_r_norm):
    #        return True
    #    else:
    #        return False
        

    def copy_simulation_norm_to_dir(self,dir):
        #copies norm in addition to filename
        sim_copy=normalized_perfect_simulation.copy_simulation_to_dir(self,dir)
        sim_copy.__class__=type(self)
        new_n_filename=dir+"/"+self.norm_filename_nodir
        shutil.copy(self.norm_filename,new_n_filename)
        sim_copy.norm=new_n_filename
        
        return sim_copy
        

    @property
    def T(self):
        return self.TBar*self.THat

    @property
    def n(self):
        return self.nBar*self.nHat

    @property
    def normed_particle_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        return (self.particle_flux/VPrimeHat)*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_heat_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        return (self.heat_flux/VPrimeHat)*self.TBar*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_heat_source(self):
        #print self.heat_source
        #print self.masses
        #print self.heat_source/self.masses
        return self.Delta*self.nBar/(self.vBar**2*self.RBar*numpy.sqrt(self.masses))*self.heat_source

    @property
    def normed_particle_source(self):
        #print self.particle_source
        #print self.masses
        #print self.particle_source/self.masses
        return self.Delta*self.nBar/(self.vBar**2*self.RBar*numpy.sqrt(self.masses))*self.particle_source

    @property
    def normed_conductive_heat_flux(self):
        return self.normed_heat_flux-(5.0/2.0)*self.T*self.normed_particle_flux


    

