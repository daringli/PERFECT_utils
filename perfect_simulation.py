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
from get_index_range import get_index_range
from nstack import nstack



sentinel=object() #to detect default inputs

class perfect_simulation(object):

    def __init__(self,input_filename,output_filename=sentinel,species_filename=sentinel,group_name=sentinel,pedestal_points=[]):
        
        #description of simulation, for usage as legend in plots, etc.
        self.description=""
        #read species from file in simulation dir
        #relevant field can also be set by chaning self.species_list
        if species_filename is not sentinel:
            self.species=species_filename
        

        #getting input
        self.input=input_filename
        #self.input_filename=input_filename
        #self.inputs=perfectInput(self.input_filename)
        #self.input_dir = "." + "/" + self.input_filename.rsplit('/',1)[:-1][0]

        #having no output is acceptable,
        #it just means that simulation has yet to be run
        self.output=output_filename

        #where the pedestal starts and stops
        self.pedestal_points=pedestal_points

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
        self.psiAHat_name="psiAHat"
        self.theta_name="theta"
        self.particle_flux_name="particleFlux"
        self.heat_flux_name="heatFlux"
        self.momentum_flux_name="momentumFlux"
        self.VPrimeHat_name="VPrimeHat"
        self.Delta_name="Delta"
        self.omega_name="omega"
        self.THat_name="THat"
        self.dTHatdpsiN_name="d(THat)d(psi)"
        self.nHat_name="nHat"
        self.etaHat_name="etaHat"
        self.detaHatdpsiN_name="d(etaHat)d(psi)"
        self.PhiHat_name="PhiHat"
        self.dPhiHatdpsiN_name="d(PhiHat)d(psi)"
        self.U_name="U"
        self.deltaN_name="deltaN"
        self.deltaEta_name="deltaEta"
        self.deltaT_name="deltaT"
        self.num_species_name="Nspecies"
        self.heat_source_name="heatSourceProfile"
        self.particle_source_name="particleSourceProfile"

        self.flow_name="flow"
        self.kPar_name="kPar"
        self.FSAkPar_name="FSAKPar"
        self.FSAFlow_name="FSAFlow"
        self.flow_inboard_name="flowInboard"
        self.flow_outboard_name="flowOutboard"
        self.kPar_inboard_name="kParInboard"
        self.kPar_outboard_name="kParOutboard"
        self.collisionality_name="nuPrimeProfile"
        self.FSABFlow_name="FSABFlow"

        self.density_perturbation_name="densityPerturbation"

        self.local_name="makeLocalApproximation"

        print "Total source charge in domain:" + str(self.total_source_charge)
        print "final j:" + str(self.jHat[-1])
        
    @property
    def species(self):
        return self.species_list
    
    @species.setter
    def species(self,species_filename):
        #read species files on the format species1,species2,species3,...
        #should be on one line. can contain whitespaces
        self.species_filename="." + "/" + species_filename
        self.species_dir = self.species_filename.rsplit('/',1)[:-1][0]
        self.species_filename_nodir = "." + "/" + self.species_filename.rsplit('/',1)[-1]
        print open(species_filename,'r').read()
        self.species_list=[x.strip() for x in open(species_filename,'r').read().split('\n')[0].split(',')]
        
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


    def attrib_at_psi_of_theta(self,attrib,psiN):
        indices=get_index_range(self.psi,[psiN,psiN],ret_range=False)
        return getattr(self,attrib)[indices[0],:,:]
                
    def attrib_max_psi_of_theta(self,attrib,xlim=None):
        #help function to return the value of psiN that maximizes an attribute
        #for each given theta
        if xlim != None:
            print "xlim: " + str(xlim)
            indices=get_index_range(self.psi,xlim,ret_range=True)
            psi=self.psi
        else:
            psi=self.psi
            indices=range(len(psi))
        print [psi[indices[0]],psi[indices[-1]]]
        data=getattr(self,attrib)
        max_indices=numpy.zeros(data[0,:,:].shape)
        max_psiN=numpy.zeros(data[0,:,:].shape)
        for i_sp in range(len(data[0,0,:])):
            #for each species
            for i_th in range(len(data[0,:,i_sp])):
                #for each theta
                print numpy.argmax(data[indices,i_th,i_sp])
                max_indices[i_th,i_sp] = numpy.argmax(data[indices,i_th,i_sp])
                max_psiN[i_th,i_sp]=psi[max_indices[i_th,i_sp]]
        #print max_indices
        print [numpy.min(max_psiN),numpy.max(max_psiN)]
        return max_psiN
    
    @property
    def local(self):
        if self.outputs[self.group_name+self.local_name][()] == 1:
            return True
        elif self.outputs[self.group_name+self.local_name][()] == -1:
            return False
        else:
            print "perfect_simulation: error: makeLocalApproximation is neither true nor false!?"
    @property
    def Z(self):
        return numpy.array(self.inputs.charges)

    @property
    def particle_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.particle_flux_name][()]*signOfVPrimeHat

    @property
    def momentum_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.momentum_flux_name][()]*signOfVPrimeHat
    
    @property
    def heat_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.heat_flux_name][()]*signOfVPrimeHat


    @property
    def conductive_heat_flux(self):        
        return self.heat_flux-(5.0/2.0)*self.THat*self.particle_flux

    @property
    def particle_flux_over_nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.particle_flux/npeds
    
    
    @property
    def momentum_flux_over_nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.momentum_flux/npeds

    @property
    def heat_flux_over_nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.heat_flux/npeds

    @property
    def conductive_heat_flux_over_nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.conductive_heat_flux/npeds
    
    @property
    def jHat(self):
        return numpy.sum(self.Z*self.particle_flux,axis=1)

    @property
    def momentum_flux_sum(self):
        return numpy.sum(self.masses*self.momentum_flux,axis=1)

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
    def heat_source_over_m2(self):
        return self.outputs[self.group_name+self.heat_source_name][()]/(self.masses**2)

    @property
    def particle_source_over_m2(self):
        return self.outputs[self.group_name+self.particle_source_name][()]/(self.masses**2)

    @property
    def particle_source_over_m2nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.particle_source_over_m2/npeds

    @property
    def heat_source_over_m2nPed(self):
        psiN_point=0.9
        psiN_index=get_index_range(self.psi,[psiN_point,psiN_point])[1]
        npeds=[simul.nHat[i,psiN_index] for i in range(self.num_species)]
        return self.heat_source_over_m2/npeds
    
    
    #@property
    #def THat32_new_massfac_particle_source(self):
    #    #for comparrison with dGammaHat/dPsiN
    #    VPrimeHat=[self.VPrimeHat]
    #    VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
    #    VPrimeHat=numpy.fabs(VPrimeHat)
    #    return 1.6480*self.THat**(3./2.)*self.psiAHat*VPrimeHat/(self.Delta)*self.particle_source

    @property
    def GammaHat_source(self):
        #for comparrison with dGammaHat/dPsiN
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        VPrimeHat=numpy.fabs(VPrimeHat)
        return 2*(numpy.sqrt(numpy.pi)/2)*self.THat**(3./2.)*self.psiAHat*VPrimeHat/(self.Delta)*self.particle_source/(self.masses**2)

    @property
    def qHat_source(self):
        #for comparrison with dGammaHat/dPsiN
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        VPrimeHat=numpy.fabs(VPrimeHat)
        dPhiHatdpsiN=[self.dPhiHatdpsiN]
        dPhiHatdpsiN=numpy.dot(numpy.transpose(dPhiHatdpsiN),[numpy.array([1]*self.num_species)])
        
        
        return -2.0*self.THat**(5./2.)*self.psiAHat*VPrimeHat/(self.Delta)*((5*numpy.sqrt(numpy.pi)/4)*self.particle_source + (3.0*numpy.sqrt(numpy.pi)/4.0)*self.heat_source)/(self.masses**2) -((2*self.omega/self.Delta)*self.Z*dPhiHatdpsiN + (5./2.)*self.dTHatdpsiN)*self.particle_flux

    @property
    def FSAFlow(self):
        return self.outputs[self.group_name+self.FSAFlow_name][()]
    
    
    @property
    def flow_inboard(self):
        return self.outputs[self.group_name+self.flow_inboard_name][()]

    @property
    def flow_outboard(self):
        return self.outputs[self.group_name+self.flow_outboard_name][()]

    @property
    def FSAkPar(self):
        return self.outputs[self.group_name+self.FSAkPar_name][()]
    
    
    @property
    def kPar_inboard(self):
        return self.outputs[self.group_name+self.kPar_inboard_name][()]
    
    @property
    def kPar_outboard(self):
        return self.outputs[self.group_name+self.kPar_outboard_name][()]

    @property
    def FSABFlow(self):
        return self.outputs[self.group_name+self.FSABFlow_name][()]

    @property
    def flow(self):
        return self.outputs[self.group_name+self.flow_name][()]

    @property
    def flow_minus_FSAFlow(self):
        return self.flow-numpy.expand_dims(self.FSAFlow,axis=1)

    
    @property
    def flow_over_vT(self):
        return self.Delta*nstack(numpy.sqrt(self.masses/self.THat),axis=1,n=len(self.theta))*self.outputs[self.group_name+self.flow_name][()]

    @property
    def flow_difference(self):
        if self.num_species>1:
            main_index=0
            impurity_index=1
            return self.flow[:,:,main_index] - self.flow[:,:,impurity_index]
        else:
            print "perfect_simulation: flow_difference: warning: no impurity species. Difference will be just the main ion flow"
            return self.flow[:,:,main_index]

    @property
    def flow_max_psi_of_theta(self):
        return self.attrib_max_psi_of_theta("flow",[0.91,0.97463697978512787])

    @property
    def flow_at_psi_of_theta(self):
        return self.attrib_at_psi_of_theta("flow",0.955)

    
    @property
    def kPar(self):
        return self.outputs[self.group_name+self.kPar_name][()]

    @property
    def kPar_minus_FSAkPar(self):
        return self.kPar-numpy.expand_dims(self.FSAkPar,axis=1)

    
    @property
    def kPar_max_psi_of_theta(self):
        return self.attrib_max_psi_of_theta("kPar",[0.9,0.97463697978512787])

    @property
    def kPar_at_psi_of_theta(self):
        return self.attrib_at_psi_of_theta("kPar",0.955)

    

    @property
    def FSABJPar(self):
        return numpy.sum(self.Z*(self.nHat*self.FSABFlow),axis=1)

    @property
    def potential_perturbation(self):
        prefactor=numpy.sum((self.Z**2)*(self.nHat/self.THat),axis=1)**(-1)*self.Delta/(2*self.omega)
        prefactor=prefactor[:,numpy.newaxis]
        Zn=numpy.expand_dims(self.Z*self.nHat, axis=1)
        #print Zn.shape
        #print self.density_perturbation.shape
        #print prefactor.shape
        #print (prefactor*numpy.sum(Zn*self.density_perturbation,axis=2)).shape
        return prefactor*numpy.sum(Zn*self.density_perturbation,axis=2)
    
    @property
    def U(self):
        return self.outputs[self.group_name+self.U_name][()]

    @property
    def deltaN(self):
        return self.outputs[self.group_name+self.deltaN_name][()]

    @property
    def deltaEta(self):
        return self.outputs[self.group_name+self.deltaEta_name][()]

    @property
    def deltaT(self):
        return self.outputs[self.group_name+self.deltaT_name][()]

    @property
    def masses(self):
        return numpy.array(self.inputs.masses)
    
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
    def dTHatdpsiN(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"dTHatdpsis"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.dTHatdpsiN_name][()]
        except KeyError:
            print "PhiHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."


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
    def etaHat(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"etaHats"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.etaHat_name][()]
        except KeyError:
            print "etaHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def detaHatdpsiN(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"detaHatdpsis"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.detaHatdpsiN_name][()]
        except KeyError:
            print "etaHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."
    

    @property
    def PhiHat(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"PhiHat"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.PhiHat_name][()]
        except KeyError:
            print "PhiHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def dPhiHatdpsiN(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"dPhiHatdpsi"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.dPhiHatdpsiN_name][()]
        except KeyError:
            print "PhiHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def dPhiHat_dpsiN(self):
        #intentionally "crude" to use same deriv approx as dGamma_dpsiN
        psiN=self.psi
        PhiHat=self.PhiHat
        ret=numpy.zeros([len(psiN)-1,self.num_species])
        for i in range(1,len(psiN)):
            dpsiN=psiN[i]-psiN[i-1]
            dPhiHat=PhiHat[i]-PhiHat[i-1]
            #print dPhiHat/dpsiN
            ret[i-1]=dPhiHat/dpsiN
        return ret

    @property
    def dpsiN(self):
        return self.psi[1]-self.psi[0]
    
    @property
    def dTHat_dpsiN(self):
        #intentionally "crude" to use same deriv approx as dGamma_dpsiN
        psiN=self.psi
        THat=self.THat
        ret=numpy.zeros([len(psiN)-1,self.num_species])
        for i in range(1,len(psiN)):
            dpsiN=psiN[i]-psiN[i-1]
            dTHat=THat[i]-THat[i-1]
            #print dTHat/dpsiN
            ret[i-1]=dTHat/dpsiN
        return ret

    @property
    def alpha(self):
        #Zeff-Zmain
        if self.num_species>2:
            electron_index=self.species_list.index("e")
            Z=self.Z
            Z[electron_index]=0 # to exclude electrons from calculation
            return numpy.sum(Z**2*self.nHat,axis=1)/numpy.sum(Z*self.nHat,axis=1)-Z[0]
        else:
            print "perfect_simulation: alpha: warning: no impurity species. alpha=Zeff-Z_main set to zero."
            return [0]*len(self.psi)


    @property
    def source_charge(self):
        return numpy.sum(self.Z*self.GammaHat_source,axis=1)
    
    @property
    def total_source_charge(self):
        final_i=-30
        source_charge=numpy.sum(self.Z*self.GammaHat_source,axis=1)
        return scipy.integrate.simps(source_charge[0:final_i],self.psi[0:final_i])
    
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
    def psiAHat(self):
        return self.outputs[self.group_name+self.psiAHat_name][()]
        #self.inputs.psi

    @property
    def theta(self):
        return self.outputs[self.group_name+self.theta_name][()]

    @property
    def num_species(self):
        return self.outputs[self.group_name+self.num_species_name][()]

    @property
    def collisionality(self):
        return self.outputs[self.group_name+self.collisionality_name][()]

    @property
    def density_perturbation(self):
        #non-adiabatic part
        return self.outputs[self.group_name+self.density_perturbation_name][()]

    @property
    def total_density_perturbation(self):
        #non-adiabatic part
        ZoT=numpy.expand_dims(self.Z/self.THat, axis=1)
        ret=self.density_perturbation - (2*self.omega/self.Delta)*ZoT*self.potential_perturbation[:,:,numpy.newaxis]
        return ret

    @property
    def dGammaHat_dpsiN(self):
        psiN=self.psi
        GammaHat=self.particle_flux
        ret=numpy.zeros([len(psiN),self.num_species])
        for i in range(1,len(psiN)):
            dpsiN=psiN[i]-psiN[i-1]
            dGammaHat=GammaHat[i]-GammaHat[i-1]
            #print dGammaHat/dpsiN
            ret[i-1]=dGammaHat/dpsiN
        ret[0]=0
        return ret

    @property
    def djHat_dpsiN(self):
        return numpy.sum(self.Z*self.dGammaHat_dpsiN,axis=1)



    @property
    def dqHat_dpsiN(self):
        #non-adiabatic part
        psiN=self.psi
        qHat=self.conductive_heat_flux
        ret=numpy.zeros([len(psiN),self.num_species])
        for i in range(1,len(psiN)):
            dpsiN=psiN[i]-psiN[i-1]
            dqHat=qHat[i]-qHat[i-1]
            #print dGammaHat/dpsiN
            ret[i-1]=dqHat/dpsiN
        ret[0]=0
        return ret

    
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
    def __init__(self,input_filename,norm_filename,species_filename=sentinel,output_filename=sentinel,group_name=sentinel,pedestal_points=[]):
        perfect_simulation.__init__(self,input_filename,output_filename,species_filename,group_name,pedestal_points)

        self.norm=norm_filename

        tolerance=0.005
        #in the future, it might be convenient to change the quantities in the simulation to match those derived from the normfile
        if not self.verify_Delta(tolerance):
            print "Warning: Delta in simulation does not match Delta due to normalized quantities to relative precision " + str(tolerance)

        if not self.verify_omega(tolerance):
            print "Warning: omega in simulation does not match omega due to normalized quantities to relative precision " + str(tolerance)

    @property
    def norm(self):
        return self.norm_filename

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
    def eta(self):
        return self.nBar*self.etaHat

    @property
    def Phi(self):
        return self.ePhiBar*self.PhiHat/self.eBar


    @property
    def normed_particle_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.particle_flux/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_momentum_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.momentum_flux/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*(self.RBar**2)*self.BBar*self.nBar*self.vBar**2

    @property
    def normed_heat_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.heat_flux/VPrimeHat)*signOfVPrimeHat*self.TBar*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_conductive_heat_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        return (self.conductive_heat_flux/VPrimeHat)*self.TBar*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar
    #self.normed_heat_flux-(5.0/2.0)*self.T*self.normed_particle_flux

    @property
    def normed_ambipolariy(self):
        return (self.ambipolarity/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar*self.eBar

    @property
    def normed_momentum_flux_sum(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.momentum_flux_sum/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*(self.RBar**2)*self.BBar*self.nBar*self.vBar**2
    
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
    def normed_FSABFlow(self):
        return self.BBar*self.vBar*self.Delta*self.FSABFlow
    
    @property
    def normed_FSABJPar(self):
        return self.nBar*self.BBar*self.eBar*self.vBar*self.Delta*self.FSABJPar
    
    @property
    def normed_flow_inboard(self):
        return self.Delta*self.vBar*self.flow_inboard


    @property
    def normed_flow_outboard(self):
        return self.Delta*self.vBar*self.flow_outboard

    @property
    def normed_flow(self):
        return self.Delta*self.vBar*self.flow
    
    @property
    def normed_potential_perturbation(self):
        return self.PhiBar*self.potential_perturbation

    

