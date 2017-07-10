from __future__ import division

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
import scipy.ndimage
import subprocess
from subprocess import call
from get_index_range import get_index_range
from nstack import nstack
from diff_matrix import diff_matrix
from array_rank import arraylist_rank
from peakfinder import PPV,PPV_x


class perfect_simulation(object):

    def __init__(self,input_filename,output_filename=None,species_filename=None,psiN_to_psi_filename=None,global_term_multiplier_filename=None,group_name=None,pedestal_start_stop=(None,None),pedestal_point = None,core_point=None):
        
        #description of simulation, for usage as legend in plots, etc.
        self.description=""
        #read species from file in simulation dir
        #relevant field can also be set by chaning self.species_list
        if species_filename is not None:
            self.species=species_filename


        #Output intepretation and group specification
        #Different group in the same output can in principle
        #be handled by different simulation objects
        #hiding implementation details
        self.did_it_converge_name = "didItConverge"
        self.psi_name="psi"
        self.psiAHat_name="psiAHat"
        self.theta_name="theta"
        self.particle_flux_name="particleFlux"
        self.particle_flux_before_theta_integral_name="particleFluxBeforeThetaIntegral"
        self.heat_flux_name="heatFlux"
        self.heat_flux_before_theta_integral_name="heatFluxBeforeThetaIntegral"
        self.momentum_flux_name="momentumFlux"
        self.momentum_flux_before_theta_integral_name="momentumFluxBeforeThetaIntegral"
        self.VPrimeHat_name="VPrimeHat"
        self.Delta_name="Delta"
        self.omega_name="omega"
        self.THat_name="THat"
        self.dTHatdpsiN_name="d(THat)d(psi)"
        self.nHat_name="nHat"
        self.dnHatdpsiN_name="d(nHat)d(psi)"
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
        self.momentum_source_name="momentumSourceProfile"
        self.particle_source_name="particleSourceProfile"
        self.no_charge_source_name="noChargeSource"
        self.no_charge_source_momentum_source_name="noChargeSourceMomentumSourceProfile"
        self.no_charge_source_momentum_source_species_dependence_name = "momentumSourceSpeciesDependence"
        self.no_charge_source_particle_source_name="noChargeSourceParticleSourceProfile"
        self.species_indep_source_particle_source_name="speciesIndepSourceParticleSourceProfile"
        self.species_indep_source_heat_source_name="speciesIndepSourceHeatSourceProfile"
        self.species_indep_source_momentum_source_name="speciesIndepSourceMomentumSourceProfile"
        self.constant_source_particle_source_name="constantSourceParticleSourceProfile"
        self.constant_source_heat_source_name="constantSourceHeatSourceProfile"
        self.constant_source_momentum_source_name="constantSourceMomentumSourceProfile"
        
        self.flow_name="flow"
        self.toroidal_flow_name="toroidalFlow"
        self.poloidal_flow_name="poloidalFlow"
        self.FSA_flow_name="FSAFlow"
        self.FSA_toroidal_flow_name="FSAToroidalFlow"
        self.FSA_poloidal_flow_name="FSAPoloidalFlow"

        #old names, for compatibility with old output
        self.magnetization_flow_perturbation_name="magnetizationFlowPerturbation"
        self.magnetization_perturbation_name="magnetizationPerturbation"
        #new names for same thing
        self.p_perp_term_in_Vp_name="pPerpTermInVp"
        self.p_perp_term_in_Vp_underived_name="pPerpTermInVpBeforePsiDerivative"
        
        
        self.kPar_name="kPar"
        self.FSAkPar_name="FSAKPar"
        self.FSAFlow_name="FSAFlow"
        self.flow_inboard_name="flowInboard"
        self.flow_outboard_name="flowOutboard"
        self.kPar_inboard_name="kParInboard"
        self.kPar_outboard_name="kParOutboard"
        self.collisionality_name="nuPrimeProfile"
        self.FSABFlow_name="FSABFlow"

        self.RHat_name="RHat"
        self.JHat_name="JHat"
        self.IHat_name="IHat"
        self.dIHatdpsiN_name="d(IHat)d(psi)"
        self.BHat_name="BHat"
        self.dBHatdpsiN_name="d(BHat)d(psi)"
        self.dBHatdtheta_name="d(BHat)d(theta)"
        self.FSABHat2_name="FSABHat2"
        self.q_name="Miller_q"
        self.epsilon_name="epsil"

        self.density_perturbation_name="densityPerturbation"
        self.pressure_perturbation_name="pressurePerturbation"

        self.FSA_density_perturbation_name="FSADensityPerturbation"
        self.FSA_pressure_perturbation_name="FSAPressurePerturbation"

        self.Nxi_for_x_name="Nxi_for_x"
        
        self.local_name="makeLocalApproximation"
        self.include_ddpsi_name= "includeddpsiTerm"
        
        self.Npsi_sourceless_right_name = "NpsiSourcelessRight"
        self.Npsi_sourceless_left_name = "NpsiSourcelessLeft"
        
        self.version_group="/"
        self.program_mode_group="/"
        if group_name!=None:
            self.group_name=group_name
        else:
            self.group_name="/run  1/"
            
        #getting input
        self.input=input_filename
        
        #having no output is acceptable,
        #it just means that simulation has yet to be run
        self.output=output_filename

        #can have no psiN_to_psi, if uniform grid.
        self.psiN_to_psi = psiN_to_psi_filename

        #can have no global_term_multiplier, if useGlobalTermMultiplier is 0.
        self.global_term_multiplier = global_term_multiplier_filename
        
        
        #where the pedestal starts and stops
        self.pedestal_start_stop=pedestal_start_stop

        # used to evaluate quantities at a fixed radius in the pedestal
        self.pedestal_point = pedestal_point

        # used to evaluate quantities at a fixed radius in the core
        self.core_point = core_point
        

        #print "Total source charge in domain:" + str(self.total_source_charge)
        #print "final j:" + str(self.jHat[-1])
        
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
        #print open(species_filename,'r').read()
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
        if (output_filename==None) or (output_filename==None):
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

                
    @property
    def did_it_converge(self):
        conv = self.outputs[self.group_name+self.did_it_converge_name][()]
        if conv == -1:
            return False
        elif conv == 1:
            return True
        else:
            print "perfect_simulation ERROR: Invalid convergence status."
            
    @property
    def psiN_to_psi(self):
        return self.psiN_to_psi_filename

    @psiN_to_psi.setter
    def psiN_to_psi(self,psiN_to_psi_filename):
        try:
            self.inputs.psiGridType
        except AttributeError:
            self.inputs.psiGridType=0
            
        if self.inputs.psiGridType==1:
            if (psiN_to_psi_filename==None) or (psiN_to_psi_filename==None):
                self.psiN_to_psi_filename=None
                self.psiN_to_psi_profiles=None
                self.psiN_to_psi_dir=self.input_dir
                self.psiN_to_psi_filename_nodir =None
            else:
                self.psiN_to_psi_filename="." + "/" + psiN_to_psi_filename
                self.psiN_to_psi_dir = self.psiN_to_psi_filename.rsplit('/',1)[:-1][0]
                self.psiN_to_psi_filename_nodir = "." + "/" + self.psiN_to_psi_filename.rsplit('/',1)[-1]
                self.psiN_to_psi_groupname="/Npsi"+str(self.inputs.Npsi)+"/"
                if os.path.isfile(self.psiN_to_psi_filename):
                    self.psiN_to_psi_profiles=h5py.File(psiN_to_psi_filename,'r')
                    self.psiAHatArray = self.psiN_to_psi_profiles[self.psiN_to_psi_groupname+"psiAHatArray"][()]
                    self.actual_psi = self.psiN_to_psi_profiles[self.psiN_to_psi_groupname+"psiArray"][()]
                    self.actual_psiN = self.psiN_to_psi_profiles[self.psiN_to_psi_groupname+"psiArray"][()]/self.psiAHat
                else:
                    print "Error: perfect_simulation: psiN_to_psi cannot be read, file "+ psiN_to_psi_filename +" does not exist!"
                    sys.exit(1)
        elif self.inputs.psiGridType==0:
            self.psiAHatArray = numpy.array([self.psiAHat]*self.Npsi)
            self.actual_psi = self.psiAHat*self.psi
            self.actual_psiN = self.psi

    @property
    def psiAHatArray_normed(self):
        return self.psiAHatArray*self.psiAHat
            
    @property
    def global_term_multiplier(self):
        return self.global_term_multiplier_profile

    @global_term_multiplier.setter
    def global_term_multiplier(self,global_term_multiplier_filename):
        try:
            self.inputs.useGlobalTermMultiplier
        except AttributeError:
            self.inputs.useGlobalTermMultiplier=0
            
        if self.inputs.useGlobalTermMultiplier==1:
            if (global_term_multiplier_filename==None) or (global_term_multiplier_filename==None):
                self.global_term_multiplier_filename=None
                self.global_term_multiplier_profile=None
                self.global_term_multiplier_dir=self.input_dir
                self.global_term_multiplier_filename_nodir =None
            else:
                self.global_term_multiplier_filename="." + "/" + global_term_multiplier_filename
                self.global_term_multiplier_dir = self.global_term_multiplier_filename.rsplit('/',1)[:-1][0]
                self.global_term_multiplier_filename_nodir = "." + "/" + self.global_term_multiplier_filename.rsplit('/',1)[-1]
                self.global_term_multiplier_groupname="/Npsi"+str(self.inputs.Npsi)+"/"
                if os.path.isfile(self.global_term_multiplier_filename):
                    self.global_term_multiplier_file=h5py.File(global_term_multiplier_filename,'r')
                    self.global_term_multiplier_profile = self.global_term_multiplier_file[self.global_term_multiplier_groupname+"globalTermMultiplier"][()]
                else:
                    print "Error: perfect_simulation: global_term_multiplier cannot be read, file "+ global_term_multiplier_filename +" does not exist!"
                    sys.exit(1)
        elif self.inputs.useGlobalTermMultiplier==0:
            self.global_term_multiplier_profile = numpy.array([1.0]*self.Npsi)
            
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


    def dimension_check(self,arr):
        #check whether an array has appropriate dimensions
        
        pass
            
    def FSA(self,attrib,jacobian=True,periodic=True):
        if type(attrib) == str:
            a=getattr(self,attrib)
        else:
            a=attrib
        rank=arraylist_rank(a)
        if rank==3:
            #sanity check
            if len(a[1]) == self.Ntheta:
                axis = 1
                VPH=numpy.expand_dims(numpy.fabs(self.VPrimeHat),axis=1)
            else:
                print "psi, theta, species dependent quantity has axis in wrong order?? Or v-space dependence?"
            if jacobian==True:
                j=numpy.expand_dims(1/numpy.fabs(self.JHat),axis=2)
            else:
                j=1 
        elif rank==2:
            if len(a[1]) == self.Ntheta:
                print "OK!!!!!!!!!!!!!!!!!!!!"
                #psi, theta dependence
                axis = 1
                if jacobian==True:
                    j=1/numpy.fabs(self.JHat)
                else:
                    j=1
            elif len(a[0]) == self.Ntheta:
                #theta, species dependence
                #add psi dimension since jHat can be psi dependent
                axis = 1
                a=numpy.expand_dims(a,axis=0)
                if jacobian==True:
                    j=numpy.expand_dims(1/numpy.fabs(self.JHat),axis=2)
                else:
                    j=1
            else:
                #psi, species dependence (or weird thigns)
                return a
            VPH=numpy.fabs(self.VPrimeHat)
            
        elif rank==1:
            if len(a) == self.Ntheta:
                axis = 0
                a=numpy.expand_dims(a,axis=0) #since VPrime can depend on psi
                VPH=numpy.fabs(self.VPrimeHat)
            else:
                return a
            if jacobian==True:
                j=1/numpy.fabs(self.JHat)
            else:
                j=1
        else:
            print "Error, FSA can only take things with 3 indices."
        if periodic==False:
            #NOTE: should not be used, only to compare new with old buggy
            #for documentation purposes.
            #Would be approriate for non-periodic subintervals
            return numpy.trapz(a*j,dx=self.dtheta,axis=axis)/VPH
        return numpy.sum(a*j*self.dtheta,axis=axis)/VPH

    def ddpsiN(self,attrib,order=4,scale=True):
        #dds: derivative on uniform grid. supported orders are 2 and 4 (2016-06-27)
        dds=diff_matrix(self.psi[0],self.psi[-1],self.Npsi,order=order)
        if scale:
            scale_factor=(self.psiAHat/self.psiAHatArray)
        else:
            scale_factor=numpy.ones(self.Npsi) #for derivatives on internal grid
            #scale_factor=(self.psiAHat/self.psiAHatArray[0])*numpy.ones(self.Npsi)
        #get attribute array
        if type(attrib) == str:
            a = getattr(self,attrib)
        else:
            a = attrib
        # convert from dds to ddpsiN for different array ranks
        rank=arraylist_rank(a)
        if rank == 3:
            #sanity check. Should have shape Npsi,Ntheta,Nspecies
            if a.shape == (self.Npsi,self.Ntheta,self.Nspecies):
                return scale_factor[:,numpy.newaxis,numpy.newaxis]*numpy.tensordot(dds,a,([1],[0]))
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 3 array badly shaped."
                sys.exit(1)
        if rank == 2:
            #sanity check. Either (Npsi,Nspecies) or (Npsi,Ntheta) shaped
            if (a.shape == (self.Npsi,self.Ntheta)) or (a.shape == (self.Npsi,self.Nspecies)):
                return scale_factor[:,numpy.newaxis]*numpy.tensordot(dds,a,([1],[0]))
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 2 array badly shaped."
                sys.exit(1)
        if rank == 1:
            #sanity check. (Npsi,) shaped
            if a.shape == (self.Npsi,):
                return scale_factor*numpy.tensordot(dds,a,([1],[0]))
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 1 array badly shaped."
                sys.exit(1)

    def ddtheta(self,attrib,order=4,periodic=True,last=False):
        dds=diff_matrix(self.theta[0],self.theta[0]+2*numpy.pi,self.Ntheta,order=order,periodic=periodic,endpoint=False)
        if type(attrib) == str:
            a = getattr(self,attrib)
        else:
            a = attrib
        # convert from dds to ddpsiN for different array ranks
        rank=arraylist_rank(a)
        scale_factor = numpy.ones(self.Ntheta)
        if rank == 3:
            #sanity check. Should have shape Npsi,Ntheta,Nspecies
            if a.shape == (self.Npsi,self.Ntheta,self.Nspecies):
                return scale_factor[numpy.newaxis,:,numpy.newaxis]* numpy.einsum('kj,ijl',dds,a)
            
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 3 array badly shaped."
                sys.exit(1)
        if rank == 2:
            #sanity check. Either (Ntheta,Nspecies) or (Npsi,Ntheta) shaped
            if (a.shape == (self.Npsi,self.Ntheta)):
                return scale_factor[numpy.newaxis,:]*numpy.einsum('kj,ij',dds,a)
            elif (a.shape == (self.Ntheta,self.Nspecies)):
                return scale_factor[:,numpy.newaxis]*numpy.einsum('ij,jk',dds,a)
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 2 array badly shaped."
                sys.exit(1)
        if rank == 1:
            #sanity check. (Ntheta,) shaped
            if a.shape == (self.Ntheta,):
                return scale_factor*numpy.einsum('ij,j',dds,a)
            else:
                print "ERROR: perfect_simulation.ddpsiN: rank 1 array badly shaped."
                sys.exit(1)
        

    def internal_grid_fft(self,attrib,interval=None,use_internal_grid_4_interval=False):
        # this function simply ffts an attribute in interval, disregarding nonuniform sampling
        # this is useful for uniform grid simulations,
        # or to quantify grid-scale oscillations and other numerical errors on the internal uniform grid.
        if interval != None:
            s1 = interval[0]
            s2 = interval[1]
            if use_internal_grid_4_interval == True:
                #intepret interval as specified on internal grid
                s = self.psiN
            else:
                #intepret interval as specified on physical grid
                s = self.actual_psiN
            indices=get_index_range(s,[s1,s2],ret_range=False)
        else:
            indices = [0,self.Npsi]
        
        if type(attrib) == str:
            a=getattr(self,attrib)[indices[0]:indices[1]]
        else:
            a=attrib[indices[0]:indices[1]]
        fft=numpy.fft.rfft(a,axis=0)
        N=len(fft[:,0])
        fft=numpy.abs(fft)/float(N)
        return fft

    def export_attribute_to_h5(self,attrib,attrib_name,file_name):
        if type(attrib) == str:
            a=getattr(self,attrib)
        else:
            a=attrib
        h5file = h5py.File(file_name,'w')
            
        groupname = "Npsi"+str(self.Npsi)
        if groupname in h5file:
            print "Error: profiles already added for Npsi="+str(Npsi)
            exit(1)
        profilesgroup = h5file.create_group(groupname)
        profilesgroup.create_dataset(attrib_name,data=a)
        h5file.close()

    def export_attribute_to_text(self,attrib,file_name):
        #tested for psi, species dependent. No theta dependence
        if type(attrib) == str:
            a=getattr(self,attrib)
        else:
            a=attrib
        file = open(file_name, 'w')
        for line in a:
            try:
                for entry in line:
                    file.write(str(entry))
            except:
                file.write(str(line))
            file.write('\n')
        file.close()

    def get_index(self,ispecies,ix,L,itheta,ipsi,debug=False):
        if ispecies < 1:
            print "get_index Error: ispecies < 1"
            return None
        elif ispecies > self.Nspecies:
            print "get_index Error: ispecies > Nspecies"
            return None
        if ix < 1:
            print "get_index Error: ix < 1"
            return None
        elif ix > self.Nx:
            print "get_index Error: ix > Nx"
            return None
        if L < 0:
            print "get_index Error: L < 0"
            return None
        elif L > self.Nxi_for_x(ix)-1:
            print "get_index Error: L > Nxi_for_x(ix)-1"
            print "L: " + str(L)
            print "ix: " + str(ix)
            return None
        if itheta < 1:
            print "get_index Error: itheta < 1"
            return None
        elif itheta > self.Ntheta:
            print "get_index Error: itheta > Ntheta"
            return None
        if ipsi < 1:
            print "get_index Error: ipsi < 1"
            return None
        elif ipsi > self.Npsi:
            print "get_index Error: ipsi > Npsi"
            return None
        
        index = (ipsi-1)*self.local_matrix_size + (ispecies-1)*self.local_DKE_matrix_size + self.first_index_for_x(ix)*self.Ntheta + L*self.Ntheta + itheta - 1

        if debug:
            print "==get_index debug output==="
            print (ispecies,ix,L,itheta,ipsi)
            print index
            print "============end============"
        
        return index

    @property
    def boundary_indices(self):
        # calculate the grid-points that are at the boundary
        # return a numpy.array with indices
        if self.inputs.boundaryScheme == 3:
            # periodic BC, no boundary
            return numpy.array([])
        
        elif self.inputs.boundaryScheme == 1 or self.inputs.boundaryScheme == 2:
            if self.inputs.boundaryScheme == 1:
                # left BC only, where ipsi=1 (fortran indexing)
                ipsi = 1
            elif self.inputs.boundaryScheme == 2:
            # right BC only, where ipsi=Npsi (fortran indexing)
                ipsi = self.Npsi
           
            # size of matrix for one Npsi is given by local_matrix_size
            # thus also size of boundaries
            a = numpy.zeros(self.local_matrix_size,dtype=numpy.uint64) 
            i = 0
            for itheta in range(1,self.Ntheta+1):
                for ix in range(1,self.Nx+1):
                    for L in range(0,self.Nxi_for_x(ix)): #L is zero indexed
                        for ispecies in range(1,self.Nspecies+1):
                            a[i] = self.get_index(ispecies,ix,L,itheta,ipsi)
                            i = i + 1
            return a[0:i]
        elif self.inputs.boundaryScheme == 0:
            # WARNING: this is based on simple geometry with
            # dot(psiN) pointing upwards/downwards for different species
            a = numpy.zeros(self.local_matrix_size,dtype=numpy.uint64) 
            i = 0
            for itheta in range(1,self.Ntheta+1):
                for ix in range(1,self.Nx+1):
                    for L in range(0,self.Nxi_for_x(ix)+1): #L=0 OK
                        for ispecies in range(1,self.Nspecies+1):
                            if ispecies == self.Nspecies:
                                #electrons
                                if itheta <= self.Ntheta/2:
                                    ipsi=self.Npsi
                                else:
                                    ipsi=1
                            else:
                                if itheta <= self.Ntheta/2:
                                    ipsi=1
                                else:
                                    ipsi=self.Npsi
                            a[i] = self.get_index(ispecies,ix,L,itheta,ipsi)
                            i = i + 1
            return a[0:i]
        else:
            print "perfect_simulation: boundary_indices: ERROR, invalid boundaryScheme!"
            return None
                
    def attrib_at_psi_of_theta(self,attrib,psiN):
        index=get_index_range(self.actual_psiN,[psiN,psiN],ret_range=False)[0]
        if type(attrib) == str:
            a=getattr(self,attrib)
        else:
            a=attrib
        rank=arraylist_rank(a)
        if rank==3:
            return a[index,:,:]
        elif rank==2:
            return a[index,:]
        elif rank==1:
            return a[index]

    def attrib_at_theta_of_psi(self,attrib,theta):
        indices=get_index_range(self.theta,[theta,theta],ret_range=False,period=2*numpy.pi)
        if type(attrib) == str:
            a=getattr(self,attrib)
        else:
            a=attrib
        rank=arraylist_rank(a)
        if rank==3:  
            return a[:,indices[0],:]
        elif rank==2:
            if a.shape == (self.Npsi,self.Ntheta):
                return a[:,indices[0]]
            elif a.shape == (self.Ntheta,self.Nspecies):
                return a[indices[0],:]
            else:
                print "ERROR: perfect_simulation.attrib_at_theta_of_psi: rank 2 array badly shaped."
                sys.exit(1)
        elif rank==1:
            if a.shape == (self.Ntheta,):
                return a[indices[0]]
            else:
                print "ERROR: perfect_simulation.attrib_at_theta_of_psi: rank 1 array badly shaped."
                sys.exit(1)

    def shift_theta(self,new_zero,attrib=None):
        #shifts the theta grid to allow plotting attrib from -\pi,\pi, etc.
        #returns shifted theta array and attribute in a tuple
        #new_zero : value of theta to be the new zero index
        zero_index=get_index_range(self.theta,[new_zero,new_zero],ret_range=False,period=2*numpy.pi)[0]
        theta = numpy.concatenate((self.theta[zero_index:]-2*numpy.pi,self.theta[0:zero_index]))
        if attrib is not None:
            if type(attrib) == str:
                a=getattr(self,attrib)
            else:
                a=attrib
            rank=arraylist_rank(a)
            if rank==3:
                a =  numpy.concatenate((a[:,zero_index:,:],a[:,0:zero_index,:]),axis=1)
            elif rank==2:
                if a.shape == (self.Npsi,self.Ntheta):
                    a =  numpy.concatenate((a[:,zero_index:],a[:,0:zero_index]),axis=1)
                elif a.shape == (self.Ntheta,self.Nspecies):
                    a =  numpy.concatenate((a[zero_index:,:],a[0:zero_index,:]),axis=0)
                else:
                    print "ERROR: perfect_simulation.attrib_at_theta_of_psi: rank 2 array badly shaped."
                    sys.exit(1)
            elif rank==1:
                if a.shape == (self.Ntheta,):
                    a =  numpy.concatenate((a[zero_index:],a[0:zero_index]))
                else:
                    print "ERROR: perfect_simulation.attrib_at_theta_of_psi: rank 1 array badly shaped."
                    sys.exit(1)
            
            return (theta,a)
        else:
            return theta 
        
    def attrib_max_psi_of_theta(self,attrib,xlim=None):
        #help function to return the value of psiN that maximizes an attribute
        #for each given theta
        if xlim != None:
            #print "xlim: " + str(xlim)
            indices=get_index_range(self.actual_psiN,xlim,ret_range=True)
            psi=self.actual_psiN
        else:
            psi=self.actual_psiN
            indices=range(len(psi))
        #print [psi[indices[0]],psi[indices[-1]]]
        data=getattr(self,attrib)
        max_indices=numpy.zeros(data[0,:,:].shape)
        max_psiN=numpy.zeros(data[0,:,:].shape)
        for i_sp in range(len(data[0,0,:])):
            #for each species
            for i_th in range(len(data[0,:,i_sp])):
                #for each theta
                #print numpy.argmax(data[indices,i_th,i_sp])
                max_indices[i_th,i_sp] = numpy.argmax(data[indices,i_th,i_sp])
                max_psiN[i_th,i_sp]=psi[max_indices[i_th,i_sp]]
        #print max_indices
        #print [numpy.min(max_psiN),numpy.max(max_psiN)]
        return max_psiN
    
    @property
    def local(self):
        try:
            local=self.outputs[self.group_name+self.local_name][()]
        except (KeyError,TypeError):
            if self.inputs.makeLocalApproximation:
                local = 1
            else:
                local = -1
        if local == 1:
            return True
        elif local == -1:
            return False
        else:
            print "perfect_simulation: error: makeLocalApproximation is neither true nor false!?"

    @property
    def no_ddpsi(self):
        try:
            ddpsi=self.outputs[self.group_name+self.include_ddpsi_name][()]
        except (KeyError,TypeError):
            if self.inputs.includeddpsiTerm:
                ddpsi = 1
            else:
                ddpsi = -1
        if ddpsi == -1:
            return True
        elif ddpsi == 1:
            return False
        else:
            print "perfect_simulation: error: includeddpsiTerm is neither true nor false!?"

            
    @property
    def Npsi_sourceless_right(self):
        try:
            return self.outputs[self.group_name+self.Npsi_sourceless_right_name][()]
        except (KeyError,TypeError):
            return 0

    @property
    def Npsi_sourceless_left(self):
        try:
            return self.outputs[self.group_name+self.Npsi_sourceless_left_name][()]
        except (KeyError,TypeError):
            return 0

    def pad_sources(self,S,value=numpy.nan):
        #pads sources with nans on gridpoints where they are not calculated
        Nright = self.Npsi_sourceless_right
        Nleft = self.Npsi_sourceless_left
        Npsi = self.Npsi
        Nspecies =self.Nspecies
        padded_S = value*numpy.ones([Npsi,Nspecies])
        padded_S[Nleft:Npsi-Nright,:] = S
        return padded_S

 
            
    @property
    def Z(self):
        Z =numpy.array(self.inputs.charges)
        if Z.shape == ():
            Z.shape = (1,)
        return Z
    
    @property
    def particle_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.particle_flux_name][()]*signOfVPrimeHat

    @property
    def ion_particle_flux(self):
        return numpy.sum(self.particle_flux[:,:-1],axis=1)

    @property
    def particle_flux_fft(self):
        return self.internal_grid_fft(self.particle_flux)

    @property
    def particle_flux_before_theta_integral(self):
        VPrimeHat=self.VPrimeHat[:,numpy.newaxis,numpy.newaxis]
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.particle_flux_before_theta_integral_name][()]*signOfVPrimeHat

    @property
    def particle_flux_no_fsa(self):
        return numpy.expand_dims(numpy.fabs(self.JHat),axis=2)*self.particle_flux_before_theta_integral

    @property
    def pedestal_particle_flux_of_theta(self):
        psiN=self.pedestal_point
        return self.attrib_at_psi_of_theta(self.particle_flux_before_theta_integral,psiN)

    @property
    def core_particle_flux_of_theta(self):
        psiN=self.core_point
        return self.attrib_at_psi_of_theta(self.particle_flux_before_theta_integral,psiN)
    
    @property
    def particle_flux_test(self):
        #tests the external FSA routine by calculating the particle flux from the non-integrated value.
        return numpy.expand_dims(self.VPrimeHat,axis=1)*self.FSA("particle_flux_before_theta_integral",jacobian=False)

    @property
    def particle_flux_psi_unit_vector(self):
        #divide \nabla \psi with |\nabla \psi| =|RB_p| (assumes BBar, RBar positive)
        return numpy.expand_dims(self.VPrimeHat,axis=1)*self.FSA(self.particle_flux_before_theta_integral/(numpy.fabs(numpy.expand_dims(self.RHat*self.Bp,axis=2))),jacobian=False)

    @property
    def particle_flux_psi_unit_vector_no_fsa(self):
        #divide \nabla \psi with |\nabla \psi| =|RB_p| (assumes BBar, RBar positive)
        return self.particle_flux_no_fsa/(numpy.fabs(numpy.expand_dims(self.RHat*self.Bp,axis=2)))
        
    @property
    def momentum_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.momentum_flux_name][()]*signOfVPrimeHat

    @property
    def momentum_flux_before_theta_integral(self):
        VPrimeHat=self.VPrimeHat[:,numpy.newaxis,numpy.newaxis]
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.momentum_flux_before_theta_integral_name][()]*signOfVPrimeHat

    
    @property
    def momentum_flux_test(self):
        return numpy.expand_dims(self.VPrimeHat,axis=1)*self.FSA("momentum_flux_before_theta_integral",jacobian=False)

    
    @property
    def momentum_flux_psi_unit_vector(self):
        #divide \nabla \psi with |\nabla \psi| =|RB_p| (assumes BBar, RBar positive)
        return numpy.expand_dims(self.VPrimeHat,axis=1)*self.FSA(self.momentum_flux_before_theta_integral/(numpy.fabs(numpy.expand_dims(self.RHat*self.Bp,axis=2))),jacobian=False)

    
    @property
    def m_momentum_flux(self):
        return self.momentum_flux*self.masses

    @property
    def m_momentum_flux_psi_unit_vector(self):
        return self.momentum_flux_psi_unit_vector*self.masses

    @property
    def sum_m_momentum_flux(self):
        return numpy.sum(self.masses*self.momentum_flux,axis=1)

    @property
    def sum_m_momentum_flux_psi_unit_vector(self):
        return numpy.sum(self.m_momentum_flux_psi_unit_vector,axis=1)

    def Prandtl_proxy_test(self,ion_index = 0):
        #should be identical to Prandtl_proxy
        psiN_point = self.pedestal_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point],ret_range=False)[0]
        return self.momentum_flux*((1/self.conductive_heat_flux[:,ion_index])*(self.nHat[:,ion_index]*self.dTHatdpsiN[:,ion_index])/(self.Delta*self.FSA_toroidal_flow[:,ion_index]*self.dnHatdpsiN[:,ion_index]))[psiN_index]

    def Prandtl_proxy_dVtdpsiN(self,ion_index = 0):
        #more correct Prandtl proxy
        # use this!
        psiN_index=self.mid_pedestal_point_index
        return self.masses*self.momentum_flux*((1/self.conductive_heat_flux[:,ion_index])*(self.nHat[:,ion_index]*self.dTHatdpsiN[:,ion_index])/(self.Delta*self.masses[ion_index]*self.ddpsiN_n_FSA_toroidal_flow[:,ion_index]))[psiN_index]

    @property
    def Prandtli_proxy_dVtdpsiN(self):
        return self.Prandtl_proxy_dVtdpsiN(ion_index = 0)
    
    @property
    def heat_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.heat_flux_name][()]*signOfVPrimeHat

    @property
    def ion_heat_flux(self):
        return numpy.sum(self.heat_flux[:,:-1],axis=1)
    
    @property
    def heat_flux_before_theta_integral(self):
        VPrimeHat=self.VPrimeHat[:,numpy.newaxis,numpy.newaxis]
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return self.outputs[self.group_name+self.heat_flux_before_theta_integral_name][()]*signOfVPrimeHat

    @property
    def heat_flux_psi_unit_vector(self):
        #divide \nabla \psi with |\nabla \psi| =|RB_p| (assumes BBar, RBar positive)
        return numpy.expand_dims(self.VPrimeHat,axis=1)*self.FSA(self.heat_flux_before_theta_integral/(numpy.fabs(numpy.expand_dims(self.RHat*self.Bp,axis=2))),jacobian=False)

    @property
    def conductive_heat_flux(self):        
        return self.heat_flux-(5.0/2.0)*self.THat*self.particle_flux

    @property
    def ion_conductive_heat_flux(self):
        return numpy.sum(self.conductive_heat_flux[:,:-1],axis=1)
    
    @property
    def conductive_heat_flux_fft(self):
        return self.internal_grid_fft(self.conductive_heat_flux)

    @property
    def n_ped(self):
        #species wise LCFS density
        return self.nHat[self.pedestal_start_stop_indices[0]]
        # psiMid = (self.pedestal_start_stop[0] + self.pedestal_start_stop[1])/2.0
        # i = self.mid_pedestal_point_index
        # if psiMid > self.actual_psiN[i]:
        #     d1 = (psiMid - self.actual_psiN[i])
        #     d2 = (self.actual_psiN[i+1] - psiMid)
        #     d = self.actual_psiN[i+1] - self.actual_psiN[i]
        #     n0 = self.nHat[i]*(d-d1)/d + self.nHat[i+1]*(d-d2)/d
        #     ddpsiN_n0 = self.dnHatdpsiN[i]*(d-d1)/d + self.dnHatdpsiN[i+1]*(d-d2)/d
        # if psiMid < self.actual_psiN[i]:
        #     d1 = (self.actual_psiN[i] - psiMid )
        #     d2 = (psiMid - self.actual_psiN[i-1])
        #     d = self.actual_psiN[i] - self.actual_psiN[i-1]
        #     n0 = self.nHat[i]*(d-d1)/d + self.nHat[i-1]*(d-d2)/d
        #     ddpsiN_n0 = self.dnHatdpsiN[i]*(d-d1)/d + self.dnHatdpsiN[i+1]*(d-d2)/d
        # print n0
        # print ddpsiN_n0
        # ddpsiN_n0 = -12.251148545176108
        # n0 = numpy.array([0.38])
        # w = self.pedestal_width_psiN
        # return n0 - ddpsiN_n0 * w/2
    
    @property
    def n_LCFS(self):
        #species wise LCFS density
        return self.nHat[self.pedestal_start_stop_indices[1]]
        #
        #i = self.mid_pedestal_point_index
        #n0 = self.nHat[i]
        #ddpsiN_n0 = self.dnHatdpsiN[i]
        # psiMid = (self.pedestal_start_stop[0] + self.pedestal_start_stop[1])/2.0
        # i = self.mid_pedestal_point_index
        # if psiMid > self.actual_psiN[i]:
        #     d1 = (psiMid - self.actual_psiN[i])
        #     d2 = (self.actual_psiN[i+1] - psiMid)
        #     d = self.actual_psiN[i+1] - self.actual_psiN[i]
        #     n0 = self.nHat[i]*(d-d1)/d + self.nHat[i+1]*(d-d2)/d
        #     ddpsiN_n0 = self.dnHatdpsiN[i]*(d-d1)/d + self.dnHatdpsiN[i+1]*(d-d2)/d
        # if psiMid < self.actual_psiN[i]:
        #     d1 = (self.actual_psiN[i] - psiMid )
        #     d2 = (psiMid - self.actual_psiN[i-1])
        #     d = self.actual_psiN[i] - self.actual_psiN[i-1]
        #     n0 = self.nHat[i]*(d-d1)/d + self.nHat[i-1]*(d-d2)/d
        #     ddpsiN_n0 = self.dnHatdpsiN[i]*(d-d1)/d + self.dnHatdpsiN[i+1]*(d-d2)/d
        # n0 = numpy.array([0.38])
        # ddpsiN_n0 = -12.251148545176108
        # w = self.pedestal_width_psiN
        # return n0 + ddpsiN_n0 * w/2
    
    @property
    def q_max(self):
        #species wise maximum in q
        q=self.conductive_heat_flux
        imax=numpy.argmax(numpy.fabs(q),axis=0)
        return numpy.array([q[im,i] for i,im in enumerate(imax)])

    @property
    def q_max_psiN(self):
        #species wise psiN location of maximum in q
        q=self.conductive_heat_flux
        imax=numpy.argmax(numpy.fabs(q),axis=0)
        return self.actual_psiN[imax]

    @property
    def Q_PPV(self):
        return [PPV(self.heat_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 

        
    @property
    def Q_PPV_psiN(self):
        return [PPV_x(self.heat_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)]

    @property
    def Qi_PPV(self):
        return PPV(self.ion_heat_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)

    @property
    def Qi_PPV_psiN(self):
        return PPV_x(self.ion_heat_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)

    @property
    def q_PPV(self):
        return [PPV(self.conductive_heat_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 


    @property
    def q_PPV_psiN(self):
        return [PPV_x(self.conductive_heat_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)]

    @property
    def q_mid_ped(self):
        i=self.mid_pedestal_point_index
        return self.conductive_heat_flux[i]

    @property
    def qi_PPV(self):
        return PPV(self.ion_conductive_heat_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)

    @property
    def qi_PPV_psiN(self):
        return PPV_x(self.ion_conductive_heat_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)

    @property
    def qi_mid_ped(self):
        i=self.mid_pedestal_point_index
        return self.ion_conductive_heat_flux[i]

    @property
    def Pi_PPV(self):
        return [PPV(self.momentum_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 


    @property
    def Pi_PPV_psiN(self):
        return [PPV_x(self.momentum_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)]

    @property
    def Prandtl_PPV(self):
        return [PPV(self.Prandtli_proxy_dVtdpsiN[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 


    @property
    def Prandtl_PPV_psiN(self):
        return [PPV_x(self.Prandtli_proxy_dVtdpsiN[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)]

    @property
    def Gamma_PPV(self):
        return [PPV(self.particle_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 

        
    @property
    def Gamma_PPV_psiN(self):
        return [PPV_x(self.particle_flux[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)] 

    @property
    def Gammai_PPV(self):
        return PPV(self.ion_particle_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)

    @property
    def Gammai_PPV_psiN(self):
        return PPV_x(self.ion_particle_flux,self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4)
    
    @property
    def Sh_PPV(self):
        factor = 0.6
        pedestal_width = (self.pedestal_start_stop_psiN3[1]-self.pedestal_start_stop_psiN3[0])
        shift = factor*pedestal_width
        return [PPV(self.heat_source_over_m2[:,ispecies],self.psiN3,[self.pedestal_start_stop_psiN3[0]-shift,self.pedestal_start_stop_psiN3[1]],tol=1e-10,order=4,retmode="absy") for ispecies in range(self.Nspecies)]

    @property
    def Sh_PPV_psiN(self):
        factor = 0.6
        pedestal_width = (self.pedestal_start_stop_psiN3[1]-self.pedestal_start_stop_psiN3[0])
        shift = factor*pedestal_width
        return [PPV_x(self.heat_source_over_m2[:,ispecies],self.psiN3,[self.pedestal_start_stop_psiN3[0]-shift,self.pedestal_start_stop_psiN3[1]],tol=1e-10,order=4) for ispecies in range(self.Nspecies)]
        #return [PPV_x(self.heat_source_over_m2[:,ispecies],self.psiN3,self.pedestal_start_stop_psiN3,tol=1e-10,order=4) for ispecies in range(self.Nspecies)]

    @property
    def special_Pi_PPV(self):
        factor = 0.6
        pedestal_width = (self.pedestal_start_stop_psiN3[1]-self.pedestal_start_stop_psiN3[0])
        shift = factor*pedestal_width
        return [PPV(self.momentum_flux[:,ispecies],self.psiN3,[self.pedestal_start_stop_psiN3[0]-shift,self.pedestal_start_stop_psiN3[1]],tol=1e-10,order=4,retmode="absy") for ispecies in range(self.Nspecies)]

    @property
    def special_Pi_PPV_psiN(self):
        factor = 0.6
        pedestal_width = (self.pedestal_start_stop_psiN3[1]-self.pedestal_start_stop_psiN3[0])
        shift = factor*pedestal_width
        return [PPV_x(self.momentum_flux[:,ispecies],self.psiN3,[self.pedestal_start_stop_psiN3[0]-shift,self.pedestal_start_stop_psiN3[1]],tol=1e-10,order=4) for ispecies in range(self.Nspecies)]
    
    @property
    def conductive_heat_flux_psi_unit_vector(self):        
        return self.heat_flux_psi_unit_vector-(5.0/2.0)*self.THat*self.particle_flux_psi_unit_vector

    @property
    def particle_flux_over_nPed(self):
        # not really nPed
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.particle_flux/npeds
    
    
    @property
    def momentum_flux_over_nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.momentum_flux/npeds

    @property
    def m_momentum_flux_over_nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.m_momentum_flux/npeds

    @property
    def heat_flux_over_nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.heat_flux/npeds

    @property
    def conductive_heat_flux_over_nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.conductive_heat_flux/npeds
    
    @property
    def jHat(self):
        return numpy.sum(self.Z*self.particle_flux,axis=1)

    @property
    def sum_QHat(self):
        return numpy.sum(self.heat_flux,axis=1)

    @property
    def mPiHat_source_jHat(self):
        return (self.psiAHat/self.Delta)*self.jHat

    @property
    def minus_mPiHat_source_jHat(self):
        return -self.mPiHat_source_jHat

    @property
    def mPiHat_source_source(self):
        return numpy.sum(self.masses*self.PiHat_source_source,axis=1)

    @property
    def mPiHat_jHat_source(self):
        return (self.Delta/self.psiAHat)*numpy.sum(self.masses*self.PiHat_source_source,axis=1)

    @property
    def electron_particle_flux(self):
        return self.particle_flux[:,-1]
    
    @property
    def mPiHat_source_source_filtered(self):
        a=self.mPiHat_source_source
        gaussian_filter=scipy.ndimage.filters.gaussian_filter1d
        sigma = 3
        order = 0
        a[0]=0.0
        a[-1]=0.0
        filtered=gaussian_filter(a, sigma, axis=0, order=order, output=None, mode='constant', cval=0.0, truncate=4.0)
        return filtered

    @property
    def nosum_mPiHat_source_source(self):
        return self.masses*self.PiHat_source_source

    @property
    def nosum_mPiHat_source_source_filtered(self):
        a=self.nosum_mPiHat_source_source
        gaussian_filter=scipy.ndimage.filters.gaussian_filter1d
        sigma = 3
        order = 0
        a[0]=0.0
        a[-1]=0.0
        filtered=gaussian_filter(a, sigma, axis=0, order=order, output=None, mode='constant', cval=0.0, truncate=4.0)
        return filtered
        
        
    
    #for backwards compatibility, for now.
    @property
    def ambipolarity(self):
        return self.jHat

    @property
    def conductive_heat_flux_sum(self):
        return numpy.sum(self.conductive_heat_flux,axis=1)

    @property
    def conductive_heat_flux_sum_over_ne(self):
        return self.conductive_heat_flux_sum/(self.nHat[:,-1]**2)

    
    @property
    def momentum_flux_sum(self):
        return numpy.sum(self.masses*self.momentum_flux,axis=1)

    @property
    def VPrimeHat(self):
        return self.outputs[self.group_name+self.VPrimeHat_name][()]

    @property
    def heat_source(self):
        try:
            S = self.outputs[self.group_name+self.heat_source_name][()]
        except KeyError:
            S=0
        return self.pad_sources(S) + self.constant_heat_source + self.species_indep_source_heat_source

    @property
    def momentum_source(self):
        try:
            S = self.outputs[self.group_name+self.momentum_source_name][()]
        except KeyError:
            S=0
        return self.pad_sources(S) + self.constant_momentum_source + self.no_charge_source_momentum_source + self.species_indep_source_momentum_source
    
    @property
    def particle_source(self):
        # Sources may be nonexisting at the boundary
        try:
            S = self.outputs[self.group_name+self.particle_source_name][()]
        except KeyError:
            S=0
        return self.pad_sources(S) + self.constant_particle_source + self.species_indep_source_particle_source
    

    @property
    def heat_source_over_m2(self):
        return self.heat_source/(self.masses**2)

    @property
    def particle_source_over_m2(self):
        return self.particle_source/(self.masses**2)

    @property
    def particle_source_over_m2nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.particle_source_over_m2/npeds

    @property
    def heat_source_over_m2nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        return self.heat_source_over_m2/npeds
    
    @property
    def kinetic_energy_source(self):
        VPrimeHat = numpy.fabs(self.VPrimeHat)
        return self.heat_source + (2*self.omega)/(3*numpy.sqrt(numpy.pi)*self.psiAHat) * (self.Z*self.masses**2 *numpy.expand_dims(self.dPhiHatdpsiN,axis=1) *self.particle_flux)/(self.THat**(5./2.) * numpy.expand_dims(VPrimeHat,axis=1))

    @property
    def kinetic_energy_source_addition(self):
        VPrimeHat = numpy.fabs(self.VPrimeHat)
        return (2*self.omega)/(3*numpy.sqrt(numpy.pi)*self.psiAHat) * (self.Z*self.masses**2 *numpy.expand_dims(self.dPhiHatdpsiN,axis=1) *self.particle_flux)/(self.THat**(5./2.) * numpy.expand_dims(VPrimeHat,axis=1))
    
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
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        prefact=1#(self.psiAHatArray[:,numpy.newaxis]/self.psiAHat)
        #print self.psiAHat
        return prefact*2*(numpy.sqrt(numpy.pi)/2)*self.THat**(3./2.)*self.psiAHat*numpy.expand_dims(VPrimeHat,axis=1)/(self.Delta)*self.particle_source/(self.masses**2)

    @property
    def QHat_source(self):
        #for comparrison with dQHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        return -self.THat**(5./2.)*self.psiAHat*numpy.expand_dims(VPrimeHat,axis=1)*numpy.sqrt(numpy.pi)/(2*self.Delta*self.masses**2)*(3*self.heat_source) -((2*self.omega/self.Delta)*self.Z*numpy.expand_dims(self.dPhiHatdpsiN,axis=1))*self.particle_flux*self.masses**(1.0/2.0)

    @property
    def QHat_source_source(self):
        #for comparrison with dQHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        return -self.THat**(5./2.)*self.psiAHat*numpy.expand_dims(VPrimeHat,axis=1)*numpy.sqrt(numpy.pi)/(2*self.Delta*self.masses**2)*(3*self.heat_source)


    @property
    def sum_QHat_source(self):
        #for comparrison with \sum dQHat/dPsiN
        return numpy.sum(self.QHat_source,axis=1)

    
    @property
    def sum_QHat_source_minus_ddpsiN_sum_QHat(self):
        return self.sum_QHat_source - self.ddpsiN_sum_QHat

    @property
    def sum_QHat_source_source(self):
        #for comparrison with \sum dQHat/dPsiN
        return numpy.sum(self.QHat_source_source,axis=1)

    @property
    def sum_QHat_source_source_minus_ddpsiN_sum_QHat(self):
        return self.sum_QHat_source_source - self.ddpsiN_sum_QHat

        

    @property
    def QHat_source_addition(self):
        #for comparrison with dQHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        return -((2*self.omega/self.Delta)*self.Z*numpy.expand_dims(self.dPhiHatdpsiN,axis=1))*self.particle_flux*self.masses**(1.0/2.0)
    
    @property
    def qHat_source(self):
        #for comparrison with dqHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        return -self.THat**(5./2.)*self.psiAHat*numpy.expand_dims(VPrimeHat,axis=1)*numpy.sqrt(numpy.pi)/(2*self.Delta*self.masses**2)*(5*self.particle_source + 3*self.heat_source) -((2*self.omega/self.Delta)*self.Z*numpy.expand_dims(self.dPhiHatdpsiN,axis=1) + (5./2.)*self.dTHatdpsiN)*self.particle_flux

    @property
    def heat_source_Parra(self):
        #do not sum over electrons
        #assumes electron is the last index
        VPrimeHat=numpy.fabs(self.VPrimeHat) 
        return (numpy.sqrt(numpy.pi)/4) * VPrimeHat * numpy.sum((self.THat**(5./2)/(self.Z*self.masses) * numpy.expand_dims(self.FSA(self.IHat**2/(self.BHat*self.RHat)**2),axis=1) * self.heat_source)[:,:-1],axis=1)

    @property
    def kinetic_energy_source_Parra(self):
        #do not sum over electrons
        #assumes electron is the last index
        VPrimeHat=numpy.fabs(self.VPrimeHat) 
        return (numpy.sqrt(numpy.pi)/4) * VPrimeHat * numpy.sum((self.THat**(5./2)/(self.Z*self.masses) * numpy.expand_dims(self.FSA(self.IHat**2/(self.BHat*self.RHat)**2),axis=1) * self.kinetic_energy_source)[:,:-1],axis=1)
    
    @property
    def kinetic_energy_source_Parra_addition(self):
        #do not sum over electrons
        #assumes electron is the last index
        VPrimeHat=numpy.fabs(self.VPrimeHat) 
        return (numpy.sqrt(numpy.pi)/4) * VPrimeHat * numpy.sum((self.THat**(5./2)/(self.Z*self.masses) * numpy.expand_dims(self.FSA(self.IHat**2/(self.BHat*self.RHat)**2),axis=1) * self.kinetic_energy_source_addition)[:,:-1],axis=1)
    
    @property
    def PiHat_source_source(self):
        #for comparrison with dPiHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        S = self.momentum_source
        S_nan_indices=numpy.isnan(S)
        if S_nan_indices.any():
            # nan sources
            print "Some momentum sources are NaN. Interpreting NaN as 0."
            S[S_nan_indices] = 0.0

        # NOTE: assumes that Theta_m = 1
        Theta_m = 1        
        PiHat_source = (self.psiAHat/self.Delta)*numpy.sqrt(numpy.pi)/(2*self.masses**(5./2.))*numpy.expand_dims(self.FSA(Theta_m*self.IHat/self.BHat)*VPrimeHat,axis=1)*self.THat**2*S
        return PiHat_source

    @property
    def PiHat_source_source_filtered(self):
        a=self.PiHat_source_source
        gaussian_filter=scipy.ndimage.filters.gaussian_filter1d
        sigma = 3
        order = 0
        a[0,:]=0.0
        a[-1,:]=0.0
        filtered=gaussian_filter(a, sigma, axis=0, order=order, output=None, mode='constant', cval=0.0, truncate=4.0)
        return filtered
        return PiHat_source
    
    
    @property
    def PiHat_source(self):
        #for comparrison with dPiHat/dPsiN
        PiHat_source = (self.psiAHat/self.Delta)*(self.Z/self.masses)*self.particle_flux + self.PiHat_source_source
        return PiHat_source

    @property
    def noChargeSource(self):
        try:
            return self.outputs[self.group_name+self.no_charge_source_name][()]    
        except KeyError:
            return 0

    @property
    def no_charge_source_momentum_source_species_dependence(self):
        try:
            return self.outputs[self.group_name+self.no_charge_source_momentum_source_species_dependence_name][()]
        except KeyError:
            return numpy.zeros(self.Nspecies)
        
    @property
    def no_charge_source_momentum_source(self):
        try:
            S=self.outputs[self.group_name+self.no_charge_source_momentum_source_name][()]
        except KeyError:
            S=numpy.nan*numpy.ones([self.Npsi,self.Nspecies])
        S_nan_indices=numpy.isnan(S)
        if S_nan_indices.any():
            # nan sources
            print "Some momentum sources are NaN. Interpreting NaN as 0."
            S[S_nan_indices] = 0.0
        return S

    @property
    def momentum_source_over_m52nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.nHat[psiN_index,i] for i in range(self.num_species)]
        S = self.momentum_source/(self.masses**5./2.)/npeds
        return S
    
    @property
    def no_charge_source_momentum_source_filtered(self):
        a=self.no_charge_source_momentum_source
        gaussian_filter=scipy.ndimage.filters.gaussian_filter1d
        sigma = 3
        order = 0
        a[0,:]=0.0
        a[-1,:]=0.0
        filtered=gaussian_filter(a, sigma, axis=0, order=order, output=None, mode='constant', cval=0.0, truncate=4.0)
        return filtered
        
    @property
    def no_charge_source_particle_source(self):
        try:
            return self.outputs[self.group_name+self.no_charge_source_particle_source_name][()]
        except KeyError:
            return numpy.nan*numpy.ones([self.Npsi,self.Nspecies])

    @property
    def constant_particle_source(self):
        try:
            return self.outputs[self.group_name+self.constant_source_particle_source_name][()]
        except KeyError:
            return 0

    @property
    def constant_heat_source(self):
        try:
            return self.outputs[self.group_name+self.constant_source_heat_source_name][()]
        except KeyError:
            return 0

    @property
    def constant_momentum_source(self):
        try:
            return self.outputs[self.group_name+self.constant_source_momentum_source_name][()]
        except KeyError:
            return 0

    @property
    def species_indep_source_particle_source(self):
        try:
            S=self.outputs[self.group_name+self.species_indep_source_particle_source_name][()]
        except KeyError:
            S=numpy.nan*numpy.ones([self.Npsi,self.Nspecies])
        S_nan_indices=numpy.isnan(S)
        if S_nan_indices.any():
            # nan sources
            print "Some momentum sources are NaN. Interpreting NaN as 0."
            S[S_nan_indices] = 0.0
        return S

    @property
    def species_indep_source_heat_source(self):
        try:
            S=self.outputs[self.group_name+self.species_indep_source_heat_source_name][()]
        except KeyError:
            S=numpy.nan*numpy.ones([self.Npsi,self.Nspecies])
        S_nan_indices=numpy.isnan(S)
        if S_nan_indices.any():
            # nan sources
            print "Some momentum sources are NaN. Interpreting NaN as 0."
            S[S_nan_indices] = 0.0
        return S

    @property
    def species_indep_source_momentum_source(self):
        try:
            S=self.outputs[self.group_name+self.species_indep_source_momentum_source_name][()]
        except KeyError:
            S=numpy.nan*numpy.ones([self.Npsi,self.Nspecies])
        S_nan_indices=numpy.isnan(S)
        if S_nan_indices.any():
            # nan sources
            print "Some momentum sources are NaN. Interpreting NaN as 0."
            S[S_nan_indices] = 0.0
        return S
        
    @property
    def FSAFlow(self):
        return self.outputs[self.group_name+self.FSAFlow_name][()]

    @property
    def parallel_current(self):
        return numpy.sum(self.Z*self.FSAFlow,axis=1)
    
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
    def lboundary_flow(self):
        psiN_point=self.actual_psiN[0]
        return self.attrib_at_psi_of_theta("flow",psiN_point)

    @property
    def rboundary_flow(self):
        psiN_point=self.actual_psiN[-1]
        return self.attrib_at_psi_of_theta("flow",psiN_point)
    
    
    @property
    def flow_minus_FSAFlow(self):
        return self.flow-numpy.expand_dims(self.FSAFlow,axis=1)

    
    @property
    def flow_over_vT(self):
        return self.Delta*nstack(numpy.sqrt(self.masses/self.THat),axis=1,n=len(self.theta))*self.flow

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
    def toroidal_flow(self):
        return self.toroidal_flow_test
        #return self.outputs[self.group_name+self.toroidal_flow_name][()]

    @property
    def toroidal_flow_at_psi_of_theta(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta("toroidal_flow",psiN_point)

    @property
    def toroidal_flow_outboard(self):
        return self.attrib_at_theta_of_psi("toroidal_flow",0)

    @property
    def toroidal_flow_inboard(self):
        return self.attrib_at_theta_of_psi("toroidal_flow",numpy.pi)

    @property
    def FSA_flow(self):
        return self.FSAFlow
    
    @property
    def FSA_toroidal_flow(self):
        return self.outputs[self.group_name+self.FSA_toroidal_flow_name][()]

    @property
    def FSA_poloidal_flow(self):
        return self.outputs[self.group_name+self.FSA_poloidal_flow_name][()]

    @property
    def ddpsiN_n_FSA_toroidal_flow(self):
        return self.ddpsiN(self.nHat*self.FSA_toroidal_flow)

    @property
    def FSA_toroidal_flow_test(self):
        return self.FSA("toroidal_flow",jacobian=True)

    @property
    def toroidal_flow_over_vT(self):
        return self.Delta*nstack(numpy.sqrt(self.masses/self.THat),axis=1,n=len(self.theta))*self.toroidal_flow

    @property
    def FSA_toroidal_flow_over_vT(self):
        return self.Delta*self.FSA_toroidal_flow/numpy.sqrt(self.THat/self.masses)
        
    @property
    def toroidal_flow_difference(self):
        if self.num_species>1:
            main_index=0
            impurity_index=1
            return self.toroidal_flow[:,:,main_index] - self.toroidal_flow[:,:,impurity_index]
        else:
            print "perfect_simulation: toroidal_flow_difference: warning: no impurity species. Difference will be just the main ion toroidal flow"
            return self.toroidal_flow[:,:,main_index]

    @property
    def FSA_toroidal_flow_shear(self):
        # assumes local Miller
        dRdr = self.inputs.Miller_dRdr
        epsil = self.epsilon
        q = self.q
        gamma = (epsil/q) * ((self.FSA_toroidal_flow/numpy.expand_dims(self.FSARHat,axis=1)) * dRdr - numpy.expand_dims((self.FSA(self.RHat *self.Bp)/self.psiAHat),axis=1) *self.ddpsiN(self.FSA_toroidal_flow))
        return gamma
        
    @property
    def poloidal_flow(self):
        #next line commented out since Bp seem to have the wrong sign in the code, and since we want to use higher order stencil:
        #return self.outputs[self.group_name+self.poloidal_flow_name][()]
        return self.poloidal_flow_test
        
    @property
    def poloidal_flow_at_psi_of_theta(self):
        psiN_point = self.pedestal_point
        return self.attrib_at_psi_of_theta("poloidal_flow",psiN_point)

    @property
    def poloidal_flow_outboard(self):
        return self.attrib_at_theta_of_psi("poloidal_flow",0)

    @property
    def poloidal_flow_inboard(self):
        return self.attrib_at_theta_of_psi("poloidal_flow",numpy.pi)

    @property
    def poloidal_flow_inboard_ExB(self):
        return self.attrib_at_theta_of_psi("poloidal_flow_ExB",numpy.pi)

    @property
    def poloidal_flow_inboard_gradP(self):
        return self.attrib_at_theta_of_psi("poloidal_flow_gradP",numpy.pi)

    
    @property
    def poloidal_flow_inboard_inputs_parallel(self):
        return self.attrib_at_theta_of_psi("poloidal_flow_inputs_parallel",numpy.pi)

    @property
    def poloidal_flow_to_k_poloidal_factor(self):
        return 2*self.psiAHat*self.Z*self.FSABHat2[:,numpy.newaxis,numpy.newaxis]/(numpy.expand_dims(self.Bp*self.IHat,axis=2)*numpy.expand_dims(self.dTHatdpsiN,axis=1))

    @property
    def toroidal_flow_to_k_toroidal_factor(self):
        return 2*self.psiAHat*self.Z*self.FSABHat2[:,numpy.newaxis,numpy.newaxis]/(numpy.expand_dims(self.Bt*self.IHat,axis=2)*numpy.expand_dims(self.dTHatdpsiN,axis=1))

    @property
    def k_toroidal(self):
        return self.toroidal_flow*self.toroidal_flow_to_k_toroidal_factor
    
    @property
    def k_poloidal(self):
        return self.poloidal_flow*self.poloidal_flow_to_k_poloidal_factor

    @property
    def k_poloidal_at_psi_of_theta(self):
        psiN_point = self.pedestal_point
        return self.attrib_at_psi_of_theta("k_poloidal",psiN_point)

    @property
    def k_poloidal_outboard(self):
        return self.attrib_at_theta_of_psi("k_poloidal",0)

    @property
    def k_poloidal_outboard_ExB(self):
        return self.attrib_at_theta_of_psi("k_poloidal_ExB",0)

    @property
    def k_poloidal_outboard_gradP(self):
        return self.attrib_at_theta_of_psi("k_poloidal_gradP",0)

    @property
    def k_poloidal_outboard_gradP_test(self):
        return self.attrib_at_theta_of_psi("k_poloidal_gradP_test",0)

    
    @property
    def k_poloidal_outboard_inputs_parallel(self):
        return self.attrib_at_theta_of_psi("k_poloidal_inputs_parallel",0)

    @property
    def k_poloidal_outboard_inputs(self):
        return self.attrib_at_theta_of_psi("k_poloidal_inputs",0)

    @property
    def k_poloidal_outboard_parallel(self):
        return self.attrib_at_theta_of_psi("k_poloidal_parallel",0)
    
    @property
    def k_poloidal_inboard(self):
        return self.attrib_at_theta_of_psi("k_poloidal",numpy.pi)
    
    @property
    def k_poloidal_inboard_ExB(self):
        return self.attrib_at_theta_of_psi("k_poloidal_ExB",numpy.pi)

    @property
    def k_poloidal_inboard_gradP(self):
        return self.attrib_at_theta_of_psi("k_poloidal_gradP",numpy.pi)

    
    @property
    def k_poloidal_inboard_gradP_test(self):
        return self.attrib_at_theta_of_psi("k_poloidal_gradP_test",numpy.pi)

    
    @property
    def k_poloidal_inboard_inputs_parallel(self):
        return self.attrib_at_theta_of_psi("k_poloidal_inputs_parallel",numpy.pi)

    @property
    def k_poloidal_inboard_inputs(self):
        return self.attrib_at_theta_of_psi("k_poloidal_inputs",numpy.pi)

    @property
    def k_poloidal_inboard_parallel(self):
        return self.attrib_at_theta_of_psi("k_poloidal_parallel",numpy.pi)

    @property
    def magnetization_perturbation_inboard(self):
        return self.attrib_at_theta_of_psi("magnetization_perturbation",numpy.pi)

    @property
    def magnetization_perturbation_outboard(self):
        return self.attrib_at_theta_of_psi("magnetization_perturbation",0.0)

    
    @property
    def magnetization_flow_perturbation_inboard(self):
        return self.attrib_at_theta_of_psi("magnetization_flow_perturbation",numpy.pi)

    @property
    def magnetization_flow_perturbation_outboard(self):
        return self.attrib_at_theta_of_psi("magnetization_flow_perturbation",0.0)

    @property
    def perpendicular_flow_inboard(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow",numpy.pi)

    @property
    def perpendicular_flow_outboard(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow",0)

    @property
    def perpendicular_flow_inboard_ExB(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_ExB",numpy.pi)

    @property
    def perpendicular_flow_outboard_ExB(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_ExB",0)

    @property
    def perpendicular_flow_inboard_gradP(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_gradP",numpy.pi)

    @property
    def perpendicular_flow_outboard_gradP(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_gradP",0)

    @property
    def perpendicular_flow_inboard_inputs(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs",numpy.pi)

    @property
    def perpendicular_flow_outboard_inputs(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs",0)

    @property
    def perpendicular_flow_inboard_inputs_ExB(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs_ExB",numpy.pi)

    @property
    def perpendicular_flow_outboard_inputs_ExB(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs_ExB",0)

    @property
    def perpendicular_flow_inboard_inputs_diamagnetic(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs_diamagnetic",numpy.pi)

    @property
    def perpendicular_flow_outboard_inputs_diamagnetic(self):
        return self.attrib_at_theta_of_psi("perpendicular_flow_inputs_diamagnetic",0)

    
    
    
    @property
    def poloidal_flow_over_vT(self):
        return self.Delta*nstack(numpy.sqrt(self.masses/self.THat),axis=1,n=len(self.theta))*self.poloidal_flow

    @property
    def poloidal_flow_difference(self):
        if self.num_species>1:
            main_index=0
            impurity_index=1
            return self.poloidal_flow[:,:,main_index] - self.poloidal_flow[:,:,impurity_index]
        else:
            print "perfect_simulation: poloidal_flow_difference: warning: no impurity species. Difference will be just the main ion poloidal flow"
            return self.poloidal_flow[:,:,main_index]

    
    @property
    def magnetization_flow_perturbation(self):
        try:
            return self.outputs[self.group_name+self.magnetization_flow_perturbation_name][()]
        except KeyError:
            return self.outputs[self.group_name+self.p_perp_term_in_Vp_name][()]
        #return self.magnetization_flow_perturbation_test
        
    @property
    def magnetization_perturbation(self):
        try:
            return self.outputs[self.group_name+self.magnetization_perturbation_name][()]
        except KeyError:
            return self.outputs[self.group_name+self.p_perp_term_in_Vp_underived_name][()]
        #return numpy.expand_dims(3*self.nHat*self.THat/(self.Delta*self.masses),axis=1)*self.pressure_perturbation - self.outputs[self.group_name+self.magnetization_perturbation_name][()]

    @property
    def flow_max_psi_of_theta(self):
        range=[self.core_point,self.pedestal_start_stop[-1]]
        return self.attrib_max_psi_of_theta("flow",range)

    @property
    def flow_at_psi_of_theta(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta("flow",psiN_point)

    
    @property
    def kPar(self):
        return self.outputs[self.group_name+self.kPar_name][()]

    @property
    def kPar_minus_FSAkPar(self):
        return self.kPar-numpy.expand_dims(self.FSAkPar,axis=1)

    
    @property
    def kPar_max_psi_of_theta(self):
        range=[self.core_point,self.pedestal_start_stop[-1]]
        return self.attrib_max_psi_of_theta("kPar",range)

    @property
    def kPar_at_psi_of_theta(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta("kPar",psiN_point)

    

    @property
    def FSABJPar(self):
        return numpy.sum(self.Z*(self.nHat*self.FSABFlow),axis=1)

    @property
    def potential_perturbation(self):
        prefactor=numpy.sum((self.Z**2)*(self.nHat/self.THat),axis=1)**(-1)*self.Delta/(2*self.omega)
        prefactor=prefactor[:,numpy.newaxis]
        return prefactor*self.charge_perturbation

    @property
    def charge_perturbation(self):
        Zn=numpy.expand_dims(self.Z*self.nHat, axis=1)
        return numpy.sum(Zn*self.density_perturbation,axis=2)
    
    @property
    def U(self):
        return self.outputs[self.group_name+self.U_name][()]

    @property
    def ddpsilog_to_delta_factor(self):
        return self.Delta*numpy.sqrt(self.masses*self.THat/self.FSABHat2[:,numpy.newaxis])*self.IHat[:,0,numpy.newaxis]/(self.psiAHat*self.Z)
    
    @property
    def deltaN(self):
        try:
            return self.outputs[self.group_name+self.deltaN_name][()]
        except KeyError:
            return numpy.abs(self.ddpsilog_to_delta_factor * self.dnHatdpsiN/self.nHat)

    @property
    def deltaN_ped(self):
        return self.deltaN[self.mid_pedestal_point_index]
        
    @property
    def deltaEta(self):
        return numpy.abs(self.ddpsilog_to_delta_factor * self.detaHatdpsiN/self.etaHat)
        try:
            return self.outputs[self.group_name+self.deltaEta_name][()]
        except KeyError:
            return numpy.abs(self.ddpsilog_to_delta_factor * self.detaHatdpsiN/self.etaHat)

    @property
    def deltaT(self):
        try:
            return self.outputs[self.group_name+self.deltaT_name][()]
        except KeyError:
            return  numpy.abs(self.ddpsilog_to_delta_factor *self.dlogTHatdpsiN)
            

    @property
    def masses(self):
        m= numpy.array(self.inputs.masses)
        if m.shape == ():
            m.shape = (1,)
        return m
        
    @property
    def THat(self):
        try:
            return self.input_profiles[self.input_profiles_groupname+"THats"][()]
        except AttributeError:
            print "perfect_simulatio: Warning: THat not found in input profiles"
        try:
            return self.outputs[self.group_name+self.THat_name][()]
        except KeyError:
            print "THat could not be obtained since no external profiles have been specified and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def dTHatdpsiN(self):
        try:
            scale_factor=(self.psiAHat/self.psiAHatArray)
            return scale_factor[:,numpy.newaxis]*self.input_profiles[self.input_profiles_groupname+"dTHatdpsis"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.dTHatdpsiN_name][()]
        except KeyError:
            print "PhiHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def numerical_dTHatdpsiN(self):
        return self.ddpsiN(self.THat)

    @property
    def dlogTHatdpsiN(self):
        return self.dTHatdpsiN/self.THat

    @property
    def mdlogTHatdpsiN(self):
        return -self.dlogTHatdpsiN


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
    def dnHatdpsiN(self):
        try:
            return numpy.true_divide(self.psiAHat,self.psiAHatArray)[:,numpy.newaxis]*self.input_profiles[self.input_profiles_groupname+"dnHatdpsis"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.dnHatdpsiN_name][()]
        except KeyError:
            print "nHat could not be obtained since no external profiles have been specified and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

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
            return (self.psiAHat/self.psiAHatArray)[:,numpy.newaxis]*self.input_profiles[self.input_profiles_groupname+"detaHatdpsis"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.detaHatdpsiN_name][()]
        except KeyError:
            print "etaHat could not be obtained since no external profiles have been specified and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def numerical_detaHatdpsiN(self):
        return self.ddpsiN(self.etaHat)
        

    @property
    def pHat(self):
        return self.nHat*self.THat

    @property
    def dpHatdpsiN(self):
        return self.dnHatdpsiN*self.THat + self.nHat*self.dTHattdpsiN

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
            return (self.psiAHat/self.psiAHatArray)*self.input_profiles[self.input_profiles_groupname+"dPhiHatdpsi"][()]
        except AttributeError:
            pass
        try:
            return self.outputs[self.group_name+self.dPhiHatdpsiN_name][()]
        except KeyError:
            print "PhiHat could not be obtained since no external profiles have been speciied and simulation output probably does not exist. Try running perfect with solveSystem=.false. to generate the inputs."

    @property
    def dPhiHat_dpsiN(self):
        #intentionally "crude" to use same deriv approx as dGamma_dpsiN
        #psiN=self.actual_psiN
        #PhiHat=self.PhiHat
        #ret=numpy.zeros([len(psiN)-1,self.num_species])
        #for i in range(1,len(psiN)):
        #    dpsiN=psiN[i]-psiN[i-1]
        #    dPhiHat=PhiHat[i]-PhiHat[i-1]
        #    #print dPhiHat/dpsiN
        #    ret[i-1]=dPhiHat/dpsiN
        #return ret
        return self.ddpsiN(self.PhiHat,order=4)

    @property
    def dtheta(self):
        #corresponding dpsiN does not exist atm, since non-uniform grid.
        return self.theta[1]-self.theta[0]

    @property
    def Npsi(self):
        return len(self.psi)

    @property
    def Ntheta(self):
        return len(self.theta)

    @property
    def Nxi(self):
        return self.inputs.Nxi

    def Nxi_for_x(self,ix,f=True):
        if f==True:
            # ix input is fortran indexed
            # subtract by one to translate to python indexing
            ix = ix-1
        return self.outputs[self.group_name+self.Nxi_for_x_name][()][ix]

    def first_index_for_x(self,ix,f=True):
        if f:
            # ix input is fortran indexed
            # subtract by one to translate to python indexing
            ix = ix-1
        if ix == 0:
            return 0
        else:
            return self.first_index_for_x(ix-1,f=False) + self.Nxi_for_x(ix-1,f=False)

    def min_x_for_L(self,L):
        if L==0:
            return 1
        a=numpy.ones(self.Nxi,dtype=numpy.uint64)
        for j in range(1,self.Nx):
            # print self.Nxi_for_x(j)
            a[self.Nxi_for_x(j):] = j+1
        return a[L]
  
    @property
    def Nx(self):
        return self.inputs.Nx

    @property
    def NL(self):
        # L for rosenbluth potential
        return self.inputs.NL

    @property
    def local_matrix_size(self):
        # 2017-02-06: matches interneally computed in local simulation
        return self.local_DKE_matrix_size*self.Nspecies

    @property
    def local_DKE_matrix_size(self):
        # 2017-02-06: matches interneally computed in local simulation
        return self.Ntheta * sum([self.Nxi_for_x(i) for i in range(self.Nx)]) 

    @property
    def Nspecies(self):
        #stupidity wrapper since I did not realise I had a num_species and made a new one.
        return self.num_species 
    
    @property
    def dTHat_dpsiN(self):
        #intentionally "crude" to use same deriv approx as dGamma_dpsiN
        #psiN=self.actual_psiN
        #THat=self.THat
        #ret=numpy.zeros([len(psiN)-1,self.num_species])
        #for i in range(1,len(psiN)):
        #    dpsiN=psiN[i]-psiN[i-1]
        #    dTHat=THat[i]-THat[i-1]
        #    #print dTHat/dpsiN
        #    ret[i-1]=dTHat/dpsiN
        return self.ddpsiN(self.THat,order=4)

    @property
    def Zeff(self):
        electron_index=self.species_list.index("e")
        Z=self.Z
        Z[electron_index]=0 # to exclude electrons from calculation
        return numpy.sum(Z**2*self.nHat,axis=1)/numpy.sum(Z*self.nHat,axis=1)

    @property
    def alpha(self,main_index=0):
        print "perfect_simulation: warning: alpha: Assumes main species has index " + str(main_index)
        self.Zeff - self.Z[main_index]
        
    @property
    def source_charge(self):
        return numpy.sum(self.Z*self.GammaHat_source,axis=1)

    @property
    def charge_source(self):
        return self.source_charge

    @property
    def GammaHat_source_for_export(self):
        #for comparrison with dGammaHat/dPsiN
        VPrimeHat=numpy.fabs(self.VPrimeHat)
        prefact=1#(self.psiAHatArray[:,numpy.newaxis]/self.psiAHat)
        #print self.psiAHat
        return self.THat**(3./2.)*self.particle_source/(self.masses**2)
    
    @property
    def charge_source_for_export(self):
        return numpy.sum(self.Z*self.GammaHat_source_for_export,axis=1)
    
    @property
    def total_source_charge(self):
        final_i=-30
        source_charge=numpy.sum(self.Z*self.GammaHat_source,axis=1)
        #probably OK with non-uniform grid...
        return scipy.integrate.simps(source_charge[0:final_i],self.actual_psiN[0:final_i])
    
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
        try:
            return self.outputs[self.group_name+self.psi_name][()]
        except (KeyError,TypeError):
            #print "cannot get internal psi from output, will generate uniform from input."
            return self.inputs.psi

    @property
    def pedestal_start_stop_indices(self):
        return [get_index_range(self.actual_psiN,[point,point])[1] for point in  self.pedestal_start_stop]

    @property
    def mid_pedestal_point_index(self):
        point = (self.pedestal_start_stop[0] + self.pedestal_start_stop[1])/2.0
        return get_index_range(self.actual_psiN,[point,point])[1]

    @property
    def psi_index(self):
        return numpy.array(range(self.Npsi))
        
    @property
    def sqrtpsi(self):
        return numpy.sqrt(self.actual_psi)

    @property
    def psiAHat(self):
        try:
            return self.outputs[self.group_name+self.psiAHat_name][()]
        except (KeyError,TypeError):
            return self.inputs.psiAHat
            
    @property
    def theta(self):
        return self.outputs[self.group_name+self.theta_name][()]

    @property
    def theta_shifted(self,shift=numpy.pi):
        return self.shift_theta(shift)

    @property
    def num_species(self):
        try:
            return self.outputs[self.group_name+self.num_species_name][()]
        except (KeyError,TypeError):
            return len(self.masses)

    @property
    def collisionality(self):
        return self.outputs[self.group_name+self.collisionality_name][()]

    @property
    def density_perturbation(self):
        #non-adiabatic part
        return self.outputs[self.group_name+self.density_perturbation_name][()]

    @property
    def FSA_density_perturbation(self):
        #non-adiabatic part
        return self.outputs[self.group_name+self.FSA_density_perturbation_name][()]

    
    @property
    def density_perturbation_shifted(self,shift=numpy.pi):
        #non-adiabatic part
        (theta,density_perturbation_shifted) = self.shift_theta(shift,self.density_perturbation)
        return density_perturbation_shifted
    
    @property
    def pressure_perturbation(self):
        #non-adiabatic part
        return self.outputs[self.group_name+self.pressure_perturbation_name][()]

    @property
    def FSA_pressure_perturbation(self):
        #non-adiabatic part
        return self.outputs[self.group_name+self.FSA_pressure_perturbation_name][()]

    
    @property
    def total_density_perturbation(self):
        # n_1/n_0
        ZoT=numpy.expand_dims(self.Z/self.THat, axis=1)
        ret=self.density_perturbation - (2*self.omega/self.Delta)*ZoT*self.potential_perturbation[:,:,numpy.newaxis]
        return ret

    @property
    def dGammaHat_dpsiN(self):
        GammaHat=self.particle_flux
        return self.ddpsiN(GammaHat,order=4)

    @property
    def dPiHat_dpsiN(self):
        PiHat=self.momentum_flux
        return self.ddpsiN(PiHat,order=4)

    @property
    def jprefac_dPi_dpsiN(self):
        #psiN=self.actual_psiN
        #Pi=self.sum_m_momentum_flux
        #ret=numpy.zeros(len(psiN))
        #for i in range(1,len(psiN)-1):
        #    dpsiN=psiN[i]-psiN[i-1]
        #    dPi=Pi[i+1]-Pi[i-1]     
        #    #print dGammaHat/dpsiN
        #    ret[i]=dPi/(2.0*dpsiN)
        Pi=self.sum_m_momentum_flux
        return self.ddpsiN(Pi,order=4)*self.Delta/self.psiAHat
    

    @property
    def djHat_dpsiN(self):
        return numpy.sum(self.Z*self.dGammaHat_dpsiN,axis=1)



    @property
    def dqHat_dpsiN(self):
        #non-adiabatic part
        #psiN=self.actual_psiN
        #qHat=self.conductive_heat_flux
        #ret=numpy.zeros([len(psiN),self.num_species])
        #for i in range(1,len(psiN)):
        #    dpsiN=psiN[i]-psiN[i-1]
        #    dqHat=qHat[i]-qHat[i-1]
        #    #print dGammaHat/dpsiN
        #    ret[i-1]=dqHat/dpsiN
        #ret[0]=0
        #return ret
        qHat=self.conductive_heat_flux
        return self.ddpsiN(qHat,order=4)

    @property
    def dQHat_dpsiN(self):
        QHat=self.heat_flux
        return self.ddpsiN(QHat,order=4)

    @property
    def ddpsiN_sum_QHat(self):
        return self.ddpsiN(self.sum_QHat,order=4)

    @property
    def IHat(self):
        return nstack(self.outputs[self.group_name+self.IHat_name][()],axis=1,n=len(self.theta)) #add theta axis
    
    @property
    def FSA_IHat(self):
        return self.outputs[self.group_name+self.IHat_name][()]
    
    @property
    def ddpsiN_IHat(self):
        return self.outputs[self.group_name+self.dIHatdpsiN_name][()]
    
    @property
    def RHat(self):
        RHat = self.outputs[self.group_name+self.RHat_name][()]
        if arraylist_rank(RHat) == 2:
            return RHat
        elif arraylist_rank(RHat) == 1:
            #for compatibility with old 1D RHat output
            return nstack(RHat,axis=0,n=len(self.psi)) #add psi axis

    @property
    def RHat_mid_ped(self):
        self.RHat[self.mid_pedestal_point_index]
        
    @property
    def FSARHat(self):
        return self.FSA(self.RHat)

    @property
    def JHat(self):
        return self.outputs[self.group_name+self.JHat_name][()]

    @property
    def FSA_JHat(self):
        return self.FSA(self.JHat)

    @property
    def JHat_mid_ped(self):
        return self.JHat[self.mid_pedestal_point_index]
    
    @property
    def BHat(self):
        return self.outputs[self.group_name+self.BHat_name][()]

    @property
    def ddpsiN_BHat(self):
        return self.outputs[self.group_name+self.dBHatdpsiN_name][()]

    @property
    def ddtheta_BHat(self):
        return self.outputs[self.group_name+self.dBHatdtheta_name][()]
    
    
    @property
    def BHat_mid_ped(self):
        return self.BHat[self.mid_pedestal_point_index]

    @property
    def ddpsiN_BHat_mid_ped(self):
        return self.ddpsiN_BHat[self.mid_pedestal_point_index]
    
    @property
    def ddtheta_BHat_mid_ped(self):
        return self.ddtheta_BHat[self.mid_pedestal_point_index]
    
    
    @property
    def FSABHat2(self):
        return self.outputs[self.group_name+self.FSABHat2_name][()]

    @property
    def Bt(self):
        #will always be positive for PERFECT Miller (in toroidal direction)
        return self.IHat/self.RHat

    @property
    def q(self):
        return self.outputs[self.group_name+self.q_name][()]

    @property
    def epsilon(self):
        return self.outputs[self.group_name+self.epsilon_name][()]

    @property
    def Bp(self):
        #will typically be negative for PERFECT Miller (opposite to poloidal direction)
        return numpy.sign(self.JHat)*numpy.sqrt(self.BHat**2-self.Bt**2)

    @property
    def Bp_mid_ped(self):
        self.Bp[self.mid_pedestal_point_index]

    
    @property
    def FSABp(self):
        return self.FSA(self.Bp)

    @property
    def sqrtFSABp2(self):
        #another typical Bp
        return numpy.sqrt(self.FSA(self.Bp**2))

    @property
    def FSABHat2_test(self):
        #to compare external FSA to PERFECT internal
        #2016-06-03: differ at 0.01
        return self.FSA(self.BHat**2)
    
    @property
    def FSABp_old(self):
        #a typical Bp
        return numpy.trapz(numpy.fabs(self.Bp)/numpy.fabs(self.JHat),dx=self.dtheta)/numpy.fabs(self.VPrimeHat)

    @property
    def sqrtFSABp2_old(self):
        #another typical Bp
        return numpy.sqrt(numpy.trapz(numpy.fabs(self.Bp**2)/numpy.fabs(self.JHat),dx=self.dtheta)/numpy.fabs(self.VPrimeHat))

    @property
    def FSABHat2_test_old(self):
        #to compare external FSA to PERFECT internal
        #2016-06-03: differ at 0.01
        return numpy.trapz(self.BHat**2/numpy.fabs(self.JHat),dx=self.dtheta)/numpy.fabs(self.VPrimeHat)

    @property
    def Bp_at_psi_of_theta(self):
        psiN_point=self.pedestal_point
        #for PERFECT Mille geometry, Bp does not depend on psi
        return self.attrib_at_psi_of_theta("Bp",psiN_point)

    @property
    def poloidal_flow_parallel(self):
        return numpy.expand_dims(self.Bp/self.BHat,axis=2)*self.flow

    @property
    def poloidal_flow_ExB(self):
        return self.omega/(self.Delta*self.psiAHat)*numpy.expand_dims(self.Bp*self.IHat/self.BHat**2,axis=2)*self.density_perturbation*numpy.expand_dims(numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=2)

    @property
    def poloidal_flow_gradP(self):
        return self.Delta/(2*self.psiAHat)*(self.masses/self.Z)*(numpy.expand_dims(self.Bp*self.IHat/(self.BHat**2),axis=2)/numpy.expand_dims(self.nHat,axis=1))*self.magnetization_flow_perturbation

    @property
    def magnetization_flow_perturbation_test(self):
        return self.magnetization_flow_perturbation_test_order(order=4)

    @property
    def magnetization_flow_perturbation_test_new(self):
        #to test new external derivative routine
        return self.ddpsiN(self.magnetization_perturbation,order=4)

    
    
    def magnetization_flow_perturbation_test_order(self,order=4):
        #for testing different derivative approximations
        if order == 1:
            #forward except at the first point, where it is backward
            mp=self.magnetization_perturbation
            psiN=self.actual_psiN
            for i in range(0,len(psiN)):
                dpsiN=psiN[i]-psiN[i-1]
                if i==0:
                    mpf=[(mp[i+1]-mp[i])/dpsiN]
                else:
                    mpf=mpf+[(mp[i]-mp[i-1])/dpsiN]
            mpf=numpy.array(mpf)
            #print mpf
            return mpf
        elif order == 2:
            #centered 3-point, forward/backward 3-point at the ends
            mp=self.magnetization_perturbation
            ddpsiN=diff_matrix(self.psi[0],self.psi[-1],self.Npsi,order=2)
            return (self.psiAHat/self.psiAHatArray)[:,numpy.newaxis,numpy.newaxis]*numpy.tensordot(ddpsiN,mp,([1],[0])) 
        elif order == 4:
            #centered 5-point, forward/backward 5-point at the ends
            mp=self.magnetization_perturbation
            ddpsiN=diff_matrix(self.psi[0],self.psi[-1],self.Npsi,order=4)
            return (self.psiAHat/self.psiAHatArray)[:,numpy.newaxis,numpy.newaxis]*numpy.tensordot(ddpsiN,mp,([1],[0]))
        
    
    @property
    def poloidal_flow_gradP_test(self):
        return self.Delta/(2*self.psiAHat)*(self.masses/self.Z)*(numpy.expand_dims(self.Bp*self.IHat/(self.BHat**2),axis=2)/numpy.expand_dims(self.nHat,axis=1))*self.magnetization_flow_perturbation_test

    
    @property
    def poloidal_flow_inputs(self):
        return 1/(2*self.psiAHat)*(numpy.expand_dims(self.THat,axis=1)/self.Z)*numpy.expand_dims(self.Bp*self.IHat/(self.BHat**2),axis=2)*numpy.expand_dims(self.dnHatdpsiN/self.nHat + self.dTHatdpsiN/self.THat +(2*self.omega*self.Z/(self.Delta*self.THat))*numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=1)

    @property
    def poloidal_flow_inputs_parallel(self):
        return self.poloidal_flow_parallel + self.poloidal_flow_inputs
    
    @property
    def poloidal_flow_test(self):
        if self.local:
            #local
            return self.poloidal_flow_parallel + self.poloidal_flow_inputs
        elif self.no_ddpsi:
            # no v_m \cdot \nabla g
            return self.poloidal_flow_parallel + self.poloidal_flow_ExB + self.poloidal_flow_inputs
        else:
            #global
            return self.poloidal_flow_parallel + self.poloidal_flow_ExB + self.poloidal_flow_gradP_test + self.poloidal_flow_inputs

    @property
    def k_poloidal_parallel(self):
        return self.poloidal_flow_parallel*self.poloidal_flow_to_k_poloidal_factor
        
    @property
    def k_poloidal_ExB(self):
        return self.poloidal_flow_ExB*self.poloidal_flow_to_k_poloidal_factor
    
    @property
    def k_poloidal_gradP(self):
        return self.poloidal_flow_gradP*self.poloidal_flow_to_k_poloidal_factor
    

    @property
    def k_poloidal_gradP_test(self):
        return self.poloidal_flow_gradP_test*self.poloidal_flow_to_k_poloidal_factor
    
    @property
    def k_poloidal_inputs(self):
        return self.poloidal_flow_inputs*self.poloidal_flow_to_k_poloidal_factor
    
    @property
    def k_poloidal_inputs_parallel(self):
        return self.k_poloidal_inputs + self.k_poloidal_parallel
    
    @property
    def k_poloidal_test(self):
        return self.k_poloidal_parallel + self.k_poloidal_ExB + self.k_poloidal_gradP_test + self.k_poloidal_inputs

    
    @property
    def toroidal_flow_parallel(self):
        return numpy.expand_dims(self.Bt/self.BHat,axis=2)*self.flow

    @property
    def toroidal_flow_ExB(self):
        return self.omega/(self.Delta*self.psiAHat)*numpy.expand_dims(self.Bp**2*self.RHat/self.BHat**2,axis=2)*self.density_perturbation*numpy.expand_dims(numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=2)

    @property
    def toroidal_flow_gradP(self):
        return (self.Delta/(2*self.psiAHat))*(self.masses/self.Z)*(numpy.expand_dims(self.Bp**2*self.RHat/(self.BHat**2),axis=2)/numpy.expand_dims(self.nHat,axis=1))*self.magnetization_flow_perturbation

    @property
    def toroidal_flow_gradP_test(self):
        return (self.Delta/(2*self.psiAHat))*(self.masses/self.Z)*(numpy.expand_dims(self.Bp**2*self.RHat/(self.BHat**2),axis=2)/numpy.expand_dims(self.nHat,axis=1))*self.magnetization_flow_perturbation_test

    @property
    def toroidal_flow_inputs(self):
        return 1/(2*self.psiAHat)*(numpy.expand_dims(self.THat,axis=1)/self.Z)*numpy.expand_dims(self.Bp**2*self.RHat/(self.BHat**2),axis=2)*numpy.expand_dims(self.dnHatdpsiN/self.nHat + self.dTHatdpsiN/self.THat +(2*self.omega*self.Z/(self.Delta*self.THat))*numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=1)
    
    @property
    def toroidal_flow_test(self):
        if self.local:
            #local
            return self.toroidal_flow_parallel - self.toroidal_flow_inputs
        elif self.no_ddpsi:
            # no v_m \cdot \nabla g
            return self.toroidal_flow_parallel - self.toroidal_flow_ExB  - self.toroidal_flow_inputs
        else:
            #global
            return self.toroidal_flow_parallel - self.toroidal_flow_ExB - self.toroidal_flow_gradP_test - self.toroidal_flow_inputs
        
    @property
    def perpendicular_flow_inputs(self):
        return self.perpendicular_flow_inputs_diamagnetic + self.perpendicular_flow_inputs_ExB

    @property
    def perpendicular_flow_inputs_ExB(self):
        return (self.omega/(self.Delta*self.psiAHat*self.Z))*numpy.expand_dims(self.Bp*self.RHat/self.BHat,axis=2)*numpy.expand_dims(self.Z*numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=1)
    #self.dPhiHatdpsiN[:,numpy.newaxis,numpy.newaxis]

    @property
    def perpendicular_flow_inputs_diamagnetic(self):
        return 1/(2*self.psiAHat)*(numpy.expand_dims(self.THat,axis=1)/self.Z)*numpy.expand_dims(self.Bp*self.RHat/self.BHat,axis=2)*numpy.expand_dims(self.dnHatdpsiN/self.nHat + self.dTHatdpsiN/self.THat,axis=1)

    @property
    def perpendicular_flow_inputs_test(self):
        #to see if
        #perpendicular_flow_inputs = perpendicular_flow_inputs_ExB + perpendicular_flow_inputs_diamagnetic
        #which it is 2016-05-27
        return self.perpendicular_flow_inputs_diamagnetic + self.perpendicular_flow_inputs_ExB

    @property
    def perpendicular_flow_ExB(self):
        #equal to flow ExB sqrt(Vp^2+Vt^2).
        return self.omega/(self.Delta*self.psiAHat)*numpy.expand_dims(self.Bp*self.RHat/self.BHat,axis=2)*self.density_perturbation*numpy.expand_dims(numpy.expand_dims(self.dPhiHatdpsiN,axis=1),axis=2)

    @property
    def perpendicular_flow_gradP(self):
        #equal to flow gradP sqrt(Vp^2+Vt^2).
        return (self.Delta/(2*self.psiAHat))*(self.masses/self.Z)*(numpy.expand_dims(self.Bp*self.RHat/(self.BHat),axis=2)/numpy.expand_dims(self.nHat,axis=1))*self.magnetization_flow_perturbation

    @property
    def perpendicular_flow(self):
        return self.perpendicular_flow_inputs + self.perpendicular_flow_ExB + self.perpendicular_flow_gradP

    @property
    def perpendicular_flow_over_vT(self):
        self.Delta*nstack(numpy.sqrt(self.masses/self.THat),axis=1,n=len(self.theta))*self.perpendicular_flow


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
    def __init__(self,input_filename,norm_filename,species_filename=None,output_filename=None,psiN_to_psi_filename=None,global_term_multiplier_filename=None,group_name=None,pedestal_start_stop=(None,None),pedestal_point = None,core_point=None):
        perfect_simulation.__init__(self,input_filename,output_filename,species_filename,psiN_to_psi_filename,global_term_multiplier_filename,group_name,pedestal_start_stop,pedestal_point,core_point)

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
    def a(self):
        return self.epsilon*self.RBar

    @property
    def a_theta(self):
        return self.a*self.theta

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
    def p(self):
        return self.n*self.T

    @property
    def vT(self):
        return self.vBar*numpy.sqrt(self.THat/self.masses)
         
        
    @property
    def rho_p(self):
        return self.mBar*self.masses*self.vT/(self.eBar*self.Z*numpy.expand_dims(self.FSABp,axis=1)*self.BBar)


    @property
    def normed_particle_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.particle_flux/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_particle_flux_psi_unit_vector(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.particle_flux_psi_unit_vector/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.nBar*self.vBar

    @property
    def normed_particle_flux_psi_unit_vector_no_fsa(self):
        return numpy.pi*self.Delta**2*self.nBar*self.vBar*self.particle_flux_psi_unit_vector_no_fsa

    @property
    def pw_normed_particle_flux_psi_unit_vector_no_fsa(self,i_s=0):
        # pedestal width normed flux
        return self.normed_particle_flux_psi_unit_vector_no_fsa/self.pedestal_width

    @property
    def pw_nbar_normed_particle_flux_psi_unit_vector_no_fsa(self,i_s=0):
        # pedestal width normed flux over nBar
        return numpy.pi*self.Delta**2*self.vBar*self.particle_flux_psi_unit_vector_no_fsa/self.pedestal_width

    
    @property
    def normed_particle_flow_psi_unit_vector_no_fsa(self):
        return numpy.pi*self.Delta**2*self.vBar*self.particle_flux_psi_unit_vector_no_fsa/numpy.expand_dims(self.nHat,axis=1)

    @property
    def normed_particle_flow_psi_unit_vector_no_fsa_km_s(self):
        return 1e-3*self.normed_particle_flow_psi_unit_vector_no_fsa
    
    @property
    def ow_normed_particle_flow_psi_unit_vector_no_fsa(self,i_s=0):
        # i_s: index of species used for orbit width calculation 
        return self.normed_particle_flow_psi_unit_vector_no_fsa/self.FSA_orbit_width[:,i_s,numpy.newaxis,numpy.newaxis]

    @property
    def pw_normed_particle_flow_psi_unit_vector_no_fsa(self,i_s=0):
        # pedestal width normed
        return self.normed_particle_flow_psi_unit_vector_no_fsa/self.pedestal_width

    @property
    def normed_momentum_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.momentum_flux/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*(self.RBar**2)*self.BBar*self.nBar*self.vBar**2

    @property
    def normed_momentum_flux_psi_unit_vector(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.momentum_flux_psi_unit_vector/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.RBar*self.nBar*self.vBar**2

    @property
    def normed_m_momentum_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.mBar*self.m_momentum_flux/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*(self.RBar**2)*self.BBar*self.nBar*self.vBar**2

    @property
    def normed_m_momentum_flux_psi_unit_vector(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.mBar*self.m_momentum_flux_psi_unit_vector/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*self.RBar*self.nBar*self.vBar**2

    @property
    def normed_heat_flux(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.heat_flux/VPrimeHat)*signOfVPrimeHat*self.TBar*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar

    @property
    def normed_heat_flux_psi_unit_vector(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.heat_flux_psi_unit_vector/VPrimeHat)*signOfVPrimeHat*self.TBar*(numpy.pi*self.Delta**2)*self.nBar*self.vBar

    @property
    def normed_conductive_heat_flux(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.conductive_heat_flux/VPrimeHat)*signOfVPrimeHat*self.TBar*(numpy.pi*self.Delta**2)*self.RBar*self.BBar*self.nBar*self.vBar
    #self.normed_heat_flux-(5.0/2.0)*self.T*self.normed_particle_flux

    @property
    def normed_conductive_heat_flux_psi_unit_vector(self):
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.conductive_heat_flux_psi_unit_vector/VPrimeHat)*signOfVPrimeHat*self.TBar*(numpy.pi*self.Delta**2)*self.nBar*self.vBar

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
    def normed_m_momentum_flux_sum(self):
        #to make appropriate size to divide the particle flux with
        VPrimeHat=[self.VPrimeHat]
        VPrimeHat=numpy.dot(numpy.transpose(VPrimeHat),[numpy.array([1]*self.num_species)])
        signOfVPrimeHat=numpy.sign(VPrimeHat)
        return (self.mBar*self.m_momentum_flux_sum/VPrimeHat)*signOfVPrimeHat*(numpy.pi*self.Delta**2)*(self.RBar**2)*self.BBar*self.nBar*self.vBar**2
    
    @property
    def normed_heat_source(self):
        return self.Delta*self.nBar/(self.vBar**2*self.RBar*numpy.sqrt(self.masses))*self.heat_source

    @property
    def normed_particle_source(self):
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

    @property
    def normed_poloidal_flow(self):
        return self.Delta*self.vBar*self.poloidal_flow

    @property
    def normed_poloidal_flow_km_s(self):
        return 1e-3*self.normed_poloidal_flow

    @property
    def a_normed_poloidal_flow(self):
        return self.normed_poloidal_flow/(self.a*2*numpy.pi)

    @property
    def a_normed_poloidal_flux(self):
        return self.a_normed_poloidal_flow*numpy.expand_dims(self.n,axis=1)

    @property
    def vorticity(self):
        # psiN derivative
        #return self.ddpsiN(self.normed_poloidal_flux) - self.ddtheta(self.normed_particle_flux_psi_unit_vector_no_fsa)
        return self.ddr(self.normed_poloidal_flux) - self.ddtheta(self.normed_particle_flux_psi_unit_vector_no_fsa)

    @property
    def vorticity_over_nPed(self):
        psiN_point=self.core_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        npeds=[self.n[psiN_index,i] for i in range(self.num_species)]
        return self.vorticity/npeds
        
    
    @property
    def a_nbar_normed_poloidal_flux(self):
        return self.a_normed_poloidal_flow*numpy.expand_dims(self.nHat,axis=1)
    
    @property
    def normed_poloidal_flux(self):
        return self.Delta*self.vBar*numpy.expand_dims(self.n,axis=1)*self.poloidal_flow
    
    @property
    def normed_perpendicular_flow(self):
        return self.Delta*self.vBar*self.perpendicular_flow

    @property
    def normed_perpendicular_flow_inboard(self):
        return self.Delta*self.vBar*self.perpendicular_flow_inboard

    @property
    def normed_perpendicular_flow_outboard(self):
        return self.Delta*self.vBar*self.perpendicular_flow_outboard

    @property
    def normed_toroidal_flow(self):
        return self.Delta*self.vBar*self.toroidal_flow

    @property
    def normed_toroidal_flow_km_s(self):
        return 1e-3*self.normed_toroidal_flow

    @property
    def normed_FSA_toroidal_flow(self):
        return self.Delta*self.vBar*self.FSA_toroidal_flow

    @property
    def normed_n_FSA_toroidal_flow(self):
        return self.n*self.normed_FSA_toroidal_flow

    @property
    def normed_n_FSA_toroidal_mass_flow(self):
        return self.mBar*self.masses*self.normed_n_FSA_toroidal_flow

    @property
    def sum_normed_n_FSA_toroidal_mass_flow(self):
        return numpy.sum(self.normed_n_FSA_toroidal_mass_flow,axis=1)

    @property
    def RBar_nabla_psiN(self):
        # old hardcoded value: RBar_dpsiN_over_dr =4.92
        return numpy.fabs(self.RHat*self.Bp/self.psiAHat)

    @property
    def RBar_nabla_psiN_of_theta(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta(self.RBar_nabla_psiN,psiN_point)

    @property
    def FSA_RBar_nabla_psiN(self):
        return self.FSA(self.RBar_nabla_psiN)
    

    @property
    def orbit_width(self):
        return numpy.fabs((self.Delta/self.psiAHat)*numpy.expand_dims(self.RHat,axis=2)*numpy.expand_dims(numpy.sqrt(self.epsilon*self.THat*self.masses)/self.Z,axis=1))
    
    @property
    def orbit_width_of_theta(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta(self.orbit_width,psiN_point)

    @property
    def FSA_orbit_width(self):
        return self.FSA(self.orbit_width)

    @property
    def ped_FSA_orbit_width(self):
        psiN_point=(self.pedestal_start_stop[-1] - self.pedestal_start_stop[0])/2.0
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        return self.FSA_orbit_width[psiN_point]


    @property
    def FSA_orbit_width_over_ped_width(self):
        return self.FSA_orbit_width/(self.pedestal_start_stop[-1] - self.pedestal_start_stop[0])

    @property
    def ped_FSA_ion_orbit_width_over_ped_width(self):
        ii = 0 # ion species index
        psiN_point=(self.pedestal_start_stop[-1] - self.pedestal_start_stop[0])/2.0
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point])[1]
        return self.FSA_orbit_width_over_ped_width[ii,psiN_index]

    @property
    def orbit_width_over_rN(self):
        return numpy.fabs(self.FSA_orbit_width*self.dnHatdpsiN/self.nHat)

    @property
    def orbit_width_old(self):
        RBar_dpsiN_over_dr =4.92
        return numpy.fabs(numpy.expand_dims(numpy.sqrt(self.epsilon*self.THat*self.masses),axis=1)/(self.Z*numpy.expand_dims(self.Bp,axis=2))*self.Delta*RBar_dpsiN_over_dr)

    @property
    def orbit_width_of_theta_old(self):
        psiN_point=self.pedestal_point
        return self.attrib_at_psi_of_theta(self.orbit_width_old,psiN_point)

    @property
    def FSA_orbit_width_old(self):
        return self.FSA(self.orbit_width_old)


    @property
    def psiN2(self):
        return self.actual_psiN/self.pedestal_start_stop[-1]

    @property
    def psiN3(self):
        #shift is positive if LCFS needs to be shifted upwards
        # i.e. pedestal stop is below 1
        shift = 1 - self.pedestal_start_stop[-1]
        return (self.actual_psiN + shift)
    
    @property
    def psi_o(self,i_s=0):
        #i_s : index of species for which to calculate orbit width for
        ped_stop = self.pedestal_start_stop[-1]
        psi_o = (self.actual_psiN-ped_stop)/self.FSA_orbit_width[:,i_s]
        return psi_o

    def psi_o_to_psiN(self,psi_o):
        #psi_o : psi_o value to get psi_N for
        i=get_index_range(self.psi_o,[psi_o,psi_o],ret_range=False)[0]
        psiN=self.actual_psiN[i]
        return psiN

    def psiN_to_psi_o(self,psiN):
        #psiN : psiN value to get psi_o for
        i=get_index_range(self.actual_psiN,[psiN,psiN],ret_range=False)[0]
        psi_o=self.psi_o[i]
        return psi_o

    @property
    def pedestal_point_psi_o(self):
        return self.psiN_to_psi_o(self.pedestal_point)
        

    @property
    def r(self):
        #dr_dpsiN : dr/dpsiN
        return (self.actual_psiN - self.pedestal_start_stop[-1])*self.FSA_RBar_nabla_psiN/self.RBar

    def ddr(self,attrib,order=4,scale=True):
        return self.ddpsiN(attrib,order=4,scale=True)*numpy.expand_dims(self.RBar_nabla_psiN/self.RBar,axis=2)
        

    @property
    def pedestal_width(self):
        #in meters
        return numpy.fabs((self.pedestal_start_stop[0] - self.pedestal_start_stop[-1])*self.FSA_RBar_nabla_psiN[0]/self.RBar)

    @property
    def pedestal_width_psiN(self):
        return numpy.fabs(self.pedestal_start_stop[-1] - self.pedestal_start_stop[0])

    @property
    def pedestal_width_psiN2(self):
        return numpy.fabs(self.pedestal_start_stop_psiN2[-1] - self.pedestal_start_stop_psiN2[0])

    @property
    def pedestal_width_psiN3(self):
        return numpy.fabs(self.pedestal_start_stop_psiN3[-1] - self.pedestal_start_stop_psiN3[0])
    
    @property
    def r_old(self, dr_dpsiN=0.591412216405):
        #dr_dpsiN : dr/dpsiN
        return (self.actual_psiN - self.pedestal_start_stop[-1])*dr_dpsiN

    @property
    def pedestal_start_stop_psiN(self):
        return self.pedestal_start_stop

    
    @property
    def pedestal_start_stop_psiN2(self):
        return [point/self.pedestal_start_stop[-1] for point in self.pedestal_start_stop]

    @property
    def pedestal_start_stop_psiN3(self):
        shift = 1 - self.pedestal_start_stop[-1]
        return [point + shift for point in self.pedestal_start_stop]
    
    @property
    def pedestal_start_stop_psi_o(self):
        return [self.psi_o[point] for point in self.pedestal_start_stop_indices]

    @property
    def pedestal_start_stop_r(self):
        return [self.r[point] for point in self.pedestal_start_stop_indices]
    

    def Prandtl_proxy(self,ion_index = 0):
        psiN_point = self.pedestal_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point],ret_range=False)[0]
        return self.normed_m_momentum_flux*((1/self.normed_conductive_heat_flux[:,ion_index])*(self.p[:,ion_index])/(self.RBar*self.normed_n_FSA_toroidal_mass_flow[:,ion_index])*(self.deltaT[:,ion_index]/self.deltaN[:,ion_index]))[psiN_index]

    def Prandtl_proxy_unit_vector(self,ion_index = 0):
        psiN_point = self.pedestal_point
        print psiN_point
        psiN_index=get_index_range(self.actual_psiN,[psiN_point,psiN_point],ret_range=False)[0]
        return self.normed_m_momentum_flux_psi_unit_vector*((1/self.normed_conductive_heat_flux_psi_unit_vector[:,ion_index])*(self.p[:,ion_index])/(self.RBar*self.normed_n_FSA_toroidal_mass_flow[:,ion_index])*(self.deltaT[:,ion_index]/self.deltaN[:,ion_index]))[psiN_index]

    @property
    def Prandtl_proxy1(self):
        #over Delta(massflow_toroidal)/Delta(psi)
        return self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*2.54*numpy.expand_dims(numpy.fabs(self.FSABp*self.BBar*self.RBar),axis=1))/(2**(3.0/2.0)*self.T**(3.0/2.0)*(self.masses*self.mBar)**(1.0/2.0)/(self.eBar**2*self.BBar**2*0.7))

    @property
    def Prandtl_proxy2(self):
        #assume \nabla n_a V_t = \delta_n n_a V_ta/\rho_pa
        nabla_toroidal_mass_flow=numpy.sum(numpy.expand_dims(self.sum_normed_n_FSA_toroidal_mass_flow,axis=1)*self.deltaN/self.rho_p,axis=1)
        return self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*numpy.expand_dims(nabla_toroidal_mass_flow,axis=1))/(2**(3.0/2.0)*self.T**(3.0/2.0)*(self.masses*self.mBar)**(1.0/2.0)/(self.eBar**2*self.BBar**2*0.7))

    @property
    def Prandtl_proxy2_1(self):
        #assume \nabla n_a V_t = \delta_n n_a V_ta/\rho_pa
        nabla_toroidal_mass_flow=numpy.max(numpy.sum(numpy.expand_dims(self.sum_normed_n_FSA_toroidal_mass_flow,axis=1)*self.deltaN/self.rho_p,axis=1))
        return self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*numpy.expand_dims(nabla_toroidal_mass_flow,axis=1))/(2**(3.0/2.0)*self.T**(3.0/2.0)*(self.masses*self.mBar)**(1.0/2.0)/(self.eBar**2*self.BBar**2*0.7))

    
    @property
    def Prandtl_proxy3(self):
        print numpy.fabs(self.FSABp*self.BBar*self.RBar)
        #over Delta(massflow_toroidal)/Delta(psi)
        q_max=numpy.max(numpy.sum(self.normed_conductive_heat_flux_psi_unit_vector,axis=1))
        #print "q_max in keV: " + str(q_max/(1000*self.eBar))
        return self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*2.54*numpy.expand_dims(numpy.fabs(self.FSABp*self.BBar*self.RBar),axis=1))/(q_max/(self.p*self.deltaT/self.rho_p))

    @property
    def Prandtl_proxy3_1(self):
        print numpy.fabs(self.FSABp*self.BBar*self.RBar)
        #over Delta(massflow_toroidal)/Delta(psi)
        q_max=numpy.max(numpy.sum(self.normed_conductive_heat_flux_psi_unit_vector,axis=1))
        #print "q_max in keV: " + str(q_max/(1000*self.eBar))
        return self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*2.54*numpy.expand_dims(numpy.fabs(self.FSABp*self.BBar*self.RBar),axis=1))/numpy.expand_dims(q_max/(numpy.sum(self.p*self.deltaT/self.rho_p,axis=1)),axis=1)

    @property
    def Prandtl_proxy4(self):
        return (self.normed_m_momentum_flux_psi_unit_vector/numpy.max(self.normed_conductive_heat_flux_psi_unit_vector,axis=0))*(self.TBar/(self.mBar*self.masses))**(1.0/2.0)/numpy.max(self.deltaN,axis=0)

    @property
    def Prandtl_proxy4_1(self):
        chi_pi = self.normed_m_momentum_flux_psi_unit_vector/(self.RBar*self.normed_n_FSA_toroidal_mass_flow*self.deltaN/self.rho_p)
        chi_q = self.normed_conductive_heat_flux_psi_unit_vector/(self.p*self.deltaT/self.rho_p)
        return chi_pi/chi_q

    @property
    def Prandtl_proxy4_2(self):
        ion_index=0
        return self.normed_m_momentum_flux_psi_unit_vector*numpy.expand_dims((1/self.normed_conductive_heat_flux_psi_unit_vector[:,ion_index])*(self.p[:,ion_index])/(self.RBar*self.normed_n_FSA_toroidal_mass_flow[:,ion_index])*(self.deltaT[:,ion_index]/self.deltaN[:,ion_index]),axis=1)

    @property
    def Prandtl_proxy4_4(self):
        ion_index=0
        psi_index=90
        dpsi_index=1
        dmnVtdr=numpy.fabs((self.normed_n_FSA_toroidal_mass_flow[psi_index+dpsi_index,ion_index]-self.normed_n_FSA_toroidal_mass_flow[psi_index-dpsi_index,ion_index])/(self.actual_psi[psi_index+dpsi_index] - self.actual_psi[psi_index-dpsi_index])*numpy.fabs(self.FSABp[psi_index]/self.RBar))
        print "dmnVtdr: " + str(dmnVtdr)
        print "mVtn deltaN/rho_p " + str(self.deltaN[psi_index,ion_index]*self.normed_n_FSA_toroidal_mass_flow[psi_index,ion_index]/self.rho_p[psi_index,ion_index])
        print "mVtn/n dn/dr " + str((self.normed_n_FSA_toroidal_mass_flow[psi_index,ion_index]/self.n[psi_index,ion_index])*(self.n[psi_index+dpsi_index,ion_index] - self.n[psi_index-dpsi_index,ion_index])/(self.actual_psi[psi_index+dpsi_index] - self.actual_psi[psi_index-dpsi_index])*numpy.fabs(self.FSABp[psi_index]/self.RBar))
        
        return self.normed_m_momentum_flux_psi_unit_vector*((1/self.normed_conductive_heat_flux_psi_unit_vector[:,ion_index])*(self.p[:,ion_index])/(self.RBar)*(self.deltaT[:,ion_index]/(dmnVtdr*self.rho_p[:,ion_index])))[psi_index]


