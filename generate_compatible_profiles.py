#this class generates n_i, n_z and Phi for a given n_e, T_e, T_i and T_z
#so that that the orderings in PERFECT are satisfied
from __future__ import division

from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simulation
import h5py #to write profiles
import numpy #matrix things, etc
from bezier_transition import bezier_transition,derivative_bezier_transition #for smooth transition between 2 functions
import scipy.constants #constants
import scipy.optimize
import os
from nu_r import nu_r #nu_r
from coulomb_logarithms import lambda_ee as coulombLog #lambda=ln(Lambda)
from perfectProfilesFile import perfectProfiles
from perfectGlobalMultiplierProfileFile import create_globalMultiplier_of_Npsi

from set_species_param import set_species_param

from mtanh import generate_m2tanh_profile, extrapolate_m2tanh_sections, match_heat_flux_proxy, derivative_m3tanh_transition
from sinusoidal import generate_sin_profile
from generate_nonuniform_grid import generate_nonuniform_grid_Int_arctan
from global_multiplier import generate_global_multiplier

#for debugging purposes
import matplotlib.pyplot as plt

#print things
verbose=False

###############################
# GLOBAL VARIABLES
psiMinPed = None #pedestal start in psiN
psiMaxPed = None #pedestal end in psiN
psiMin = None #domain start in psiN
psiMax = None #domain end in psiN
Nspecies = None
mtanh = None
species = None
dxdpsiN_at_a = None
psiList = None
pairList = None
Npsi = None
main_index = None
e_index = None
imp_index = None
Zs = None
ms = None
#Needed to extrapolate to final constant Phi
TiHatPed = None
etaiHatPed = None
#Needed for eta_i
TiHatAft = None
dTiHatAftdpsi = None
specialeta = None
niPed = None
niCoreGrad = None
dniHatAftdpsi = None
niHatPed = None
niHatAft = None
# Needed for constant n_e case
neHatPre = None
dneHatPredpsi = None


# END GLOBAL VARIABLES
##############################

def sample_funclist(f_i,x):
    #samples a list of 1d functions f_i
    #returns a 2D array A_ij = f_i(x_j)
    if isinstance(f_i, (list, tuple, numpy.ndarray)):
        f_ij = [f(x) for f in f_i]
    else:
        #if not a list, we just sample
        f_ij = f_i(x)
    return numpy.array(f_ij)
    

def write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi,dpsi_ds,C):
    s = simul.psi
    #psi(s)
    psiN = psi(s)
    psi = psiN/simul.psiAHat
    #dpsi_ds(s)
    dhds = dpsi_ds(s)
    
    #tranpose
    THats = numpy.transpose(sample_funclist(THats,psiN))
    dTHatdss = numpy.transpose(sample_funclist(dTHatdss,psiN))
    nHats = numpy.transpose(sample_funclist(nHats,psiN))
    dnHatdss = numpy.transpose(sample_funclist(dnHatdss,psiN))
    etaHats = numpy.transpose(sample_funclist(etaHats,psiN))
    detaHatdss = numpy.transpose(sample_funclist(detaHatdss,psiN))
    PhiHat = sample_funclist(PhiHat,psiN)
    dPhiHatds = sample_funclist(dPhiHatds,psiN)

    #convert back to s derivatives
    dPhiHatds =  dPhiHatds * dhds
    dhds = dhds[:,numpy.newaxis]
    dTHatdss =dTHatdss * dhds
    dnHatdss = dnHatdss * dhds
    detaHatdss = detaHatdss * dhds

    #if we need to generate a global term multiplier, we do so here.
    if simul.inputs.useGlobalTermMultiplier == 1:
        C = numpy.transpose(sample_funclist(C,psiN))
        create_globalMultiplier_of_Npsi("globalTermMultiplier.h5",Npsi,C)  
    
    # Write profiles
    #print simul.input_dir+"/"+simul.inputs.profilesFilename
    try:
        outputfile = perfectProfiles(simul.input_dir+"/"+simul.inputs.profilesFilename)
    except IOError:
        print "#######################################################"
        print simul.input_dir+"/"+simul.inputs.profilesFilename + " already exists, cannot generate input profiles hdf5 file."
        print "#######################################################"
    else:
        if Nspecies == 1:
            outputfile.create_profiles_for_Npsi(Npsi,1,PhiHat,dPhiHatds,THats[:,[main_index]],dTHatdss[:,[main_index]],nHats[:,[main_index]],dnHatdss[:,[main_index]],etaHats[:,[main_index]],detaHatdss[:,[main_index]])
        elif Nspecies == 2:
            outputfile.create_profiles_for_Npsi(Npsi,2,PhiHat,dPhiHatds,THats[:,[main_index,e_index]],dTHatdss[:,[main_index,e_index]],nHats[:,[main_index,e_index]],dnHatdss[:,[main_index,e_index]],etaHats[:,[main_index,e_index]],detaHatdss[:,[main_index,e_index]])
        elif Nspecies == 3:
            outputfile.create_profiles_for_Npsi(Npsi,Nspecies,PhiHat,dPhiHatds,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss)
        else:
            print "#######################################################"
            print "Nspecies > 3 in profile writing function. How could this happen?"
            print "#######################################################"
 

def get_psiN(nonuniform,simul,**kwargs):
    #get new psiN coordinate from input
    if nonuniform:
        grid_type = kwargs["grid_type"]
        if grid_type == "Int_arctan":
            a = kwargs["transition_length"]
            b = kwargs["pedestal_grid_density"]
            c=1
            
            s = simul.inputs.psi
            psiAHat=simul.inputs.psiAHat
        
            # by design of non-uniform grid, s[0]=psiN[0], s[-1]=psiN[-1]
            # for the physical psiN. Not true for indices between
            if "leftshift" in kwargs.keys():
                leftshift = (psiMaxPed-psiMinPed)*kwargs["leftshift"]
            else:
                leftshift = 0

            if "rightshift" in kwargs.keys():
                rightshift = (psiMaxPed-psiMinPed)*kwargs["rightshift"]
            else:
                rightshift = 0
            
            psiN1= psiMinPed - leftshift
            psiN2= psiMaxPed + rightshift

            actual_psi,dactual_psi_ds = generate_nonuniform_grid_Int_arctan(simul.input,a,b,c,psiN1,psiN2)
        else:
            "generate_compatible_profiles: ERROR: unrecognized grid type!"
            exit(1)

        #function psiN(s) (here: psi) in terms of psi(s) (here: actual_psi)
        psi=lambda x: actual_psi(x)/psiAHat #nonuniform psiN
        dpsi_ds = lambda x: dactual_psi_ds(x)/psiAHat
    else:
        psi= lambda x: x
        dpsi_ds = lambda x:  1.0 + 0*x #trivial mapping
    return (psi,dpsi_ds)
            
def update_domain_size(simul,psiN_ped_width,midShift,psiDiamFact,leftBoundaryShift,rightBoundaryShift):
    psiMid=1-psiN_ped_width/2.0+midShift
    psiDiameter=psiDiamFact*psiN_ped_width
    
    simul.inputs.psiDiameter = psiDiameter
    simul.inputs.psiMid = psiMid
    simul.inputs.leftBoundaryShift = leftBoundaryShift
    simul.inputs.rightBoundaryShift = rightBoundaryShift

def update_species_parameters(species_list,species_filename,norm_filename,simulation):
    set_species_param(species_list,species_filename,norm_filename,simulation)
    Zs=simulation.inputs.charges
    if not isinstance(Zs, (list, tuple, numpy.ndarray)):
        Zs=numpy.array([Zs])
    ms=simulation.inputs.masses
    if not isinstance(Zs, (list, tuple, numpy.ndarray)):
        ms=numpy.array([ms])
    return (Zs,ms)

def update_Delta_omega(simul):
    #print "Delta before:" + str(simul.Delta)
    if simul.units=="SI":
        Delta = simul.mBar*simul.vBar/(simul.eBar*simul.BBar*simul.RBar)
        omega = simul.ePhiBar/(simul.eBar*simul.BBar*simul.RBar*simul.vBar)
        simul.inputs.Delta=Delta
        simul.inputs.omega=omega
    else:
        print "Only SI units are currently supported"
    return (Delta,omega)
    #print "Delta after:" + str(Delta)

def update_nu_r(simul,n,T):
    if simul.units=="SI":
        logLambda=coulombLog(n,T) #n,T not used beyond this point
        print "logLambda: "+str(logLambda)
        nur = nu_r(simul.RBar,simul.nBar,simul.TBar,logLambda)
        simul.inputs.nu_r=nur
        return nur
    else:
        print "Only SI units are currently supported"

def generate_ne_profile(simul,**kwargs):
    global neHatPre
    global dneHatPredpsi
    # for use in d-He scan where the electron profile is fixed
    if "nScale_"+species[e_index] in kwargs.keys():
        neScale=kwargs["nScale_"+species[e_index]]
    else:
        neScale=1.0
    nePed=neScale*kwargs["nped_"+species[e_index]]
    neCoreGrad=neScale*kwargs["dnCoredx_"+species[e_index]]*dxdpsiN_at_a
    nepedGrad=neScale*kwargs["dnpeddx_"+species[e_index]]*dxdpsiN_at_a
    neSOLGrad=neScale*kwargs["dnSOLdx_"+species[e_index]]*dxdpsiN_at_a

    if mtanh:
        (neHat,dneHatds) = generate_m2tanh_profile(nePed,neCoreGrad,nepedGrad,neSOLGrad,psiMaxPed-psiMinPed,psiMinPed)
        (core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(neHat,dneHatds,psiMin,psiMinPed,psiMaxPed,psiMax)
        neHatPre = core
        dneHatPredpsi = ddx_core
    else:
        neHatPre =(lambda psiN: (nePed-neCoreGrad*(psiMinPed-psiMin) +  neCoreGrad* (psiN-psiMin)))
        neHatPed =(lambda psiN: (nePed +  nepedGrad* (psiN-psiMinPed)))
        neHatAft =(lambda psiN: (nePed+nepedGrad*(psiMaxPed-psiMinPed) +  neSOLGrad* (psiN-psiMaxPed)))
    
        dneHatPredpsi = (lambda psiN:  neCoreGrad + psiN*0)
        dneHatPeddpsi = (lambda psiN:  nepedGrad + psiN*0)
        dneHatAftdpsi = (lambda psiN:  neSOLGrad + psiN*0)
    
        nelist=[neHatPre,neHatPed,neHatAft]
        dnedpsiList = [dneHatPredpsi,dneHatPeddpsi,dneHatAftdpsi] 
        (neHat,dneHatdpsi)=derivative_bezier_transition(nelist,dnedpsiList,psiList,pairList)
        dneHatds = lambda x : dneHatdpsi(x)

    
    
    return (neHat,dneHatds)
        
def generate_ni_profile(**kwargs):
    # for use in d-He scan where the electron profile varies
    global niPed
    global niCoreGrad
    global niHatPed
    global niHatAft 
    global dniHatAftdpsi
    if "nScale_"+species[main_index] in kwargs.keys():
        niScale=kwargs["nScale_"+species[main_index]]
    else:
        niScale=1.0
    niPed=niScale*kwargs["nped_"+species[main_index]]
    niCoreGrad=niScale*kwargs["dnCoredx_"+species[main_index]]*dxdpsiN_at_a
    nipedGrad=niScale*kwargs["dnpeddx_"+species[main_index]]*dxdpsiN_at_a
    niSOLGrad=niScale*kwargs["dnSOLdx_"+species[main_index]]*dxdpsiN_at_a
    

    # generate n_i
    if mtanh:
        (niHat,dniHatds) = generate_m2tanh_profile(niPed,niCoreGrad,nipedGrad,niSOLGrad,psiMaxPed-psiMinPed,psiMinPed)
        (core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(niHat,dniHatds,psiMin,psiMinPed,psiMaxPed,psiMax)
        niHatPed = ped
        niHatAft = sol
        dniHatAftdpsi = ddx_sol
    else:
        niHatPre =(lambda psiN: (niPed + niCoreGrad*(psiN-psiMinPed)))
        niHatPed =(lambda psiN: (niPed +  nipedGrad* (psiN-psiMinPed)))
        niHatAft =(lambda psiN: (niPed+nipedGrad*(psiMaxPed-psiMinPed) +  niSOLGrad* (psiN-psiMaxPed)))

        dniHatPredpsi = (lambda psiN:  niCoreGrad + 0*psiN)
        dniHatPeddpsi = (lambda psiN:  nipedGrad + 0*psiN)
        dniHatAftdpsi = (lambda psiN:  niSOLGrad + 0*psiN)

        nilist=[niHatPre,niHatPed,niHatAft]
        dnidpsiList = [dniHatPredpsi,dniHatPeddpsi,dniHatAftdpsi] 
        
        (niHat,dniHatds)=derivative_bezier_transition(nilist,dnidpsiList,psiList,pairList)
    
    return (niHat,dniHatds)

def generate_T_profiles_sameflux(nt,nb,breakpoint_shift,T_transition_length,**kwargs):
    global TiHatPed
    global TiHatAft
    global dTiHatAftdpsi
    
    #THats = numpy.zeros((Nspecies,Npsi))
    #dTHatdpsis = numpy.zeros((Nspecies,Npsi))
    THats = [None]*Nspecies
    dTHatdpsis = [None]*Nspecies
    dTHatdss = [None]*Nspecies
    
    TScale=[0]*Nspecies
    Tpeds=[0]*Nspecies
    TCoreGrads=[0]*Nspecies
    TpedGrads=[0]*Nspecies
    TSOLGrads=[0]*Nspecies

    #Pre,Ped,Aft will contain functions describe T in core, ped and outwards.
    THatPre=[0]*Nspecies
    THatPed=[0]*Nspecies
    THatAft=[0]*Nspecies

    dTHatPredpsi=[0]*Nspecies
    dTHatPeddpsi=[0]*Nspecies
    dTHatAftdpsi=[0]*Nspecies
    
    #reading T related parameters
    for i in range(Nspecies):
        if "TScale_"+species[i] in kwargs.keys():
            TScale[i]=kwargs["TScale_"+species[i]]
        else:
            TScale[i]=1.0
        Tpeds[i]=TScale[i]*kwargs["Tped_"+species[i]]
        TCoreGrads[i]=TScale[i]*kwargs["dTCoredx_"+species[i]]*dxdpsiN_at_a
        TpedGrads[i]=TScale[i]*kwargs["dTpeddx_"+species[i]]*dxdpsiN_at_a
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        
    TScale=numpy.array(TScale)
    Tpeds=numpy.array(Tpeds)
    TCoreGrads=numpy.array(TCoreGrads)
    TSOLGrads=numpy.array(TSOLGrads)
    TpedGrads=numpy.array(TpedGrads)
    
    if mtanh:
        transition_length = T_transition_length*(psiMaxPed-psiMinPed)
        # we always want the transition to end at the same place and start at the same density...
        #... thus we shift the start when we change the transition length
        transition_start = psiMinPed - (T_transition_length -1)*(psiMaxPed-psiMinPed)
        print "!!!!!!!!!!!!!!!!!!"
        print psiMin
        print transition_start
        
        #... and change the Tped
        Tped = Tpeds[main_index] - TpedGrads[main_index]*(T_transition_length -1)*(psiMaxPed-psiMinPed)
        print Tped
        
        (THats[main_index],dTHatdss[main_index]) = match_heat_flux_proxy(Tped,TpedGrads[main_index],TpedGrads[main_index],TSOLGrads[main_index],transition_length,transition_start,nt,nb,psiMin,psiMax)
        (core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(THats[main_index],dTHatdss[main_index],psiMin,psiMinPed,psiMaxPed,psiMax)
        
        TiHatPed = ped
        TiHatAft = sol
        dTiHatAftdpsi = ddx_sol
        
        if Nspecies >= 2:
            (THats[e_index],dTHatdss[e_index]) = generate_m2tanh_profile(Tpeds[e_index],TCoreGrads[e_index],TpedGrads[e_index],TSOLGrads[e_index],psiMaxPed-psiMinPed,psiMinPed)
    else:
        breakpoint=psiMinPed + breakpoint_shift
        Tp=Tpeds[main_index] + TCoreGrads[main_index] *(breakpoint-psiMinPed)

        if (TCoreGrads[main_index] != TSOLGrads[main_index]) and (TPedGrads[main_index] != TSOLGrads[main_index]):
            print "generate_compatible_profiles: Warning: sameflux option uses the Core T gradients of the main species everywhere"
        dTdpsi=TCoreGrads[main_index]
        d1=breakpoint - psiMin
        d2=psiMax - breakpoint

        f=lambda x : x*(Tp-d1*x)**(3.0/2.0) - (nb*1.0/nt)*(Tp + dTdpsi*d2)**(3.0/2.0)*dTdpsi
        dTdpsiTop=scipy.optimize.fsolve(f,0)[0] #[0] since returns an numpy.ndarray
        dTdpsiBot=dTdpsi

        THatPre[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiTop*(psiN-breakpoint)))
        THatPed[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
        THatAft[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
        dTHatPredpsi[main_index] =(lambda psiN: dTdpsiTop + 0*psiN)
        dTHatPeddpsi[main_index] =(lambda psiN: dTdpsiBot + 0*psiN)
        dTHatAftdpsi[main_index] =(lambda psiN: dTdpsiBot + 0*psiN)


        Tlist=[THatPre[main_index],THatPed[main_index]]
        dTdpsilist=[dTHatPredpsi[main_index],dTHatPeddpsi[main_index]]

        (THats[main_index],dTHatdpsis[main_index])=derivative_bezier_transition(Tlist,dTdpsilist,[breakpoint],pairList[:-1])    
        dTHatdss[main_index]= lambda x : dTHatdpsis[main_index](x)

        #T2=simul.TBar*bezier_transition(Tlist,[breakpoint],pairList[:-1],numpy.array([psiMidPed]))[0]

        TiHatPed = THatPed[main_index]
        TiHatAft = THatAft[main_index]
        dTiHatAftdpsi = dTHatAftdpsi[main_index]
        
        if Nspecies >= 2:
            #T_e:
            THatPre[e_index] =(lambda psiN: (Tpeds[e_index] + TCoreGrads[e_index]*(psiN-psiMinPed)))
            THatPed[e_index] =(lambda psiN: (Tpeds[e_index] + TpedGrads[e_index]*(psiN-psiMinPed)))
            THatAft[e_index] =(lambda psiN: (Tpeds[e_index] + TpedGrads[e_index]*(psiMaxPed-psiMinPed) + TSOLGrads[e_index]*(psiN-psiMaxPed)))
            dTHatPredpsi[e_index] = (lambda psiN: TCoreGrads[e_index] + 0*psiN)
            dTHatPeddpsi[e_index] = (lambda psiN: TpedGrads[e_index] + 0*psiN)
            dTHatAftdpsi[e_index] = (lambda psiN: TSOLGrads[e_index] + 0*psiN)
            Tlist=[THatPre[e_index],THatPed[e_index],THatAft[e_index]]
            dTdpsiList = [dTHatPredpsi[e_index],dTHatPeddpsi[e_index],dTHatAftdpsi[e_index]]

            (THats[e_index],dTHatdpsis[e_index])=derivative_bezier_transition(Tlist,dTdpsiList,psiList,pairList)
            dTHatdss[e_index]=lambda x : dTHatdpsis[e_index](x)

    if Nspecies == 3:
        THats[imp_index]=THats[main_index]
        #dTHatdpsis[imp_index]=dTHatdpsis[main_index]
        dTHatdss[imp_index]=dTHatdss[main_index]


    
    
    return (THats,dTHatdss)
    
def generate_T_profiles(**kwargs):
    global TiHatPed
    global TiHatAft
    global dTiHatAftdpsi
    
    #THats = numpy.zeros((Nspecies,Npsi))
    #dTHatdpsis = numpy.zeros((Nspecies,Npsi))
    THats = [None]*Nspecies
    dTHatdpsis = [None]*Nspecies
    dTHatdss = [None]*Nspecies
    
    TScale=[0]*Nspecies
    Tpeds=[0]*Nspecies
    TCoreGrads=[0]*Nspecies
    TpedGrads=[0]*Nspecies
    TSOLGrads=[0]*Nspecies
    
    #Pre,Ped,Aft will contain functions describe T in core, ped and outwards.
    THatPre=[0]*Nspecies
    THatPed=[0]*Nspecies
    THatAft=[0]*Nspecies

    dTHatPredpsi=[0]*Nspecies
    dTHatPeddpsi=[0]*Nspecies
    dTHatAftdpsi=[0]*Nspecies
    
    #reading T related parameters
    for i in range(Nspecies):
        if "TScale_"+species[i] in kwargs.keys():
            TScale[i]=kwargs["TScale_"+species[i]]
        else:
            TScale[i]=1.0
        Tpeds[i]=TScale[i]*kwargs["Tped_"+species[i]]
        TCoreGrads[i]=TScale[i]*kwargs["dTCoredx_"+species[i]]*dxdpsiN_at_a
        TpedGrads[i]=TScale[i]*kwargs["dTpeddx_"+species[i]]*dxdpsiN_at_a
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        
    TScale=numpy.array(TScale)
    Tpeds=numpy.array(Tpeds)
    TCoreGrads=numpy.array(TCoreGrads)
    TSOLGrads=numpy.array(TSOLGrads)
    TpedGrads=numpy.array(TpedGrads)
    
    if mtanh:
        for i in range(Nspecies):
            (THats[i],dTHatdss[i]) = generate_m2tanh_profile(Tpeds[i],TCoreGrads[i],TpedGrads[i],TSOLGrads[i],psiMaxPed-psiMinPed,psiMinPed)
        (core,ped,sol,ddx_core,ddx_ped,ddx_sol) = extrapolate_m2tanh_sections(THats[main_index],dTHatdss[main_index],psiMin,psiMinPed,psiMaxPed,psiMax)

        TiHatPed = ped
        TiHatAft = sol
        dTiHatAftdpsi = ddx_sol
    else:
        #Generating functions from the parameters
        for i in range(Nspecies):
            THatPre[i] =(lambda psiN,i=i: (Tpeds[i] + TCoreGrads[i]*(psiN-psiMinPed)))
            THatPed[i] =(lambda psiN,i=i: (Tpeds[i] + TpedGrads[i]*(psiN-psiMinPed)))
            THatAft[i] =(lambda psiN,i=i: (Tpeds[i] + TpedGrads[i]*(psiMaxPed-psiMinPed) + TSOLGrads[i]*(psiN-psiMaxPed)))
            dTHatPredpsi[i] = (lambda psiN,i=i: TCoreGrads[i])
            dTHatPeddpsi[i] = (lambda psiN,i=i: TpedGrads[i])
            dTHatAftdpsi[i] = (lambda psiN,i=i: TSOLGrads[i])
            Tlist=[THatPre[i],THatPed[i],THatAft[i]]
            dTdpsiList = [dTHatPredpsi[i],dTHatPeddpsi[i],dTHatAftdpsi[i]]
            
            (THats[i],dTHatdpsis[i])=derivative_bezier_transition(Tlist,dTdpsiList,psiList,pairList)
            dTHatdss[i]=lambda x,i=i: dTHatdpsis[i](x)
        
    
        TiHatPed = THatPed[main_index]
        TiHatAft = THatAft[main_index]
        dTiHatAftdpsi = dTHatAftdpsi[main_index]
    
    return (THats,dTHatdss)

def generate_etai_profile(Delta,omega,**kwargs):    
    etaiHatPre =(lambda psiN: (niPed+  niCoreGrad* (psiN-psiMinPed)))
    detaiHatPredpsi =(lambda psiN: niCoreGrad + 0*psiN)

    if mtanh:
        #PhiTopPoint=psiMax/2.0 + psiMaxPed/2.0
        PhiTopPoint=psiMaxPed
        width = (psiMaxPed - psiMinPed)
        
        prefactor=1.0
        PhiTop=prefactor*(TiHatPed(PhiTopPoint)/Zs[main_index])*numpy.log(etaiHatPre(PhiTopPoint)*1.0/niHatAft(PhiTopPoint))*(2.0*omega/Delta)

        etaiHatAft =(lambda psiN: niHatAft(psiN)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))
        detaiHatAftdpsi = (lambda psiN: (dniHatAftdpsi(psiN) - niHatAft(psiN)*PhiTop*Zs[main_index]*dTiHatAftdpsi(psiN)/(TiHatAft(psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))

       
        (etaiHat,detaiHatds) = derivative_m3tanh_transition([etaiHatPre,etaiHatAft],[detaiHatPredpsi,detaiHatAftdpsi],psiMaxPed,width)

        #DEBUG PLOT
        if verbose:
            plot_psi = numpy.linspace(psiMin,psiMax)
            plt.hold(True)
            plt.plot(plot_psi,etaiHat(plot_psi))
            plt.plot(plot_psi,etaiHatPre(plot_psi))
            plt.plot(plot_psi,etaiHatAft(plot_psi))
            plt.show()
    else:
        if specialeta:
            etaiHatPed =(lambda psiN: (niPed+2*niCoreGrad* (psiN-psiMinPed)))
            detaiHatPeddpsi =(lambda psiN: 2*niCoreGrad + 0*psiN)
        else:
            etaiHatPed =(lambda psiN: (niPed+  niCoreGrad* (psiN-psiMinPed)))
            detaiHatPeddpsi =(lambda psiN: niCoreGrad + 0*psiN)

        
        PhiTopPoint=psiMaxPed

        prefactor=1.0
        PhiTop=prefactor*(TiHatPed(PhiTopPoint)/Zs[main_index])*numpy.log(etaiHatPed(PhiTopPoint)*1.0/niHatPed(PhiTopPoint))*(2.0*omega/Delta)

        etaiHatAft =(lambda psiN: niHatAft(psiN)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))
        detaiHatAftdpsi = (lambda psiN: (dniHatAftdpsi(psiN) - niHatAft(psiN)*PhiTop*Zs[main_index]*dTiHatAftdpsi(psiN)/(TiHatAft(psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))
    
        etailist=[etaiHatPre,etaiHatPed,etaiHatAft]
        detaidpsilist=[detaiHatPredpsi,detaiHatPeddpsi,detaiHatAftdpsi]

        (etaiHat,detaiHatds) =derivative_bezier_transition(etailist,detaidpsilist,psiList,pairList)

    return (etaiHat,detaiHatds)

def generate_etaz_profile(Delta,omega,**kwargs):
    imp_conc=kwargs["imp_conc"]
    if specialeta==True:
        gradScale=4.0725
        etazHatInner=(lambda psiN: imp_conc*(niPed+niCoreGrad*(psiN-psiMinPed)))
        etazHatMiddle=(lambda psiN: imp_conc*(niPed-gradScale*niCoreGrad*(psiN-psiMinPed)))
        etazHatOuter=(lambda psiN: imp_conc*(niPed-gradScale*niCoreGrad*(psiMaxPed-psiMinPed) +  niCoreGrad* (psiN-psiMaxPed)))

        detazHatInnerdpsi=(lambda psiN: imp_conc*niCoreGrad)
        detazHatMiddledpsi=(lambda psiN: -imp_conc*gradScale*niCoreGrad)
        detazHatOuterdpsi=(lambda psiN: imp_conc*niCoreGrad)

        etailist=[etazHatInner,etazHatMiddle,etazHatOuter]
        detaidpsilist = [detazHatInnerdpsi,detazHatMiddledpsi,detazHatOuterdpsi]

        (etazHat,detazHatdpsi) = derivative_bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:])
        detazHatds = lambda x : detazHatdpsi(x)
        
    else:
        etazHat=lambda x : imp_conc*(niPed+niCoreGrad*(x-psiMinPed))
        detazHatds = lambda x: imp_conc*niCoreGrad + 0*x #add 0*x to make it return ndarray when x is one.
    return (etazHat,detazHatds)

def generate_Phi_profile(etaiHat,niHat,TiHat,detaiHatdpsi,dniHatdpsi,dTiHatdpsi,Zi,Delta,omega):
    PhiHat=lambda x : numpy.log(etaiHat(x)/niHat(x))*TiHat(x)/Zi*Delta/(2*omega)
    dPhiHatdpsi=lambda x : (-etaiHat(x)*dniHatdpsi(x)/niHat(x)**2 + detaiHatdpsi(x)/niHat(x))*TiHat(x)*niHat(x)/etaiHat(x) + numpy.log(etaiHat(x)/niHat(x))*dTiHatdpsi(x)
    return (PhiHat,dPhiHatdpsi)

def generate_eta_from_n_Phi_profile(nHat,PhiHat,THat,dnHatdpsi,dPhiHatdpsi,dTHatdpsi,Z,Delta,omega):
    etaHat=lambda x: nHat(x)*numpy.exp((Z*omega*2/Delta)*PhiHat(x)/THat(x))
    detaHatdpsi=lambda x: (dnHatdpsi(x) + (2*omega/Delta)*(nHat(x)*Z/THat(x))*(dPhiHatdpsi(x)-PhiHat(x)*dTHatdpsi(x)/THat(x)))*numpy.exp((Z*omega*2/Delta)*PhiHat(x)/THat(x))
    return (etaHat,detaHatdpsi)

def generate_n_from_eta_Phi_profile(etaHat,PhiHat,THat,detaHatdpsi,dPhiHatdpsi,dTHatdpsi,Z,Delta,omega):
    nHat=lambda x : etaHat(x)*numpy.exp(-(Z*omega*2/Delta)*PhiHat(x)/THat(x))
    dnHatdpsi=lambda x : (detaHatdpsi(x) - (2*omega/Delta)*(etaHat(x)*Z/THat(x))*(dPhiHatdpsi(x)-PhiHat(x)*dTHatdpsi(x)/THat(x)))*numpy.exp(-(Z*omega*2/Delta)*PhiHat(x)/THat(x))
    return (nHat,dnHatdpsi)

def generate_n_from_eta_X_profile(etaHat,X,detaHatds,dXds,Z):
    nHat=lambda x : etaHat(x)*X(x)**Z
    dnHatdpsi=lambda x : detaHatds(x)*X(x)**Z + etaHat(x)*Z*X(x)**(Z-1)*dXds(x)
    return (nHat,dnHatdpsi)

def generate_compatible_profiles(simul,xwidth,nonuniform=False,sameflux=False,oldsameflux=False,sameeta=False,samefluxshift=0,specialEta=False,psiDiamFact=5,transitionFact=0.2,dxdpsiN=1,midShift=0,upShift_denom=3.0,m2tanh=False,mode="const_Phi",leftBoundaryShift=0.0,rightBoundaryShift=0.0,T_transition_length=1.0,**kwargs):
    #NOTE: uses the dx/dpsi at minor radius for both core and SOL, which assumes that
    #simulated region is small enough that it doesn't change much.
    if mode == "periodic":
        generate_periodic_profiles(simul,xwidth,nonuniform,dxdpsiN,psiDiamFact,midShift,upShift_denom,leftBoundaryShift,rightBoundaryShift,**kwargs)
        return 
    
    e=scipy.constants.e
    log=numpy.log
    exp=numpy.exp
    sqrt=numpy.sqrt
    
    global psiMinPed
    global psiMaxPed
    global psiMin
    global psiMax
    global Nspecies
    global species
    global dxdpsiN_at_a
    global psiList
    global pairList
    global Npsi
    global main_index
    global e_index
    global imp_index
    global specialeta
    global Zs
    global ms
    global mtanh
    mtanh = m2tanh
    specialeta = specialEta

    species=simul.species
    Nspecies=len(species)

    #intepret possible keywords in the keyword-argument dict kwargs
    if Nspecies == 1:
        main_index=kwargs["mI"]
    
    elif Nspecies == 2:
        main_index=kwargs["mI"]
        e_index=kwargs["eI"]
    
    elif Nspecies == 3:
        main_index=kwargs["mI"]
        imp_index=kwargs["zI"]
        e_index=kwargs["eI"]
    else:
        print "Script is too poor to handle an arbitrary number of species!"
        return -1

    if dxdpsiN == 1:
        print "Intepreting gradients as being given in psi_N."
    else:
        print "Intepreting gradients as being given in x space"
        print "Will use function specified by dx_dpsi to convert" 

    dxdpsiN_at_a = dxdpsiN
    x_ped_width=xwidth
    psiN_ped_width=x_ped_width/dxdpsiN_at_a
    
    #calculate new psiMid and diameter for this profile
    midShift=midShift*psiN_ped_width

    #2015-12-16: -psiN_ped_width/3.0 was the old default.
    upShift=-psiN_ped_width/upShift_denom

    
    update_domain_size(simul,psiN_ped_width,midShift,psiDiamFact,leftBoundaryShift,rightBoundaryShift)
    (Delta,omega) = update_Delta_omega(simul)

    #these things do not depend on whether the grid is uniform or nonuniform
    # if intepreted as midpoint of physical psiN coordinate, etc.
    # this is by design, since s[0] and s[-1] coincide with psiN start and stop
    Npsi=simul.inputs.Npsi
    psiMid = simul.inputs.psiMid
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax
    
    #start and stop pedestal in physical psiN
    psiMinPed=psiMid-psiN_ped_width/2.0
    psiMaxPed=psiMid+psiN_ped_width/2.0+upShift
    psiMidPed=(psiMinPed+psiMaxPed)/2.0
    
    
    #list of psi where our profiles change slope
    psiList=[psiMinPed,psiMaxPed]
    print "pedestal start,stop: " + str(psiList)


    if specialEta and (mode == "const_ne"):
        print "generate_compatible_profiles: WARNING: specialEta and mode == const_ne is incompatible. specialEta will be ignored"

    # not too relevant for mtanh
    offset=(psiMaxPed-psiMinPed)*transitionFact
    pairList=[[offset,offset],[offset,offset]]
        
    #allocate arrays
    THats = [None]*Nspecies
    dTHatdss = [None]*Nspecies
    nHats = [None]*Nspecies
    dnHatdss = [None]*Nspecies
    etaHats = [None]*Nspecies
    detaHatdss = [None]*Nspecies

    #this modifies the charge and mass of the species in the simulation to match species file, and reads the new parameters back
    species_filename= os.path.join(os.path.dirname(__file__), 'species_database.namelist')
    (Zs,ms)=update_species_parameters(simul.species,species_filename,simul.norm,simul) 

    if mode == "const_Phi":
        #generate n_i profile
        (nHats[main_index],dnHatdss[main_index]) = generate_ni_profile(**kwargs)
        
    elif mode == "const_ne":
        #create n_e:
        (nHats[e_index],dnHatdss[e_index]) = generate_ne_profile(simul,**kwargs)
        # generate ion etas
        c=float(kwargs["imp_conc"])
        etaHats[main_index] = lambda x: neHatPre(x)/(Zs[main_index]*(1+2*c))
        detaHatdss[main_index] =  lambda x: dneHatPredpsi(x)/(Zs[main_index]*(1+2*c))
        etaHats[imp_index] = lambda x: c*etaHats[main_index](x)
        detaHatdss[imp_index] = lambda x: c*detaHatdss[main_index](x)

        # generate ion densities
        X = lambda x: sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Zs[main_index]*etaHats[main_index](x)))-1/(4*c)
        dXds = lambda x: 1/(sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Zs[main_index]*etaHats[main_index](x)))*4*Zs[main_index]*c*etaHats[main_index](x))*(dnHatdss[e_index](x) - nHats[e_index](x)*(detaHatdss[main_index](x)/etaHats[main_index](x)))

        (nHats[main_index],dnHatdss[main_index]) = generate_n_from_eta_X_profile(etaHats[main_index],X,detaHatdss[main_index],dXds,Zs[main_index])

        (nHats[imp_index],dnHatdss[imp_index]) = generate_n_from_eta_X_profile(etaHats[imp_index],X,detaHatdss[imp_index],dXds,Zs[imp_index])

    
        
    if (sameflux==True):
        nt=nHats[main_index](psiMin)
        nb=nHats[main_index](psiMax)
        (THats,dTHatdss) = generate_T_profiles_sameflux(nt,nb,samefluxshift,T_transition_length,**kwargs)
    elif (oldsameflux==True):
        #this setting uses the old ni profile generation to generatea  temporary
        # n_i profile to generate sameflux T profiles that are the same as in
        # the constant_Phi case.
        # Useful for comparing const_Phi and const_ne with the same T
        #NOTE: need inputs specifying n_i to generate the temporary profile
        # does not really match the heat flux proxy for the actual n_i profiles
        # used in the const n_e case.
        (nH,dnHds) = generate_ni_profile(**kwargs)
        nt=nH(psiMin)
        nb=nH(psiMax)
        (THats,dTHatdss) = generate_T_profiles_sameflux(nt,nb,samefluxshift,T_transition_length,**kwargs)
    else:
        (THats,dTHatdss) = generate_T_profiles(**kwargs)

    #with n_i and T_i generated, we can evaluate logLambda at a suitable point
    #which is done to calculate reference collisionality nu_r
    T=simul.TBar*THats[main_index](psiMidPed)
    n=simul.nBar*nHats[main_index](psiMidPed)
    print "T: "+str(T)
    print "n: "+str(n)
    nu_r=update_nu_r(simul,n,T)
    print "nu_r: " + str(nu_r)

    if mode == "const_Phi":
        (etaHats[main_index],detaHatdss[main_index]) = generate_etai_profile(Delta,omega,**kwargs)

        (PhiHat,dPhiHatds) = generate_Phi_profile(etaHats[main_index],nHats[main_index],THats[main_index],detaHatdss[main_index],dnHatdss[main_index],dTHatdss[main_index],Zs[main_index],Delta,omega)

        #if Phi=0 at top of the pedestal, this gives the top of the n_z pedestal.
        #To make n_z and n_i same at points, those points should satisfy
        # eta_z=n_i (n_i/eta_i)^(-[Zz/Zi] Ti/Tz)
        #etaHats[imp_index] = 0.01*nHats[main_index][psiMinPedIndex]
        if (Nspecies>2) and (sameeta==False):
            if sameeta==True:
                etaHats[imp_index]=etaHats[main_index]
                detaHatdss[imp_index]=detaHatdss[main_index]
            else:
                (etaHats[imp_index],detaHatdss[imp_index]) = generate_etaz_profile(Delta,omega,**kwargs)

        if (Nspecies > 2):
            (nHats[imp_index],dnHatdss[imp_index]) = generate_n_from_eta_Phi_profile(etaHats[imp_index],PhiHat,THats[imp_index],detaHatdss[imp_index],dPhiHatds,dTHatdss[imp_index],Zs[imp_index],Delta,omega)

            
        if (Nspecies > 1):
            #n_e from quasi-neutrality
            if (Nspecies == 2):
                nHats[e_index]=lambda x : Zs[main_index]*nHats[main_index](x)
                dnHatdss[e_index]= lambda x : Zs[main_index]*dnHatdss[main_index](x)
            elif (Nspecies == 3):
                nHats[e_index]=lambda x : Zs[imp_index]*nHats[imp_index](x) + Zs[main_index]*nHats[main_index](x)
                dnHatdss[e_index]=lambda x: Zs[imp_index]*dnHatdss[imp_index](x) + Zs[main_index]*dnHatdss[main_index](x)

    elif mode == "const_ne":
        # generate Phi
        PhiHat = lambda x: -(Delta/(Zs[main_index]*2*omega))*THats[main_index](x)*log(sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Zs[main_index]*etaHats[main_index](x)))-1/(4*c))
        dPhiHatds = lambda x: -(Delta/(Zs[main_index]*2*omega))*(dTHatdss[main_index](x)*log(X(x)) + THats[main_index](x)*dXds(x)/X(x))

        
    if Nspecies > 1:
        #eta_e from n_e
        (etaHats[e_index], detaHatdss[e_index]) = generate_eta_from_n_Phi_profile(nHats[e_index],PhiHat,THats[e_index],dnHatdss[e_index],dPhiHatds,dTHatdss[e_index],Zs[e_index],Delta,omega)

    #get psiN, dpsiN_ds that maps from uniform grid to physical psi
    (psi,dpsi_ds) = get_psiN(nonuniform,simul,**kwargs)
    if simul.inputs.useGlobalTermMultiplier == 1:
        multiplier_delta_a = kwargs["multiplier_transition_shift"]
        multiplier_delta_b = multiplier_delta_a
        multiplier_c = kwargs["multiplier_edge_value"]
        multiplier_Delta = kwargs["multiplier_transition_length"]
        #multiplier_delta_a = 0.04
        #multiplier_delta_b = 0.04
        #multiplier_c =0.1
        #multiplier_Delta = 1/200.0
        C=generate_global_multiplier(psiMin,psiMax,Delta=multiplier_Delta,delta_a=multiplier_delta_a,delta_b=multiplier_delta_b,c=multiplier_c)
    else:
        C=lambda x: 1.0 + 0*x
    
    
    #sample profiles at physical psi
    #rescale derivatives from physical psi to internal uniform grid
    write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi,dpsi_ds,C)

##########################################################
# generate_compatible_profiles_constant_ne
#########################################################

def generate_compatible_profiles_constant_ne(simul,xwidth,nonuniform=False,sameflux=False,oldsameflux=False,sameeta=False,samefluxshift=0,specialEta=False,psiDiamFact=5,transitionFact=0.2,dxdpsiN=1,midShift=0,upShift_denom=3.0,m2tanh=False,T_transition_length=1.0,**kwargs):
    generate_compatible_profiles(simul,xwidth,nonuniform,sameflux,oldsameflux,sameeta,samefluxshift,specialEta,psiDiamFact,transitionFact,dxdpsiN,midShift,upShift_denom,m2tanh,mode="const_ne",T_transition_length=1.0,**kwargs)

def generate_periodic_profiles(simul,xwidth,nonuniform,dxdpsiN,psiDiamFact,midShift,upShift_denom,leftBoundaryShift,rightBoundaryShift,**kwargs):
    global psiMin
    global psiMax
    global Nspecies
    global main_index
    global Npsi

    species=simul.species
    Nspecies=len(species)
    if Nspecies == 1:
        main_index=kwargs["mI"]
    else:
        print "Error: periodic only does one species right now!"
        return -1
    
    if dxdpsiN == 1:
        print "Intepreting gradients as being given in psi_N."
    else:
        print "Intepreting gradients as being given in x space"
        print "Will use function specified by dx_dpsi to convert" 
    
    dxdpsiN_at_a = dxdpsiN
    x_ped_width=xwidth
    psiN_ped_width=x_ped_width/dxdpsiN_at_a

    #old domain calculation
    #calculate new psiMid and diameter for this profile
    midShift=midShift*psiN_ped_width

    #2015-12-16: -psiN_ped_width/3.0 was the old default.
    upShift=-psiN_ped_width/upShift_denom
    
    update_domain_size(simul,psiN_ped_width,midShift,psiDiamFact,leftBoundaryShift,rightBoundaryShift)
    (Delta,omega) = update_Delta_omega(simul)
    psiMid = simul.inputs.psiMid
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax

    psiMinPed=psiMid-psiN_ped_width/2.0
    psiMaxPed=psiMid+psiN_ped_width/2.0+upShift

    # update simulation domain to in fact only include pedestal and
    # a mirror image of the pedestal
    simul.inputs.psiDiameter = 2*(psiMaxPed - psiMinPed)
    simul.inputs.psiMid = psiMaxPed
    simul.inputs.leftBoundaryShift = 0.0
    simul.inputs.rightBoundaryShift = 0.0

    psiMid = simul.inputs.psiMid
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax
    Npsi=simul.inputs.Npsi
    print Npsi
    print psiMin
    print psiMax

    # update species parameters
    species_filename= os.path.join(os.path.dirname(__file__), 'species_database.namelist')
    (Zs,ms)=update_species_parameters(simul.species,species_filename,simul.norm,simul) 

    # generate profiles
    if "nScale_"+species[main_index] in kwargs.keys():
        niScale=kwargs["nScale_"+species[main_index]]
    else:
        niScale=1.0
    niPed=niScale*kwargs["nped_"+species[main_index]]
    niPedGrad=niScale*kwargs["dnpeddx_"+species[main_index]]*dxdpsiN_at_a
    (niHat,dniHatdpsi) = generate_sin_profile(niPed,niPedGrad,psiMinPed,psiMaxPed)

    if "TScale_"+species[main_index] in kwargs.keys():
        TiScale=kwargs["TScale_"+species[main_index]]
    else:
        TiScale = 1.0
    TiPed = TiScale*kwargs["Tped_"+species[main_index]]
    TiPedGrad = TiScale*kwargs["dTpeddx_"+species[main_index]]*dxdpsiN_at_a
    (TiHat,dTiHatdpsi) = generate_sin_profile(TiPed,TiPedGrad,psiMinPed,psiMaxPed)

    PhiPed = 0.0
    An=-niPedGrad*(psiMaxPed - psiMinPed)
    AT=-TiPedGrad*(psiMaxPed - psiMinPed)
    k = 2*numpy.pi/(psiMaxPed - psiMinPed)
    PhiPedGrad = (Delta/(2*omega))*k*An/(2*Zs[main_index])*((TiPed - AT)/(niPed - An))/6.0
    (PhiHat,dPhiHatdpsi) = generate_sin_profile(PhiPed,PhiPedGrad,psiMinPed,psiMaxPed)
    (etaiHat,detaiHatdpsi)=generate_eta_from_n_Phi_profile(niHat,PhiHat,TiHat,dniHatdpsi,dPhiHatdpsi,dTiHatdpsi,Zs[main_index],Delta,omega)


    THats = [TiHat]
    nHats = [niHat]
    etaHats = [etaiHat]
    dTHatdss = [dTiHatdpsi]
    dnHatdss = [dniHatdpsi]
    detaHatdss = [detaiHatdpsi]
    dPhiHatds = dPhiHatdpsi
    
    if simul.inputs.useGlobalTermMultiplier == 1:
        print "useGlobalTermMultiplier == 1 when using periodic boundary conditions. This is currently not supported and does not make sense. Will generate a trivial (constant 1) C profile."
    C=lambda x: 1.0 + 0*x

    if nonuniform == True:
        print "nonuniform grid is not supported with periodic boundary conditions yet: the default nonuniform grid would not make sense."
    (psi,dpsi_ds) = get_psiN(nonuniform=False,simul=simul,**kwargs)
    write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi,dpsi_ds,C)
 

    
# KEPT SINCE A PAIN TO CALCULATE
#nHats[imp_index] = etaHats[imp_index] *(nHats[main_index]/etaHats[main_index])**((Zs[imp_index]/Zs[main_index])*(THats[main_index]/THats[imp_index]))
#derivative calculate from sympy
#dnHatdpsis[imp_index]=(nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*((-Zs[imp_index]*THats[main_index]*dTHatdpsis[imp_index]/(Zs[main_index]*THats[imp_index]**2) + Zs[imp_index]*dTHatdpsis[main_index]/(Zs[main_index]*THats[imp_index]))*numpy.log(nHats[main_index]/etaHats[main_index]) + Zs[imp_index]*(dnHatdpsis[main_index]/etaHats[main_index] - nHats[main_index]*detaHatdpsis[main_index]/etaHats[main_index]**2)*THats[main_index]*etaHats[main_index]/(Zs[main_index]*THats[imp_index]*nHats[main_index]))*etaHats[imp_index] + (nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*detaHatdpsis[imp_index]
        
   


