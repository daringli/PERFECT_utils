#this class generates n_i, n_z and Phi for a given n_e, T_e, T_i and T_z
#so that that the orderings in PERFECT are satisfied

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

from set_species_param import set_species_param


#print things
verbose=False

###############################
# GLOBAL VARIABLES
psiMinPed = None #pedestal start in psiN
psiMaxPed = None #pedestal end in psiN
psiMinNotFlat = None #for allflat
psiMaxNotFlat = None #for allflat
psi = None
dpsi_ds = None
Nspecies = None
allflat = None
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
niHatPed = None
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
    

def write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi):
    #tranpose
    THats = numpy.transpose(sample_funclist(THats,psi))
    dTHatdss = numpy.transpose(sample_funclist(dTHatdss,psi))
    nHats = numpy.transpose(sample_funclist(nHats,psi))
    dnHatdss = numpy.transpose(sample_funclist(dnHatdss,psi))
    etaHats = numpy.transpose(sample_funclist(etaHats,psi))
    detaHatdss = numpy.transpose(sample_funclist(detaHatdss,psi))
    PhiHat = sample_funclist(PhiHat,psi)
    dPhiHatds = sample_funclist(dPhiHatds,psi)

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
 

def get_psiN(nonuniform,simul):
    #get new psiN coordinate from input
    if nonuniform:
        actual_psi = simul.actual_psi
        dactual_psi_ds = simul.psiAHatArray
        #maps between psiN (here: psi) and psi (here: actual_psi)
        psiAHat= simul.psiAHat
        psi=actual_psi/psiAHat #nonuniform psiN
        dpsi_ds = dactual_psi_ds/psiAHat
    else:
        psi=simul.inputs.psi
        dpsi_ds = 1 #trivial mapping
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
    
    neHatPre =(lambda psiN: (nePed-neCoreGrad*(psiMinPed-psi[0]) +  neCoreGrad* (psiN-psi[0])))
    neHatPed =(lambda psiN: (nePed +  nepedGrad* (psiN-psiMinPed)))
    neHatAft =(lambda psiN: (nePed+nepedGrad*(psiMaxPed-psiMinPed) +  neSOLGrad* (psiN-psiMaxPed)))
    
    dneHatPredpsi = (lambda psiN:  neCoreGrad + psiN*0)
    dneHatPeddpsi = (lambda psiN:  nepedGrad + psiN*0)
    dneHatAftdpsi = (lambda psiN:  neSOLGrad + psiN*0)
    
    nelist=[neHatPre,neHatPed,neHatAft]
    dnedpsiList = [dneHatPredpsi,dneHatPeddpsi,dneHatAftdpsi] 
    (neHat,dneHatdpsi)=derivative_bezier_transition(nelist,dnedpsiList,psiList,pairList,psi)
    dneHatds = lambda x : dneHatdpsi(x)*dpsi_ds
    
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
    if allflat==True:
        niinnerGrad=niScale*0
        niouterGrad=niScale*0

    # generate n_i
    niHatPre =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
    niHatPed =(lambda psiN: (niPed +  nipedGrad* (psiN-psiMinPed)))
    niHatAft =(lambda psiN: (niPed+nipedGrad*(psiMaxPed-psiMinPed) +  niSOLGrad* (psiN-psiMaxPed)))
    
    dniHatPredpsi = (lambda psiN:  niCoreGrad)
    dniHatPeddpsi = (lambda psiN:  nipedGrad)
    dniHatAftdpsi = (lambda psiN:  niSOLGrad)
    
    nilist=[niHatPre,niHatPed,niHatAft]
    dnidpsiList = [dniHatPredpsi,dniHatPeddpsi,dniHatAftdpsi] 
    if allflat==True:
        niHatinner =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMinNotFlat-psi[0])) + niinnerGrad*(psiN - psiMinNotFlat))
        niHatouter =(lambda psiN: (niPed+nipedGrad*(psiMaxPed-psiMinPed) +  niSOLGrad* (psiMaxNotFlat-psiMaxPed)) + niouterGrad*(psiN - psiMaxNotFlat))
        dniHatinnerdpsi = (lambda psiN: niinnerGrad)
        dniHatouterdpsi = (lambda psiN: niouterGrad)
        nilist=[niHatinner,niHatPre,niHatPed,niHatAft,niHatouter]
        dnidpsiList = [dniHatinnerdpsi,dniHatPredpsi,dniHatPeddpsi,dniHatAftdpsi,dniHatouterdpsi]

    (niHat,dniHatdpsi)=derivative_bezier_transition(nilist,dnidpsiList,psiList,pairList,psi)
    dniHatds = lambda x : dniHatdpsi(x)*dpsi_ds

    #n2=simul.nBar*bezier_transition(nilist,psiList,pairList,numpy.array([psiMidPed]))[0]

    return (niHat,dniHatds)

def generate_T_profiles_sameflux(nt,nb,breakpoint_shift,**kwargs):
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
    if allflat==True:
        TinnerGrad=[0]*Nspecies
        TouterGrad=[0]*Nspecies
    #Pre,Ped,Aft will contain functions describe T in core, ped and outwards.
    THatPre=[0]*Nspecies
    THatPed=[0]*Nspecies
    THatAft=[0]*Nspecies

    dTHatPredpsi=[0]*Nspecies
    dTHatPeddpsi=[0]*Nspecies
    dTHatAftdpsi=[0]*Nspecies
    if allflat==True:
        THatinner=[0]*Nspecies
        THatouter=[0]*Nspecies
        dTHatinnerdpsi=[0]*Nspecies
        dTHatouterdpsi=[0]*Nspecies
    #reading T related parameters
    for i in range(Nspecies):
        if "TScale_"+species[i] in kwargs.keys():
            TScale[i]=kwargs["TScale_"+species[i]]
        else:
            TScale[i]=1.0
        Tpeds[i]=TScale[i]*kwargs["Tped_"+species[i]]
        TCoreGrads[i]=TScale[i]*kwargs["dTCoredx_"+species[i]]*dxdpsiN_at_a
        TpedGrads[i]=TScale[i]*kwargs["dTpeddx_"+species[i]]*dxdpsiN_at_a
        #print "-------------------------"
        #print "TpedGrads["+str(i)+"] = " + str(TpedGrads[i])
        #print "-------------------------"
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        if allflat==True:
            TinnerGrad=[0]*Nspecies
            TouterGrad=[0]*Nspecies
    TScale=numpy.array(TScale)
    Tpeds=numpy.array(Tpeds)
    TCoreGrads=numpy.array(TCoreGrads)
    TSOLGrads=numpy.array(TSOLGrads)
    TpedGrads=numpy.array(TpedGrads)
    if allflat==True:
        TinnerGrad=numpy.array(TinnerGrad)
        TouterGrad=numpy.array(TouterGrad)
    #############################################
    #generate T profiles
    #############################################
    
    breakpoint=psiMinPed + breakpoint_shift
    Tp=Tpeds[main_index] + TCoreGrads[main_index] *(breakpoint-psiMinPed)

    if (TCoreGrads[main_index] != TSOLGrads[main_index]) and (TPedGrads[main_index] != TSOLGrads[main_index]):
        print "generate_compatible_profiles: Warning: sameflux option uses the Core T gradients of the main species everywhere"
    dTdpsi=TCoreGrads[main_index]
    d1=breakpoint - psi[0]
    d2=psi[-1] - breakpoint
    
    f=lambda x : x*(Tp-d1*x)**(3.0/2.0) - (nb*1.0/nt)*(Tp + dTdpsi*d2)**(3.0/2.0)*dTdpsi
    dTdpsiTop=scipy.optimize.fsolve(f,0)[0] #[0] since returns an numpy.ndarray
    dTdpsiBot=dTdpsi

    THatPre[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiTop*(psiN-breakpoint)))
    THatPed[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
    THatAft[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
    dTHatPredpsi[main_index] =(lambda psiN: dTdpsiTop)
    dTHatPeddpsi[main_index] =(lambda psiN: dTdpsiBot)
    dTHatAftdpsi[main_index] =(lambda psiN: dTdpsiBot)

    
    Tlist=[THatPre[main_index],THatPed[main_index]]
    dTdpsilist=[dTHatPredpsi[main_index],dTHatPeddpsi[main_index]]

    (THats[main_index],dTHatdpsis[main_index])=derivative_bezier_transition(Tlist,dTdpsilist,[breakpoint],pairList[:-1],psi)    
    dTHatdss[main_index]= lambda x : dTHatdpsis[main_index](x)*dpsi_ds

    #T2=simul.TBar*bezier_transition(Tlist,[breakpoint],pairList[:-1],numpy.array([psiMidPed]))[0]


    if Nspecies >= 2:
        #T_e:
        THatPre[e_index] =(lambda psiN: (Tpeds[e_index] + TCoreGrads[e_index]*(psiN-psiMinPed)))
        THatPed[e_index] =(lambda psiN: (Tpeds[e_index] + TpedGrads[e_index]*(psiN-psiMinPed)))
        THatAft[e_index] =(lambda psiN: (Tpeds[e_index] + TpedGrads[e_index]*(psiMaxPed-psiMinPed) + TSOLGrads[e_index]*(psiN-psiMaxPed)))
        dTHatPredpsi[e_index] = (lambda psiN: TCoreGrads[e_index])
        dTHatPeddpsi[e_index] = (lambda psiN: TpedGrads[e_index])
        dTHatAftdpsi[e_index] = (lambda psiN: TSOLGrads[e_index])
        Tlist=[THatPre[e_index],THatPed[e_index],THatAft[e_index]]
        dTdpsiList = [dTHatPredpsi[e_index],dTHatPeddpsi[e_index],dTHatAftdpsi[e_index]]
        #PROBABLY NOT COMPATIBLE WITH SAMEFLUX
        #if allflat==True:
        #    THatinner[e_index]=(lambda psiN: (Tpeds[e_index] + TCoreGrads[e_index]*(psiMinNotFlat-psiMinPed)) + TinnerGrad[e_index]*(psiN - psiMinNotFlat))
        #    THatouter[e_index]=(lambda psiN: (Tpeds[e_index] + TpedGrads[e_index]*(psiMaxPed-psiMinPed) + TSOLGrads[e_index]*(psiMaxNotFlat-psiMaxPed)) + TouterGrad[e_index]*(psiN - psiMaxNotFlat))
        #    dTHatinnerdpsi[e_index] = (lambda psiN: TouterGrad[e_index])
        #    dTHatouterdpsi[e_index]=(lambda psiN: TouterGrad[e_index])
        #    Tlist=[THatinner[e_index],THatPre[e_index],THatPed[e_index],THatAft[e_index],THatouter[e_index]]
        #    dTdpsiList = [dTHatinnerdpsi[e_index],dTHatPredpsi[e_index],dTHatPeddpsi[e_index],dTHatAftdpsi[e_index],dTHatouterdpsi[e_index]]

        (THats[e_index],dTHatdpsis[e_index])=derivative_bezier_transition(Tlist,dTdpsiList,psiList,pairList,psi)
        dTHatdss[e_index]=lambda x : dTHatdpsis[e_index](x)*dpsi_ds
        
    if Nspecies == 3:
        THats[imp_index]=THats[main_index]
        dTHatdpsis[imp_index]=dTHatdpsis[main_index]
        dTHatdss[imp_index]=dTHatdss[main_index]


    global TiHatPed
    global TiHatAft
    global dTiHatAftdpsi
    TiHatPed = THatPed[main_index]
    TiHatAft = THatAft[main_index]
    dTiHatAftdpsi = dTHatAftdpsi[main_index]
    return (THats,dTHatdss)
    
def generate_T_profiles(**kwargs):
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
    if allflat==True:
        TinnerGrad=[0]*Nspecies
        TouterGrad=[0]*Nspecies
    #Pre,Ped,Aft will contain functions describe T in core, ped and outwards.
    THatPre=[0]*Nspecies
    THatPed=[0]*Nspecies
    THatAft=[0]*Nspecies

    dTHatPredpsi=[0]*Nspecies
    dTHatPeddpsi=[0]*Nspecies
    dTHatAftdpsi=[0]*Nspecies
    if allflat==True:
        THatinner=[0]*Nspecies
        THatouter=[0]*Nspecies
        dTHatinnerdpsi=[0]*Nspecies
        dTHatouterdpsi=[0]*Nspecies
    #reading T related parameters
    for i in range(Nspecies):
        if "TScale_"+species[i] in kwargs.keys():
            TScale[i]=kwargs["TScale_"+species[i]]
        else:
            TScale[i]=1.0
        Tpeds[i]=TScale[i]*kwargs["Tped_"+species[i]]
        TCoreGrads[i]=TScale[i]*kwargs["dTCoredx_"+species[i]]*dxdpsiN_at_a
        TpedGrads[i]=TScale[i]*kwargs["dTpeddx_"+species[i]]*dxdpsiN_at_a
        #print "-------------------------"
        #print "TpedGrads["+str(i)+"] = " + str(TpedGrads[i])
        #print "-------------------------"
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        TSOLGrads[i]=TScale[i]*kwargs["dTSOLdx_"+species[i]]*dxdpsiN_at_a
        if allflat==True:
            TinnerGrad=[0]*Nspecies
            TouterGrad=[0]*Nspecies
    TScale=numpy.array(TScale)
    Tpeds=numpy.array(Tpeds)
    TCoreGrads=numpy.array(TCoreGrads)
    TSOLGrads=numpy.array(TSOLGrads)
    TpedGrads=numpy.array(TpedGrads)
    if allflat==True:
        TinnerGrad=numpy.array(TinnerGrad)
        TouterGrad=numpy.array(TouterGrad)


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
        if allflat==True:
            THatinner[i]=(lambda psiN,i=i: (Tpeds[i] + TCoreGrads[i]*(psiMinNotFlat-psiMinPed)) + TinnerGrad[i]*(psiN - psiMinNotFlat))
            THatouter[i]=(lambda psiN,i=i: (Tpeds[i] + TpedGrads[i]*(psiMaxPed-psiMinPed) + TSOLGrads[i]*(psiMaxNotFlat-psiMaxPed)) + TouterGrad[i]*(psiN - psiMaxNotFlat))
            dTHatinnerdpsi[i] = (lambda psiN,i=i: TouterGrad[i])
            dTHatouterdpsi[i]=(lambda psiN,i=i: TouterGrad[i])
            Tlist=[THatinner[i],THatPre[i],THatPed[i],THatAft[i],THatouter[i]]
            dTdpsiList = [dTHatinnerdpsi[i],dTHatPredpsi[i],dTHatPeddpsi[i],dTHatAftdpsi[i],dTHatouterdpsi[i]]

        (THats[i],dTHatdpsis[i])=derivative_bezier_transition(Tlist,dTdpsiList,psiList,pairList,psi)
        dTHatdss[i]=lambda x,i=i: dTHatdpsis[i](x)*dpsi_ds
        
    global TiHatPed
    global TiHatAft
    global dTiHatAftdpsi
    TiHatPed = THatPed[main_index]
    TiHatAft = THatAft[main_index]
    dTiHatAftdpsi = dTHatAftdpsi[main_index]
    
    return (THats,dTHatdss)

def generate_etai_profile(Delta,omega,**kwargs):
    etaiHatPre =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
    detaiHatPredpsi =(lambda psiN: niCoreGrad)
    
    
    if specialeta:
       etaiHatPed =(lambda psiN: (niPed+2*niCoreGrad* (psiN-psiMinPed)))
       detaiHatPeddpsi =(lambda psiN: 2*niCoreGrad)
    else:
        etaiHatPed =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
        detaiHatPeddpsi =(lambda psiN: niCoreGrad)

    i=main_index #to get the right lambda functions
    PhiTopPoint=psiMaxPed

    prefactor=1.0
    PhiTop=prefactor*(TiHatPed(PhiTopPoint)/Zs[main_index])*numpy.log(etaiHatPed(PhiTopPoint)*1.0/niHatPed(PhiTopPoint))*(2.0*omega/Delta)

    i=main_index #to get the right lambda functions
    etaiHatAft =(lambda psiN: niHatAft(psiN)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))
    detaiHatAftdpsi = (lambda psiN: (dniHatAftdpsi(psiN) - niHatAft(psiN)*PhiTop*Zs[main_index]*dTiHatAftdpsi(psiN)/(TiHatAft(psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/TiHatAft(psiN)))
    
    etailist=[etaiHatPre,etaiHatPed,etaiHatAft]
    detaidpsilist=[detaiHatPredpsi,detaiHatPeddpsi,detaiHatAftdpsi]
    
    if allflat==True:
        etaiHatinner =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMinNotFlat-psi[0]) + niinnerGrad*(psiN - psiMinNotFlat)))
        etaiHatouter = (lambda psiN: (niHatouter(psiN)*numpy.exp(PhiTop*Zs[main_index]/THatouter[main_index](psiN))))
        detaiHatinnerdpsi =(lambda psiN: niinnerGrad)
        
        detaiHatouterdpsi = (lambda psiN: (dniHatouterdpsi(psiN) - niHatouter(psiN)*PhiTop*Zs[main_index]*dTHatouterdpsi[main_index](psiN)/(THatouter[main_index](psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/THatouter[main_index](psiN)))
        
        etailist=[etaiHatinner,etaiHatPre,etaiHatPed,etaiHatAft,etaiHatouter]
        detaidpsilist=[detaiHatinnerdpsi,detaiHatPredpsi,detaiHatPeddpsi,detaiHatAftdpsi,detaiHatouterdpsi]
    
    (etaiHat,detaiHatdpsi) =derivative_bezier_transition(etailist,detaidpsilist,psiList,pairList,psi)
    detaiHatds= lambda x: detaiHatdpsi(x)*dpsi_ds

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

        (etazHat,detazHatdpsi) = derivative_bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],psi)
        detazHatds = lambda x : detazHatdpsi(x) * dpsi_ds
        
    elif allflat==True:
        etazHatInner=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMinNotFlat-psi[0]) + niinnerGrad*(psiN - psiMinNotFlat)))
        etazHatMiddle=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
        etazHatOuter=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMaxNotFlat-psi[0])))

        detazHatInnerdpsi=(lambda psiN: imp_conc*niinnerGrad)
        detazHatMiddledpsi=(lambda psiN: imp_conc*niCoreGrad)
        detazHatOuterdpsi=(lambda psiN: 0)


        etailist=[etazHatInner,etazHatMiddle,etazHatOuter]
        detaidpsilist=[detazHatInnerdpsi,detazHatMiddledpsi,detazHatOuterdpsi]

        (etazHat,detazHatdpsi) = derivative_bezier_transition(etailist,detaidpsilist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],psi)
        detazHatds = lambda x : detazHatdpsi(x) * dpsi_ds
        
    else:
        etazHat=lambda x : imp_conc*(niPed+niCoreGrad*(x-psiMinPed))
        detazHatds = lambda x: imp_conc*niCoreGrad*dpsi_ds + 0*x #add 0*x to make it return ndarray when x is one.
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

def generate_compatible_profiles(simul,xwidth,is_allflat=False,nonuniform=False,sameflux=False,sameeta=False,samefluxshift=0,zeroPhi=False,specialEta=False,psiDiamFact=5,transitionFact=0.2,dxdpsiN=1,midShift=0,upShift_denom=3.0,**kwargs):
    #NOTE: uses the dx/dpsi at minor radius for both core and SOL, which assumes that
    #simulated region is small enough that it doesn't change much.

    global psiMinPed
    global psiMaxPed
    global psi
    global dpsi_ds
    global Nspecies
    global species
    global allflat
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
    allflat=is_allflat
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

    if allflat==True:
        leftBoundaryShift=-psiN_ped_width
    else:
        leftBoundaryShift=0
        
    update_domain_size(simul,psiN_ped_width,midShift,psiDiamFact,leftBoundaryShift,rightBoundaryShift=0)
    (Delta,omega) = update_Delta_omega(simul)
    psiMid = simul.inputs.psiMid
    
    #get new psiN coordinate from input
    (psi,dpsi_ds) = get_psiN(nonuniform,simul)
    
    Npsi=simul.inputs.Npsi
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax

    #determine gridpoints to start and stop pedestal
    psiMinPed=psiMid-psiN_ped_width/2.0
    psiMaxPed=psiMid+psiN_ped_width/2.0+upShift
    psiMidPed=(psiMinPed+psiMaxPed)/2.0
    
    if allflat==True:
        global psiMinNotFlat
        global psiMaxNotFlat
        psiMinNotFlat=psiMid-2.5*psiN_ped_width
        psiMaxNotFlat=psiMid+1.5*psiN_ped_width
    #list of psi where our profiles change slope
    psiList=[psiMinPed,psiMaxPed]
    print "pedestal start,stop: " + str(psiList)
    
    if allflat==True:
        psiList=[psiMinNotFlat,psiMinPed,psiMaxPed,psiMaxNotFlat]
    #distance in psi until our smooth transition between gradients is finished
    offset=(psiMaxPed-psiMinPed)*transitionFact
    pairList=[[offset,offset],[offset,offset]]
    if allflat==True:
        pairList=[[offset,offset],[offset,offset],[offset,offset],[offset,offset]]

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

    #generate n_i profile
    (nHats[main_index],dnHatdss[main_index]) = generate_ni_profile(**kwargs)
       
    if (sameflux==True):
        nt=nHats[main_index](psi[0])
        nb=nHats[main_index](psi[-1])
        (THats,dTHatdss) = generate_T_profiles_sameflux(nt,nb,samefluxshift,**kwargs)
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

    (etaHats[main_index],detaHatdss[main_index]) = generate_etai_profile(Delta,omega,**kwargs)

    if zeroPhi:
        PhiHat= lambda x : 0*x
        dPhiHatds= lambda x : 0*x
    else:
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

        # KEPT SINCE A PAIN TO CALCULATE
        #nHats[imp_index] = etaHats[imp_index] *(nHats[main_index]/etaHats[main_index])**((Zs[imp_index]/Zs[main_index])*(THats[main_index]/THats[imp_index]))
        #derivative calculate from sympy
        #dnHatdpsis[imp_index]=(nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*((-Zs[imp_index]*THats[main_index]*dTHatdpsis[imp_index]/(Zs[main_index]*THats[imp_index]**2) + Zs[imp_index]*dTHatdpsis[main_index]/(Zs[main_index]*THats[imp_index]))*numpy.log(nHats[main_index]/etaHats[main_index]) + Zs[imp_index]*(dnHatdpsis[main_index]/etaHats[main_index] - nHats[main_index]*detaHatdpsis[main_index]/etaHats[main_index]**2)*THats[main_index]*etaHats[main_index]/(Zs[main_index]*THats[imp_index]*nHats[main_index]))*etaHats[imp_index] + (nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*detaHatdpsis[imp_index]
    
    if (Nspecies > 1):
        #n_e from quasi-neutrality
        if (Nspecies == 2):
            nHats[e_index]=lambda x : Zs[main_index]*nHats[main_index](x)
            dnHatdss[e_index]= lambda x : Zs[main_index]*dnHatdss[main_index](x)
        elif (Nspecies == 3):
            nHats[e_index]=lambda x : Zs[imp_index]*nHats[imp_index](x) + Zs[main_index]*nHats[main_index](x)
            dnHatdss[e_index]=lambda x: Zs[imp_index]*dnHatdss[imp_index](x) + Zs[main_index]*dnHatdss[main_index](x)
        #eta_e from n_e
        (etaHats[e_index], detaHatdss[e_index]) = generate_eta_from_n_Phi_profile(nHats[e_index],PhiHat,THats[e_index],dnHatdss[e_index],dPhiHatds,dTHatdss[e_index],Zs[e_index],Delta,omega)
        
    write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi)

##########################################################
# generate_compatible_profiles_constant_ne
#########################################################

def generate_compatible_profiles_constant_ne(simul,xwidth,is_allflat=False,nonuniform=False,sameflux=False,oldsameflux=False,sameeta=False,samefluxshift=0,zeroPhi=False,specialEta=False,psiDiamFact=5,transitionFact=0.2,dxdpsiN=1,midShift=0,upShift_denom=3.0,**kwargs):
    e=scipy.constants.e
    log=numpy.log
    exp=numpy.exp
    sqrt=numpy.sqrt

    global psiMinPed
    global psiMaxPed
    global psi
    global dpsi_ds
    global Nspecies
    global species
    global allflat
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
    allflat=is_allflat
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

    if allflat==True:
        leftBoundaryShift=-psiN_ped_width
    else:
        leftBoundaryShift=0
        
    update_domain_size(simul,psiN_ped_width,midShift,psiDiamFact,leftBoundaryShift,rightBoundaryShift=0)
    (Delta,omega) = update_Delta_omega(simul)
    psiMid = simul.inputs.psiMid
    
    #get new psiN coordinate from input
    (psi,dpsi_ds) = get_psiN(nonuniform,simul)
    
    Npsi=simul.inputs.Npsi
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax

    #determine gridpoints to start and stop pedestal
    psiMinPed=psiMid-psiN_ped_width/2.0
    psiMaxPed=psiMid+psiN_ped_width/2.0+upShift
    psiMidPed=(psiMinPed+psiMaxPed)/2.0
    
    if allflat==True:
        global psiMinNotFlat
        global psiMaxNotFlat
        psiMinNotFlat=psiMid-2.5*psiN_ped_width
        psiMaxNotFlat=psiMid+1.5*psiN_ped_width
    #list of psi where our profiles change slope
    psiList=[psiMinPed,psiMaxPed]
    print "pedestal start,stop: " + str(psiList)
    
    if allflat==True:
        psiList=[psiMinNotFlat,psiMinPed,psiMaxPed,psiMaxNotFlat]
    #distance in psi until our smooth transition between gradients is finished
    offset=(psiMaxPed-psiMinPed)*transitionFact
    pairList=[[offset,offset],[offset,offset]]
    if allflat==True:
        pairList=[[offset,offset],[offset,offset],[offset,offset],[offset,offset]]

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
    Z_1=Zs[main_index]
    Z_2=Zs[imp_index]

    #create n_e:
    (nHats[e_index],dnHatdss[e_index]) = generate_ne_profile(simul,**kwargs)

    # generate ion etas
    c=float(kwargs["imp_conc"])
    etaHats[main_index] = lambda x: neHatPre(x)/(Z_1*(1+2*c))
    detaHatdss[main_index] =  lambda x: dpsi_ds*dneHatPredpsi(x)/(Z_1*(1+2*c))
    etaHats[imp_index] = lambda x: c*etaHats[main_index](x)
    detaHatdss[imp_index] = lambda x: c*detaHatdss[main_index](x)

    # generate ion densities
    X = lambda x: sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Z_1*etaHats[main_index](x)))-1/(4*c)
    dXds = lambda x: 1/(sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Z_1*etaHats[main_index](x)))*4*Z_1*c*etaHats[main_index](x))*(dnHatdss[e_index](x) - nHats[e_index](x)*(detaHatdss[main_index](x)/etaHats[main_index](x)))

    (nHats[main_index],dnHatdss[main_index]) = generate_n_from_eta_X_profile(etaHats[main_index],X,detaHatdss[main_index],dXds,Zs[main_index])

    (nHats[imp_index],dnHatdss[imp_index]) = generate_n_from_eta_X_profile(etaHats[imp_index],X,detaHatdss[imp_index],dXds,Zs[imp_index])
    
    

    #generate T
    if (sameflux==True):
        nt=nHats[main_index](psi[0])
        nb=nHats[main_index](psi[-1])
        (THats,dTHatdss) = generate_T_profiles_sameflux(nt,nb,samefluxshift,**kwargs)
    elif (oldsameflux==True):
        #this setting uses the old ni profile generation to generatea  temporary
        # n_i profile to generate sameflux T profiles that are the same as in
        # the constant_Phi case.
        # Useful for comparing const_Phi and const_ne with the same T
        #NOTE: need inputs specifying n_i to generate the temporary profile
        # does not really match the heat flux proxy for the actual n_i profiles
        # used in the const n_e case.
        (nH,dnHds) = generate_ni_profile(**kwargs)
        nt=nH(psi[0])
        nb=nH(psi[-1])
        (THats,dTHatdss) = generate_T_profiles_sameflux(nt,nb,samefluxshift,**kwargs)
    else:
        (THats,dTHatdss) = generate_T_profiles(**kwargs)
        
   

    # generate Phi
    PhiHat = lambda x: -(Delta/(Z_1*2*omega))*THats[main_index](x)*log(sqrt(1/(16*c**2) + nHats[e_index](x)/(2*c*Z_1*etaHats[main_index](x)))-1/(4*c))
    dPhiHatds = lambda x: -(Delta/(Z_1*2*omega))*(dTHatdss[main_index](x)*log(X(x)) + THats[main_index](x)*dXds(x)/X(x))
    
    T=simul.TBar*THats[main_index](psiMidPed)
    n=simul.nBar*nHats[main_index](psiMidPed)
    print "T: "+str(T)
    print "n: "+str(n)
    nu_r=update_nu_r(simul,n,T)
    print "nu_r: " + str(nu_r)

        
    # generate electron eta
    (etaHats[e_index],detaHatdss[e_index]) = generate_eta_from_n_Phi_profile(nHats[e_index],PhiHat,THats[e_index],dnHatdss[e_index],dPhiHatds,dTHatdss[e_index],Zs[e_index],Delta,omega)

    write_outputs(simul,THats,dTHatdss,nHats,dnHatdss,etaHats,detaHatdss,PhiHat,dPhiHatds,psi)

