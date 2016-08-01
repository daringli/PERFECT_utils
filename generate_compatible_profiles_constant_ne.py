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

# Notation:
# 1,2: ion species 1 and 2.
# c: pseudo-concentration of species 2 compared to 1
#    c = \eta_2/\eta_1
# n_e: electron density
# Z: ion species 1 proton number
# T: ion temperature
# ASSUMES: T_1 = T_2. Z_2 = 2 Z_1

e=scipy.constants.e
log=numpy.log
exp=numpy.exp
sqrt=numpy.sqrt

def Phi_builder(T,c,Z,n_e):
    return lambda x : -(T(x)/e)*log(sqrt(1/(16*c(x)**2) + n_e(x)/(2*c(x)*Z))-1/(4*c(x)))

def n_builder(eta,Z,T,Phi):
    return lambda x : eta(x)*exp(-Z*e*Phi/T)

def const_builder(c):
    return lambda x: c
    

def generate_compatible_profiles_constant_ne(simul,**kwargs):
    Nspecies=3
    main_index=kwargs["mI"]
    imp_index=kwargs["zI"]
    e_index=kwargs["eI"]

    Zs=simul.inputs.charges
    Z_1=Zs[main_index]

    #how sharp bezier transition between region is compared to pedestal width
    if "transitionFact" in kwargs.keys():
        offset_frac=kwargs["transitionFact"]
    else:
        offset_frac=0.2

    Npsi=simul.inputs.Npsi
    psi=simul.inputs.psi
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax
    
    psiList=[psiMinPed,psiMaxPed]
    
    offset=(psiMaxPed-psiMinPed)*offset_frac
    pairList=[[offset,offset],[offset,offset]]
    
    #allocate arrays
    THats = numpy.zeros((Nspecies,Npsi))
    dTHatdpsis = numpy.zeros((Nspecies,Npsi))
    nHats = numpy.zeros((Nspecies,Npsi))
    dnHatdpsis = numpy.zeros((Nspecies,Npsi))
    etaHats = numpy.zeros((Nspecies,Npsi))
    detaHatdpsis = numpy.zeros((Nspecies,Npsi))

    #temperature generation
    TScale=[0]*Nspecies
    Tpeds=[0]*Nspecies
    TCoreGrads=[0]*Nspecies
    TpedGrads=[0]*Nspecies
    TSOLGrads=[0]*Nspecies

    THatPre=[0]*Nspecies
    THatPed=[0]*Nspecies
    THatAft=[0]*Nspecies

    dTHatPredpsi=[0]*Nspecies
    dTHatPeddpsi=[0]*Nspecies
    dTHatAftdpsi=[0]*Nspecies
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
    TScale=numpy.array(TScale)
    Tpeds=numpy.array(Tpeds)
    TCoreGrads=numpy.array(TCoreGrads)
    TSOLGrads=numpy.array(TSOLGrads)
    TpedGrads=numpy.array(TpedGrads)

    for i in range(Nspecies):
        #Here we finally generate temperatures that should satisfy deltaT condition
        THatPre[i] =(lambda psiN: (Tpeds[i] + TCoreGrads[i]*(psiN-psiMinPed)))
        THatPed[i] =(lambda psiN: (Tpeds[i] + TpedGrads[i]*(psiN-psiMinPed)))
        THatAft[i] =(lambda psiN: (Tpeds[i] + TpedGrads[i]*(psiMaxPed-psiMinPed) + TSOLGrads[i]*(psiN-psiMaxPed)))
        dTHatPredpsi[i] = (lambda psiN: TCoreGrads[i])
        dTHatPeddpsi[i] = (lambda psiN: TpedGrads[i])
        dTHatAftdpsi[i] = (lambda psiN: TSOLGrads[i])
        Tlist=[THatPre[i],THatPed[i],THatAft[i]]
        dTdpsiList = [dTHatPredpsi[i],dTHatPeddpsi[i],dTHatAftdpsi[i]]

        THats[i]=bezier_transition(Tlist,psiList,pairList,psi)
        dTHatdpsis[i]=simul.inputs.ddpsi_accurate(THats[i])
    T=T[main_index]


    
    #create n_e:
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
    
    dneHatPredpsi = (lambda psiN:  neCoreGrad)
    dneHatPeddpsi = (lambda psiN:  nepedGrad)
    dneHatAftdpsi = (lambda psiN:  neSOLGrad)
    
    nelist=[neHatPre,neHatPed,neHatAft]
    dnedpsiList = [dneHatPredpsi,dneHatPeddpsi,dneHatAftdpsi] 
    nHats[e_index] = bezier_transition(nelist,psiList,pairList,psi)
    dnHatdpsis[e_index] = simul.inputs.ddpsi_accurate(nHats[e_index])

    # generate Phi
    
    n_e=nHats[e_index]    
    eta1=kwargs["nped_"+species[main_index]]
    c=kwargs["imp_conc"]
    eta2=c*eta1
    #arrays
    eta1=const_builder(eta1)(psi)
    eta2=const_builder(eta2)(psi)
    c=const_builder(c)(psi)
    PhiHat=-T/(e*Z_1)*log(sqrt(1/(16*c**2) + n_e/(2*c*Z_1))-1/(4*c(x)))
    dPhiHatdpsi = simul.inputs.ddpsi_accurate(PhiHat)

    n_1=eta1*exp(-Z_1*Phi/T)
    n_2=eta2*exp(-2*Z_1*Phi/T)
    
    nHats[main_index] = n_1
    dnHatdpsis[main_index] = simul.inputs.ddpsi_accurate(nHats[main_index])

    nHats[imp_index] = n_2
    dnHatdpsis[imp_index] = simul.inputs.ddpsi_accurate(nHats[imp_index])


    #generate eta profiles
    etaHats[main_index] = eta1
    detaHatdpsis[main_index] = simul.inputs.ddpsi_accurate(etaHats[main_index])

    etaHats[imp_index] = eta2
    detaHatdpsis[imp_index] = simul.inputs.ddpsi_accurate(etaHats[imp_index])
    
    etaHats[e_index] = n_e*exp(-Zs[e_index]*Phi/THats[e_index])
    detaHatdpsis[e_index] = simul.inputs.ddpsi_accurate(etaHats[e_index])

    
    #tranpose
    #profiles got messed up otherwise
    THats=numpy.transpose(THats)
    dTHatdpsis=numpy.transpose(dTHatdpsis)
    nHats=numpy.transpose(nHats) 
    dnHatdpsis=numpy.transpose(dnHatdpsis)
    etaHats=numpy.transpose(etaHats) 
    detaHatdpsis=numpy.transpose(detaHatdpsis) 


    # Write profiles
    #print simul.input_dir+"/"+simul.inputs.profilesFilename
    try:
        outputfile = perfectProfiles(simul.input_dir+"/"+simul.inputs.profilesFilename)
    except IOError:
        print "#######################################################"
        print simul.input_dir+"/"+simul.inputs.profilesFilename + " already exists, cannot generate input profiles hdf5 file."
        print "#######################################################"
    else:
        outputfile.create_profiles_for_Npsi(Npsi,Nspecies,PhiHat,dPhiHatdpsi,THats,dTHatdpsis,nHats,dnHatdpsis,etaHats,detaHatdpsis)
