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
    species=simul.species
    Nspecies=3
    main_index=kwargs["mI"]
    imp_index=kwargs["zI"]
    e_index=kwargs["eI"]

    Zs=simul.inputs.charges
    Z_1=Zs[main_index]

    #width of simulation region expressed in terms of pedestal width
    if "psiDiamFact" in kwargs.keys():
        psiDiamFact=kwargs["psiDiamFact"]
    else:
        psiDiamFact=5
    
    #how sharp bezier transition between region is compared to pedestal width
    if "transitionFact" in kwargs.keys():
        offset_frac=kwargs["transitionFact"]
    else:
        offset_frac=0.2

    
    Npsi=simul.inputs.Npsi
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax

    if "dxdpsiN" in kwargs.keys():
        print "Intepreting gradients as being given in x space"
        print "Will use function specified by dx_dpsi to convert" 
        dxdpsiN_at_a=kwargs["dxdpsiN"]
    else:
        print "Intepreting gradients as being given in psi_N."
        dxdpsiN_at_a=1 #trivial x=psi_N case

        
    x_width=kwargs["xwidth"]
    psiN_width=x_width/dxdpsiN_at_a

    #2015-12-16: -psiN_width/3.0 was the old default.
    if "upShift_denom" not in kwargs.keys():
        upShift=-psiN_width/3.0
    else:
        upShift=-psiN_width/kwargs["upShift_denom"]

    if "sameflux" in kwargs.keys():
        sameflux=kwargs["sameflux"]
    else:
        sameflux=False
    
    
    #calculate new psiMid and diameter for this profile
    if "midShift" not in kwargs.keys():
        midShift=0
    else:
        midShift=kwargs["midShift"]*psiN_width
    
    psiMid=1-psiN_width/2.0+midShift
    psiDiameter=psiDiamFact*psiN_width

    #already here, the simulation object changes
    simul.inputs.changevar("resolutionParameters","psiDiameter",psiDiameter)
    simul.inputs.changevar("physicsParameters","psiMid",psiMid)
    #this modifies the charge and mass of the species in the simulation
    species_filename= os.path.join(os.path.dirname(__file__), 'species_database.namelist')
    set_species_param(simul.species,species_filename,simul.norm,simul)
    
    simul.inputs.read(simul.input_filename)

    #get new psi coordinate from input
    psi=simul.inputs.psi
    
    
    psiMinPedIndex=numpy.argmin(numpy.abs(psi-(psiMid-psiN_width/2.0)))
    psiMaxPedIndex=numpy.argmin(numpy.abs(psi-(psiMid+psiN_width/2.0+upShift)))
    #psiMinPed=psi[psiMinPedIndex]
    #psiMaxPed=psi[psiMaxPedIndex]
    psiMinPed=psiMid-psiN_width/2.0
    psiMaxPed=psiMid+psiN_width/2.0+upShift
    
    
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
    T=THats[main_index]*simul.TBar


    
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

    # generate ion etas
    c=float(kwargs["imp_conc"])
    eta_1=neHatPre(psi)/(Z_1*(1+2*c))
    eta_2=c*eta_1

    etaHats[main_index] = eta_1
    detaHatdpsis[main_index] = simul.inputs.ddpsi_accurate(etaHats[main_index])

    etaHats[imp_index] = eta_2
    detaHatdpsis[imp_index] = simul.inputs.ddpsi_accurate(etaHats[imp_index])
      
    # generate Phi
    n_e=nHats[e_index] 
    ePhi=-T/(Z_1)*log(sqrt(1/(16*c**2) + n_e/(2*c*Z_1*eta_1))-1/(4*c))
    PhiHat=ePhi/simul.ePhiBar
    dPhiHatdpsi = simul.inputs.ddpsi_accurate(PhiHat)

    # generate ion densities
    n_1=eta_1*exp(-Z_1*ePhi/T)
    n_2=eta_2*exp(-2*Z_1*ePhi/T)
    
    nHats[main_index] = n_1
    dnHatdpsis[main_index] = simul.inputs.ddpsi_accurate(nHats[main_index])
    
    nHats[imp_index] = n_2
    dnHatdpsis[imp_index] = simul.inputs.ddpsi_accurate(nHats[imp_index])

    T2=simul.TBar*THats
    n2=simul.nBar*nHats
    if (sameflux==True):
        nt=nHats[main_index][0]
        nb=nHats[main_index][-1]
        Tp=Tpeds[main_index]
        breakpoint=psiMinPed
        dTdpsi=TCoreGrads[main_index]
        Delta0=psiN_width #pedestal width proper
        
        print "nt: " + str(nt)
        print "nb: " + str(nb)
        print "Tp: " + str(Tp)
        print "dTdpsi: " + str(dTdpsi)
        print "Delta0: " + str(Delta0)
        
        f=lambda x : x*(Tp-2*Delta0*x)**(3.0/2.0) - (nb*1.0/nt)*(Tp + 3*dTdpsi*Delta0)**(3.0/2.0)*dTdpsi
        dTdpsiTop=scipy.optimize.fsolve(f,0)
        dTdpsiBot=dTdpsi
        
        THatPre[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiTop*(psiN-breakpoint)))
        THatPed[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
        THatAft[main_index] =(lambda psiN: (Tpeds[main_index] + dTdpsiBot*(psiN-breakpoint)))
        dTHatPredpsi[main_index] =(lambda psiN: dTdpsiTop)
        dTHatPeddpsi[main_index] =(lambda psiN: dTdpsiBot)
        dTHatAftdpsi[main_index] =(lambda psiN: dTdpsiBot)
        
        Tlist=[THatPre[main_index],THatPed[main_index]]
        dTdpsilist=[dTHatPredpsi[main_index],dTHatPeddpsi[main_index]]

        THats[main_index]=bezier_transition(Tlist,[breakpoint],pairList[:-1],psi)
        dTHatdpsis[main_index]=simul.inputs.ddpsi_accurate(THats[main_index])

        #T2=simul.TBar*bezier_transition(Tlist,[breakpoint],pairList[:-1],numpy.array([(psiMinPed + psiMaxPed)/2.0]))[0]
        #n2=simul.nBar*bezier_transition(nilist,psiList,pairList,numpy.array([(psiMinPed + psiMaxPed)/2.0]))[0]
        THats[imp_index]=THats[main_index]
        dTHatdpsis[imp_index]=dTHatdpsis[main_index]

        
    # generate electron eta
    etaHats[e_index] = n_e*exp(-Zs[e_index]*PhiHat/THats[e_index])
    detaHatdpsis[e_index] = simul.inputs.ddpsi_accurate(etaHats[e_index])

    #with n_i and T_i generated, we can evaluate logLambda at a suitable point
    point=numpy.floor((psiMinPedIndex+psiMaxPedIndex)/2.0) # middle of the pedestal
    T=simul.TBar*THats[main_index][point]
    n=simul.nBar*nHats[main_index][point]
    print "T: "+str(T)
    print "n: "+str(n)
    T2=T2[0][point]
    n2=n2[0][point]
    print "T2: "+str(T2)
    print "n2: "+str(n2)

    
    print "Delta before:" + str(simul.Delta)
    if simul.units=="SI":
        logLambda=coulombLog(n2,T2)
        #print "logLambda: "+str(logLambda)
        nur=nu_r(simul.RBar,simul.nBar,simul.TBar,logLambda)
        print "nu_r: "+str(nur)
        simul.inputs.changevar("physicsParameters","nu_r",nur)
        simul.inputs.changevar("physicsParameters","Delta",simul.mBar*simul.vBar/(simul.eBar*simul.BBar*simul.RBar))
        simul.inputs.changevar("physicsParameters","omega",simul.ePhiBar/(simul.eBar*simul.BBar*simul.RBar*simul.vBar))
        simul.inputs.read(simul.input_filename)
    else:
        print "Only SI units are currently supported"
    print "Delta after:" + str(simul.Delta)

    
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
