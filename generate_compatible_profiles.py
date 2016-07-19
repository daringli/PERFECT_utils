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

def generate_compatible_profiles(simul,**kwargs):
    #simul should be a list of species names, simul.species
    #the position of a species in this list will determine the
    #index it gets in the arrays and will also be used to
    #intepret possible keywords in the keyword-argument dict kwargs
    species=simul.species
    Nspecies=len(species)

    if "allflat" in kwargs.keys():
        allflat=kwargs["allflat"]
    else:
        allflat=False

    if "nonuniform" in kwargs.keys():
        nonuniform=kwargs["nonuniform"]
    else:
        nonuniform=False

        
    if "sameflux" in kwargs.keys():
        sameflux=kwargs["sameflux"]
    else:
        sameflux=False

    if "sameeta" in kwargs.keys():
        sameeta=kwargs["sameeta"]
    else:
        sameeta=False

        
    if "samefluxshift" in kwargs.keys():
        samefluxshift=kwargs["samefluxshift"]
    else:
        samefluxshift=False

    if "oneSpecies" in kwargs.keys():
        oneSpecies=kwargs["oneSpecies"]
    else:
        oneSpecies=False
        
    if "twoSpecies" in kwargs.keys():
        twoSpecies=kwargs["twoSpecies"]
    else:
        twoSpecies=False

    if oneSpecies == True:
        main_index=kwargs["mI"]
    elif twoSpecies == True:
        main_index=kwargs["mI"]
        e_index=kwargs["eI"]
    else:
        main_index=kwargs["mI"]
        imp_index=kwargs["zI"]
        e_index=kwargs["eI"]
        
    if "zeroPhi" in kwargs.keys():
        zeroPhi=kwargs["zeroPhi"]
    else:
        zeroPhi=False

    if "specialEta" in kwargs.keys():
        specialEta=kwargs["specialEta"]
    else:
        specialEta=False

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



    #parse keywords to see how this object should be initialized
    if "dxdpsiN" in kwargs.keys():
        print "Intepreting gradients as being given in x space"
        print "Will use function specified by dx_dpsi to convert" 
        dxdpsiN_at_a=kwargs["dxdpsiN"]
    else:
        print "Intepreting gradients as being given in psi_N."
        dxdpsiN_at_a=1 #trivial x=psi_N case


    x_width=kwargs["xwidth"]
    psiN_width=x_width/dxdpsiN_at_a
    #calculate new psiMid and diameter for this profile
    if "midShift" not in kwargs.keys():
        midShift=0
    else:
        midShift=kwargs["midShift"]*psiN_width
    
    psiMid=1-psiN_width/2.0+midShift
    psiDiameter=psiDiamFact*psiN_width
    if allflat==True:
        leftBoundaryShift=-psiN_width
    else:
        leftBoundaryShift=0
    #print "psiDiameter"
    #print psiDiameter
    #print "prof gen psiDia:"
    #print type(psiDiameter)

    #already here, the simulation object changes
    simul.inputs.changevar("resolutionParameters","psiDiameter",psiDiameter)
    simul.inputs.changevar("physicsParameters","psiMid",psiMid)
    simul.inputs.changevar("physicsParameters","leftBoundaryShift",leftBoundaryShift)
    #this modifies the charge and mass of the species in the simulation
    species_filename= os.path.join(os.path.dirname(__file__), 'species_database.namelist')
    set_species_param(simul.species,species_filename,simul.norm,simul)
    
    simul.inputs.read(simul.input_filename)

    #get new psi coordinate from input
    psi=simul.inputs.psi
    if nonuniform:
        psiAHat= simul.psiAHat
        nonuniform_psi = simul.actual_psi
        d_nonuniform_psi_dpsiN = simul.psiAHatArray
    
    Npsi=simul.inputs.Npsi
    psiMin = simul.inputs.psiMin
    psiMax = simul.inputs.psiMax

    #might as well get a few other things from input now
    Zs=simul.inputs.charges
    if type(Zs) is not list:
        Zs=[Zs]
        # NO IDEA WHAT THIS IS FOR
        #if twoSpecies==True:
        #    Zs.append(0)
        #    imp_index=len(Zs)-1
    ms=simul.inputs.masses
    if not isinstance(ms, list):
        ms=[ms]

    #2015-12-16: -psiN_width/3.0 was the old default.
    if "upShift_denom" not in kwargs.keys():
        upShift=-psiN_width/3.0
    else:
        upShift=-psiN_width/kwargs["upShift_denom"]


    #determine gridpoints to start and stop pedestal
    psiMinPedIndex=numpy.argmin(numpy.abs(psi-(psiMid-psiN_width/2.0)))
    psiMaxPedIndex=numpy.argmin(numpy.abs(psi-(psiMid+psiN_width/2.0+upShift)))
    #psiMinPed=psi[psiMinPedIndex]
    #psiMaxPed=psi[psiMaxPedIndex]
    psiMinPed=psiMid-psiN_width/2.0
    psiMaxPed=psiMid+psiN_width/2.0+upShift
    if allflat==True:
        psiMinNotFlat=psiMid-2.5*psiN_width
        psiMaxNotFlat=psiMid+1.5*psiN_width
    #print "psiMinPed:"
    #print psiMinPed
    #list of psi where our profiles change slope
    psiList=[psiMinPed,psiMaxPed]
    print "pedestal start/stop: " + str(psiList)
    
    if allflat==True:
        psiList=[psiMinNotFlat,psiMinPed,psiMaxPed,psiMaxNotFlat]
    #distance in psi until our smooth transition between gradients is finished
    offset=(psiMaxPed-psiMinPed)*offset_frac
    pairList=[[offset,offset],[offset,offset]]
    if allflat==True:
        pairList=[[offset,offset],[offset,offset],[offset,offset],[offset,offset]]

    #allocate arrays
    THats = numpy.zeros((Nspecies,Npsi))
    dTHatdpsis = numpy.zeros((Nspecies,Npsi))
    nHats = numpy.zeros((Nspecies,Npsi))
    dnHatdpsis = numpy.zeros((Nspecies,Npsi))
    etaHats = numpy.zeros((Nspecies,Npsi))
    detaHatdpsis = numpy.zeros((Nspecies,Npsi))


    #get top of pedestal values
    #core gradients
    #SOL gradients
    #for each species from the keyword arguments
    #NOTE: uses the dx/dpsi at minor radius for both core and SOL, which assumes that
    #simulated region is small enough that it doesn't change much.
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
    #TpedGrads=Tpeds/psiN_width


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

    #generate T profiles
    for i in range(Nspecies):
        #Here we specify temperatures that should satisfy deltaT condition
        THatPre[i] =(lambda psiN: (Tpeds[i] + TCoreGrads[i]*(psiN-psiMinPed)))
        THatPed[i] =(lambda psiN: (Tpeds[i] + TpedGrads[i]*(psiN-psiMinPed)))
        THatAft[i] =(lambda psiN: (Tpeds[i] + TpedGrads[i]*(psiMaxPed-psiMinPed) + TSOLGrads[i]*(psiN-psiMaxPed)))
        dTHatPredpsi[i] = (lambda psiN: TCoreGrads[i])
        dTHatPeddpsi[i] = (lambda psiN: TpedGrads[i])
        dTHatAftdpsi[i] = (lambda psiN: TSOLGrads[i])
        Tlist=[THatPre[i],THatPed[i],THatAft[i]]
        dTdpsiList = [dTHatPredpsi[i],dTHatPeddpsi[i],dTHatAftdpsi[i]]
        if allflat==True:
            THatinner[i]=(lambda psiN: (Tpeds[i] + TCoreGrads[i]*(psiMinNotFlat-psiMinPed)) + TinnerGrad[i]*(psiN - psiMinNotFlat))
            THatouter[i]=(lambda psiN: (Tpeds[i] + TpedGrads[i]*(psiMaxPed-psiMinPed) + TSOLGrads[i]*(psiMaxNotFlat-psiMaxPed)) + TouterGrad[i]*(psiN - psiMaxNotFlat))
            dTHatinnerdpsi[i] = (lambda psiN: TouterGrad[i])
            dTHatouterdpsi[i]=(lambda psiN: TouterGrad[i])
            Tlist=[THatinner[i],THatPre[i],THatPed[i],THatAft[i],THatouter[i]]
            dTdpsiList = [dTHatinnerdpsi[i],dTHatPredpsi[i],dTHatPeddpsi[i],dTHatAftdpsi[i],dTHatouterdpsi[i]]

        if nonuniform:
            THats[i]=bezier_transition(Tlist,psiList,pairList,nonuniform_psi/psiAHat)
            dTHatdpsis[i]=derivative_bezier_transition(Tlist,dTdpsiList,psiList,pairList,nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
        else:
            THats[i]=bezier_transition(Tlist,psiList,pairList,psi)
            dTHatdpsis[i]=simul.inputs.ddpsi_accurate(THats[i])
        
    print "THat pedestal heights:" +str(Tpeds)
    print "THat inner boundary value:" +str(THats[:,0])
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
    if nonuniform:
        nHats[main_index]=bezier_transition(nilist,psiList,pairList,nonuniform_psi/psiAHat)
        dnHatdpsis[main_index]=derivative_bezier_transition(nilist,dnidpsiList,psiList,pairList,nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
    else:
        nHats[main_index]=bezier_transition(nilist,psiList,pairList,psi)
        dnHatdpsis[main_index] =simul.inputs.ddpsi_accurate(nHats[main_index])
    n2=simul.nBar*bezier_transition(nilist,psiList,pairList,numpy.array([(psiMinPed + psiMaxPed)/2.0]))[0]
    print "Ion nHat pedestal heights:" +str(niPed)
    print "Ion nHat inner boundary value:" +str(nHats[main_index][0])
   
    if (sameflux==True) or (samefluxshift==True):
        nt=nHats[main_index][0]
        nb=nHats[main_index][-1]
        Tp=Tpeds[main_index]
        breakpoint=psiMinPed
        if samefluxshift==True:
            breakpoint=psiMinPed-0.045
            i=main_index
            Tp=THatPre[i](breakpoint)
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

        if nonuniform:
            THats[main_index]=bezier_transition(Tlist,[breakpoint],pairList[:-1],nonuniform_psi/psiAHat)
            dTHatdpsis[main_index]=derivative_bezier_transition(Tlist,dTdpsilist,[breakpoint],pairList[:-1],nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
        else:
            THats[main_index]=bezier_transition(Tlist,[breakpoint],pairList[:-1],psi)
            dTHatdpsis[main_index]=simul.inputs.ddpsi_accurate(THats[main_index])
        T2=simul.TBar*bezier_transition(Tlist,[breakpoint],pairList[:-1],numpy.array([(psiMinPed + psiMaxPed)/2.0]))[0]
        if (twoSpecies==False and oneSpecies==False):
            THats[imp_index]=THats[main_index]
            dTHatdpsis[imp_index]=dTHatdpsis[main_index]

    #with n_i and T_i generated, we can evaluate logLambda at a suitable point
    point=numpy.floor((psiMinPedIndex+psiMaxPedIndex)/2.0) # middle of the pedestal
    T=simul.TBar*THats[main_index][point]
    n=simul.nBar*nHats[main_index][point]
    print "T: "+str(T)
    print "n: "+str(n)
    print "T2: "+str(T2)
    print "n2: "+str(n2)

    
    print "Delta before:" + str(simul.Delta)
    if simul.units=="SI":
        logLambda=coulombLog(n2,T2)
        print "logLambda: "+str(logLambda)
        nur=nu_r(simul.RBar,simul.nBar,simul.TBar,logLambda)
        simul.inputs.changevar("physicsParameters","nu_r",nur)
        simul.inputs.changevar("physicsParameters","Delta",simul.mBar*simul.vBar/(simul.eBar*simul.BBar*simul.RBar))
        simul.inputs.changevar("physicsParameters","omega",simul.ePhiBar/(simul.eBar*simul.BBar*simul.RBar*simul.vBar))
        simul.inputs.read(simul.input_filename)
    else:
        print "Only SI units are currently supported"
    print "Delta after:" + str(simul.Delta)

    etaiHatPre =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
    detaiHatPredpsi =(lambda psiN: niCoreGrad)
    
    
    if specialEta:
       etaiHatPed =(lambda psiN: (niPed+2*niCoreGrad* (psiN-psiMinPed)))
       detaiHatPeddpsi =(lambda psiN: 2*niCoreGrad)
    else:
        etaiHatPed =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
        detaiHatPeddpsi =(lambda psiN: niCoreGrad)

    i=main_index #to get the right lambda functions
    PhiTopPoint=psiMaxPed

    prefactor=1.0
    PhiTop=prefactor*(THatPed[main_index](PhiTopPoint)/Zs[main_index])*numpy.log(etaiHatPed(PhiTopPoint)*1.0/niHatPed(PhiTopPoint))*(2.0*simul.omega/simul.Delta)

    i=main_index #to get the right lambda functions
    etaiHatAft =(lambda psiN: niHatAft(psiN)*numpy.exp(PhiTop*Zs[main_index]/THatAft[main_index](psiN)))
    detaiHatAftdpsi = (lambda psiN: (dniHatAftdpsi(psiN) - niHatAft(psiN)*PhiTop*Zs[main_index]*dTHatAftdpsi[main_index](psiN)/(THatAft[main_index](psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/THatAft[main_index](psiN)))
    
    etailist=[etaiHatPre,etaiHatPed,etaiHatAft]
    detaidpsilist=[detaiHatPredpsi,detaiHatPeddpsi,detaiHatAftdpsi]
    
    if allflat==True:
        etaiHatinner =(lambda psiN: (niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMinNotFlat-psi[0]) + niinnerGrad*(psiN - psiMinNotFlat)))
        etaiHatouter = (lambda psiN: (niHatouter(psiN)*numpy.exp(PhiTop*Zs[main_index]/THatouter[main_index](psiN))))
        detaiHatinnerdpsi =(lambda psiN: niinnerGrad)
        
        detaiHatouterdpsi = (lambda psiN: (dniHatouterdpsi(psiN) - niHatouter(psiN)*PhiTop*Zs[main_index]*dTHatouterdpsi[main_index](psiN)/(THatouter[main_index](psiN))**2)*numpy.exp(PhiTop*Zs[main_index]/THatouter[main_index](psiN)))
        
        etailist=[etaiHatinner,etaiHatPre,etaiHatPed,etaiHatAft,etaiHatouter]
        detaidpsilist=[detaiHatinnerdpsi,detaiHatPredpsi,detaiHatPeddpsi,detaiHatAftdpsi,detaiHatouterdpsi]
    
    if nonuniform:
        etaHats[main_index]=bezier_transition(etailist,psiList,pairList,nonuniform_psi/psiAHat)
        detaHatdpsis[main_index] =derivative_bezier_transition(etailist,detaidpsilist,psiList,pairList,nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
    else:
        etaHats[main_index]=bezier_transition(etailist,psiList,pairList,psi)
        detaHatdpsis[main_index] =simul.inputs.ddpsi_accurate(etaHats[main_index])
    

    
    #if Phi=0 at top of the pedestal, this gives the top of the n_z pedestal.
    #To make n_z and n_i same at points, those points should satisfy
    # eta_z=n_i (n_i/eta_i)^(-[Zz/Zi] Ti/Tz)
    #etaHats[imp_index] = 0.01*nHats[main_index][psiMinPedIndex]
    if ((twoSpecies==False) and (oneSpecies==False)) and (sameeta==False):
        imp_conc=kwargs["imp_conc"]
        etaHats[imp_index]=imp_conc*(niPed+niCoreGrad*(nonuniform_psi/psiAHat-psiMinPed))
        detaHatdpsis[imp_index] = imp_conc*niCoreGrad*d_nonuniform_psi_dpsiN/psiAHat
        if allflat==True:
            etazHatInner=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMinNotFlat-psi[0]) + niinnerGrad*(psiN - psiMinNotFlat)))
            etazHatMiddle=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiN-psi[0])))
            etazHatOuter=(lambda psiN: imp_conc*(niPed-niCoreGrad*(psiMinPed-psi[0]) +  niCoreGrad* (psiMaxNotFlat-psi[0])))

            detazHatInnerdpsi=(lambda psiN: imp_conc*niinnerGrad)
            detazHatMiddledpsi=(lambda psiN: imp_conc*niCoreGrad)
            detazHatOuterdpsi=(lambda psiN: 0)
            
            
            etailist=[etazHatInner,etazHatMiddle,etazHatOuter]
            detaidpsilist=[detazHatInnerdpsi,detazHatMiddledpsi,detazHatOuterdpsi]

            if nonuniform:
                etaHats[imp_index]=bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],nonuniform_psi/psiAHat)
                detaHatdpsis[imp_index] = derivative_bezier_transition(etailist,detaidpsilist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
            else:
                etaHats[imp_index]=bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],psi)
                detaHatdpsis[imp_index] =simul.inputs.ddpsi_accurate(etaHats[imp_index])
        if specialEta==True:
            gradScale=4.0725
            etazHatInner=(lambda psiN: imp_conc*(niPed+niCoreGrad*(psiN-psiMinPed)))
            etazHatMiddle=(lambda psiN: imp_conc*(niPed-gradScale*niCoreGrad*(psiN-psiMinPed)))
            etazHatOuter=(lambda psiN: imp_conc*(niPed-gradScale*niCoreGrad*(psiMaxPed-psiMinPed) +  niCoreGrad* (psiN-psiMaxPed)))

            detazHatInnerdpsi=(lambda psiN: imp_conc*niCoreGrad)
            detazHatMiddledpsi=(lambda psiN: -imp_conc*gradScale*niCoreGrad)
            detazHatOuterdpsi=(lambda psiN: imp_conc*niCoreGrad)
            
            etailist=[etazHatInner,etazHatMiddle,etazHatOuter]
            detaidpsilist = [detazHatInnerdpsi,detazHatMiddledpsi,detazHatOuterdpsi]
            if nonuniform:
                etaHats[imp_index]=bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],nonuniform_psi/psiAHat)
                detaHatdpsis[imp_index] = derivative_bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],nonuniform_psi/psiAHat)*d_nonuniform_psi_dpsiN/psiAHat
            else:
                etaHats[imp_index]=bezier_transition(etailist,psiList[:1]+psiList[-1:],pairList[:1]+pairList[-1:],psi)
                detaHatdpsis[imp_index] =simul.inputs.ddpsi_accurate(etaHats[imp_index])

    if sameeta==True:
        etaHats[imp_index]=etaHats[main_index]
        detaHatdpsis[imp_index]=detaHatdpsis[main_index]
    #solve for Phi to make delta_etai the above value
    #delta_ni=delta_i_factor*dnHatdpsis[mI]/nHats[mI]
    #rhs=((psiAHat*B*numpy.sqrt(THats[mI]))/(Delta*I*numpy.sqrt(mi)))*(delta_etai-delta_ni)
    #A=D-numpy.diag(dlogTidpsi);
    #PhiHat=numpy.linalg.solve(A, rhs)
    PhiHat=numpy.log(etaHats[main_index]/nHats[main_index])*THats[main_index]/Zs[main_index]*simul.Delta/(2*simul.omega)
    dPhiHatdpsi=(-etaHats[main_index]*dnHatdpsis[main_index]/nHats[main_index]**2 + detaHatdpsis[main_index]/nHats[main_index])*THats[main_index]*nHats[main_index]/etaHats[main_index] + numpy.log(etaHats[main_index]/nHats[main_index])*dTHatdpsis[main_index]
    if zeroPhi:
        PhiHat=numpy.zeros(Npsi)
        dPhiHatdpsi=numpy.zeros(Npsi)

    if (twoSpecies==False) and (oneSpecies==False):
        nHats[imp_index] = etaHats[imp_index] *(nHats[main_index]/etaHats[main_index])**((Zs[imp_index]/Zs[main_index])*(THats[main_index]/THats[imp_index]))
        #derivative calculate from sympy
        dnHatdpsis[imp_index]=(nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*((-Zs[imp_index]*THats[main_index]*dTHatdpsis[imp_index]/(Zs[main_index]*THats[imp_index]**2) + Zs[imp_index]*dTHatdpsis[main_index]/(Zs[main_index]*THats[imp_index]))*numpy.log(nHats[main_index]/etaHats[main_index]) + Zs[imp_index]*(dnHatdpsis[main_index]/etaHats[main_index] - nHats[main_index]*detaHatdpsis[main_index]/etaHats[main_index]**2)*THats[main_index]*etaHats[main_index]/(Zs[main_index]*THats[imp_index]*nHats[main_index]))*etaHats[imp_index] + (nHats[main_index]/etaHats[main_index])**(Zs[imp_index]*THats[main_index]/(Zs[main_index]*THats[imp_index]))*detaHatdpsis[imp_index]
    
    if (twoSpecies==True) and (oneSpecies==False):
        nHats[e_index]=Zs[main_index]*nHats[main_index]
        dnHatdpsis[e_index]=Zs[main_index]*dnHatdpsis[main_index]
    elif (twoSpecies==False) and (oneSpecies==False):
        nHats[e_index]=Zs[imp_index]*nHats[imp_index] + Zs[main_index]*nHats[main_index]
        dnHatdpsis[e_index]=Zs[imp_index]*dnHatdpsis[imp_index]+ Zs[main_index]*dnHatdpsis[main_index]
    if (oneSpecies==False):
        etaHats[e_index]=nHats[e_index]*numpy.exp((Zs[e_index]*simul.omega*2/simul.Delta)*PhiHat/THats[e_index])
        detaHatdpsis[e_index]=(dnHatdpsis[e_index] + (2*simul.omega/simul.Delta)*(nHats[e_index]*Zs[e_index]/THats[e_index])*(dPhiHatdpsi-PhiHat*dTHatdpsis[e_index]/THats[e_index]))*numpy.exp((Zs[e_index]*simul.omega*2/simul.Delta)*PhiHat/THats[e_index])
    

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
        if twoSpecies==True:
            outputfile.create_profiles_for_Npsi(Npsi,2,PhiHat,dPhiHatdpsi,THats[:,[main_index,e_index]],dTHatdpsis[:,[main_index,e_index]],nHats[:,[main_index,e_index]],dnHatdpsis[:,[main_index,e_index]],etaHats[:,[main_index,e_index]],detaHatdpsis[:,[main_index,e_index]])
        elif oneSpecies==True:
            outputfile.create_profiles_for_Npsi(Npsi,1,PhiHat,dPhiHatdpsi,THats[:,[main_index]],dTHatdpsis[:,[main_index]],nHats[:,[main_index]],dnHatdpsis[:,[main_index]],etaHats[:,[main_index]],detaHatdpsis[:,[main_index]])
        else:
            outputfile.create_profiles_for_Npsi(Npsi,Nspecies,PhiHat,dPhiHatdpsi,THats,dTHatdpsis,nHats,dnHatdpsis,etaHats,detaHatdpsis)
