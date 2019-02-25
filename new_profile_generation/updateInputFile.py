import os, numpy
from setSpeciesParameters import speciesParams
from coulomb_logarithms import logLambda_ee, logLambda_ie
from nu_r import nu_r #nu_r
from input_psiN_params_from_geometry_file import input_psiN_params_from_geometry_file

def syncSpeciesParams(simulation,speciesFilename="species_database.namelist"):
    normFilename = simulation.norm
    speciesList = simulation.species_list
    speciesFilename= os.path.join(os.path.dirname(__file__), speciesFilename)
    [Z,mHat]=speciesParams(speciesList,speciesFilename,normFilename)
    simulation.inputs.charges=Z
    simulation.inputs.masses=mHat
    if not isinstance(Z, (list, tuple, numpy.ndarray)):
        Z=numpy.array([Z])
    if not isinstance(mHat, (list, tuple, numpy.ndarray)):
        mHat=numpy.array([mHat])
    return (Z,mHat)

def updateDeltaOmega(simul):
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

def updatePsi(simul,geometryFilename="input.geometry.h5"):
    #print "Delta before:" + str(simul.Delta)
    if simul.inputs.geometryToUse == 4:
    
        if simul.units=="SI":
            (psiMid,psiDiameter,psiA) = input_psiN_params_from_geometry_file(geometryFilename)
                
            assert (simul.RBar == 1)
            assert (simul.BBar == 1)
            if psiA is not None:
                simul.inputs.psiAHat=psiA
            simul.inputs.psiMid=psiMid
            simul.inputs.psiDiameter=psiDiameter
        else:
            print "Only SI units are currently supported when reading from input.geometry.h5"
        return (psiMid,psiDiameter,psiA)
    else:
        # nothing to do since we do not need to sync from geometry file.
        pass 
    #print "Delta after:" + str(Delta)

def updateNuR(simul,n,T,logLambda="ee"):
    if simul.units=="SI":
        if logLambda == "ee":
            logLambda=logLambda_ee(n,T)
        elif (logLambda == "ie") or logLambda == "ei":
            logLambda=logLambda_ie(n,T)
        else:
            raise ValueError("Unrecognized Coulomb logarithm calculation: '" + logLambda +"'")
        #n,T not used beyond this point
        print "logLambda: "+str(logLambda)
        print "logLambda_ee: "+str(logLambda_ee(n,T))
        print "logLambda_ie: "+str(logLambda_ie(n,T))
        
        nur = nu_r(simul.RBar,simul.nBar,simul.TBar,logLambda)
        simul.inputs.nu_r=nur
        return nur
    else:
        print "Only SI units are currently supported"
