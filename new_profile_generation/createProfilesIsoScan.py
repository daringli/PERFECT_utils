from generateConstPhiIsoScan import generateConstPhi
from updateInputFile import syncSpeciesParams, updateDeltaOmega, updateNuR
from findInProfileList import extractProfilesFromList
from writeProfileListToOutput import writeProfileListToOutput
from setPedestalParams import setPedestalStartStop
from profileGrid import ProfileGrid

def createProfiles(simul,mode, \
                   width=None,XPedStart=None,XPedStop=None, \
                   TePed=None, ddx_TePed=None, TeLCFS=None, ddx_TeCore=None, ddx_TeSOL=None, \
                   TiPed=None, ddx_TiPed=None, TiLCFS=None, ddx_TiCore=None, ddx_TiSOL=None,
                   nPed=None, ddx_nPed=None, nLCFS=None, ddx_nCore=None, ddx_nSOL=None,
                   TeKind="bezier", TiKind="bezier", nKind="bezier",
                   TeTransWidth = None, TiTransWidth = None, nTransWidth = None,
                   **gridParams):
    # calculations needed for all profile generation scenarios
    XStart = simul.inputs.psiMin
    XStop = simul.inputs.psiMax
    (XPedStart,XPedStop,width) = setPedestalStartStop(XPedStart,XPedStop,width)

    updateDeltaOmega(simul)
    Delta = simul.Delta
    omega = simul.omega
    
    syncSpeciesParams(simul)
    
    if mode == "generateConstPhi":
        assert len(simul.species_list) == 2
        DeltaOver2omega = Delta/(2*omega)
        profileList = generateConstPhi(XStart,XStop,width,XPedStart,XPedStop, \
                                       TePed, ddx_TePed, TeLCFS, ddx_TeCore, ddx_TeSOL, \
                                       TiPed, ddx_TiPed, TiLCFS, ddx_TiCore, ddx_TiSOL, \
                                       nPed, ddx_nPed, nLCFS, ddx_nCore, ddx_nSOL, \
                                       TeKind, TiKind, nKind, \
                                       TeTransWidth, TiTransWidth, nTransWidth, \
                                       simul.species_list[0], \
                                       DeltaOver2omega)
    else:
        print "generateProfiles: Error: invalid mode"
        raise ValueError("Unrecognized mode option")
    
    # update nu_r with logLambda from the mid-pedestal T and n
    # using the electron profiles, but this is arbitrary and not fully satisfactory
    n, T = extractProfilesFromList(profileList,["n","T"],["e"])
    psiPedMid = (XPedStop -XPedStart)/2.0
    n = simul.nBar * n(psiPedMid)
    T = simul.TBar * T(psiPedMid)
    updateNuR(simul,n,T)

    # generate the grid
    if simul.inputs.psiGridType == 0:
        # we are uniform no matter what we say
        if gridParams["gridType"] != "uniform":
            print "!!!!!!!!!!!!!!!!"
            print "createProfiles WARNING: psiGridType is 0, but grid-type is not-uniform. Will override with uniform grid"
            print "!!!!!!!!!!!!!!!!"
            gridParams["gridType"] = "uniform"
    genGrid = ProfileGrid(simul.inputs.psi,simul.psiAHat,XPedStart,XPedStop,**gridParams)
    psi, psiAHatArray = genGrid.generate()
    profileList.append(psi)
    profileList.append(psiAHatArray)

    # write profiles to an hdf5 file
    writeProfileListToOutput(profileList,simul)
    
