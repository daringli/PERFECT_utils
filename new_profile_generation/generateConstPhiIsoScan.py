from profile3Linear import Profile3Linear
from profilePhiFromEtaN import ProfilePhiFromEtaN
from profileEtaFromPhiN import ProfileEtaFromPhiN
from profileNFromQuasiNeutrality import ProfileNFromQuasiNeutrality

from setPedestalParams import setPedestalParams,setPedestalStartStop
from setSpeciesParameters import speciesZ

import matplotlib.pyplot as plt
import numpy

debug=False

# contant Phi, two species
def generateConstPhi(XStart=None,XStop=None,width=None,XPedStart=None,XPedStop=None, \
                     TePed=None, ddx_TePed=None, TeLCFS=None, ddx_TeCore=None, ddx_TeSOL=None, \
                     TiPed=None, ddx_TiPed=None, TiLCFS=None, ddx_TiCore=None, ddx_TiSOL=None,
                     nPed=None, ddx_nPed=None, nLCFS=None, ddx_nCore=None, ddx_nSOL=None,
                     TeKind="bezier", TiKind="bezier", nKind="bezier",
                     TeTransWidth = None, TiTransWidth = None, nTransWidth = None, 
                     ion = "D",
                     DeltaOver2omega = 1.0):

    Z = speciesZ([ion])[0]
    Zs=[Z,-1] # add e


    # OVERRIDE TEMPERATURE GENERATION TO ALWAYS YIELD THE SAME TEMPERATURE PROFILE
    # NO TEMPERATURE PEDESTAL

    Twidth = 0.05224 # base width. Should no longer correspond to a pedestal width, since SOL and ped have same gradient.
    (TiPed,ddx_TiPed, width, TiLCFS) = setPedestalParams(TiPed,ddx_TiPed, width, TiLCFS)
    (TiPed,ddx_TiPed, Twidth, TiLCFS) = setPedestalParams(TiPed,ddx_TiSOL, Twidth, None)
    TiGen = Profile3Linear(XStart,XStop,ddx_TiCore,ddx_TiSOL,Twidth,XPedStart,XPedStart +Twidth,TiPed,ddx_TiPed,TiLCFS,profile="T",species=ion,kind=TiKind, transWidth = TiTransWidth)
    
    (TePed,ddx_TePed, width, TeLCFS) = setPedestalParams(TePed, ddx_TePed, width, TeLCFS)
    (TePed,ddx_TePed, Twidth, TeLCFS) = setPedestalParams(TePed, ddx_TeSOL, Twidth, None)
    (TePed,ddx_TePed, Twidth, TeLCFS) = (0.42, -1.4, 0.05224, 0.346864)
    TeGen = Profile3Linear(XStart,XStop,ddx_TeCore,ddx_TeSOL,Twidth,XPedStart,XPedStart +Twidth,TePed,ddx_TePed,TeLCFS,profile="T",species="e",kind=TeKind, transWidth = TeTransWidth)
    
    
    niGen = Profile3Linear(XStart,XStop,ddx_nCore,ddx_nSOL,width,XPedStart,XPedStop,nPed,ddx_nPed,nLCFS,profile="n",species=ion,kind=nKind, transWidth = nTransWidth)
    neGen = ProfileNFromQuasiNeutrality(Zs,species="e")
    
    (nPed,ddx_nPed, width, nLCFS) = setPedestalParams(nPed,ddx_nPed, width,nLCFS)
    etawidth = 0.05224 # base width. Should no longer correspond to a pedestal width, since SOL and ped have same gradient.
    (etaPed,ddx_etaPed, etawidth, etaLCFS) = setPedestalParams(nPed,ddx_nSOL, etawidth,None)
    #etaTransWidth = (nTransWidth/width)*etawidth
    etaTransWidth = nTransWidth
    etaiGen = Profile3Linear(XStart,XStop,ddx_nCore,ddx_nSOL,etawidth,XPedStart,XPedStart + etawidth,etaPed,ddx_etaPed,YLCFS=etaLCFS,profile="eta",species=ion,kind=nKind, transWidth = etaTransWidth)

    PhiGen = ProfilePhiFromEtaN(Z,useSpecies=ion,DeltaOver2omega=DeltaOver2omega)
    etaeGen = ProfileEtaFromPhiN(-1,species="e",DeltaOver2omega=DeltaOver2omega)

    # gridGen = GridFromPedestalBoundaries
    
    generatorList = [TiGen,TeGen,etaiGen,niGen,neGen,PhiGen,etaeGen]

    profileList=[]
    for gen in generatorList:
        p,ddx_p = gen.generate(profileList)
        profileList.append(p)
        profileList.append(ddx_p)
    #[ni,ddx_ni,ne,ddx_ne,etai,ddx_etai,Ti,ddx_Ti,Te,ddx_Te,Phi,ddx_Phi] = profileLis

    if debug:
        # visualize
        x = numpy.linspace(XStart,XStop,50)
        N = len(profileList[0::2])
        for i,p in enumerate(profileList[0::2]):
            print p.profile + "_" + p.species

            plt.subplot(N,1,i+1)
            plt.plot(x,p(x))
            plt.ylabel(p.profile + "_" + p.species)
        plt.show()
    return profileList

if __name__ == "__main__":
    # generic
    width=0.05224
    XStart = 0.9
    XStop = 1.1
    XPedStop = 1.0

    # Te
    TePed = 0.42
    ddx_TeCore = -1.4
    TeLCFS = 0.1
    ddx_TeSOL = 0.0

    #Ti
    TiPed = 0.42
    ddx_TiCore = -1.2
    TiLCFS = 0.35
    ddx_TiSOL = 0.0

    #n
    nPed=0.7 # 10^-20 m^{-3}
    nLCFS = 0.06
    ddx_nCore = -0.40616
    ddx_nSOL = 0

    
    generateConstPhi(XStart=XStart,XStop=XStop,width=width,XPedStart=None,XPedStop=XPedStop, \
                     TePed=TePed, ddx_TePed=None, TeLCFS=TeLCFS, ddx_TeCore=ddx_TeCore, ddx_TeSOL=ddx_TeSOL, \
                     TiPed=TiPed, ddx_TiPed=None, TiLCFS=TiLCFS, ddx_TiCore=ddx_TiCore, ddx_TiSOL=ddx_TiSOL,
                     nPed=nPed, ddx_nPed=None, nLCFS=nLCFS, ddx_nCore=ddx_nCore, ddx_nSOL=ddx_nSOL,
                     TeKind="bezier", TiKind="bezier", nKind="bezier",
                     ion = "D"
    )
