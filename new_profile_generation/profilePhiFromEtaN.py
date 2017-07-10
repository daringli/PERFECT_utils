from profile import Profile
from generator import Generator
from scipy.interpolate import UnivariateSpline
from findInProfileList import extractProfilesFromList
from numpy import log

class ProfilePhiFromEtaN(Generator):
    """Profile generator that takes eta and n profile of some species and gives Phi."""
    
    def __init__(self,Z,useSpecies,DeltaOver2omega=1.0):
        self.Z = Z
        self.DeltaOver2omega = DeltaOver2omega
        self.useSpecies = useSpecies
        self.profile = "Phi"
        self.species = "none"
        
    def generate(self,profileList=[]):
        #extract needed profiles
        wantedProfiles = ["T","ddx_T","n","ddx_n","eta","ddx_eta"]
                
        [T,ddx_T,n,ddx_n,eta,ddx_eta] = extractProfilesFromList(profileList,wantedProfiles,self.useSpecies)
        Phi=lambda x : log(eta(x)/n(x))*T(x)*self.DeltaOver2omega/self.Z
        ddx_Phi=lambda x : ((ddx_eta(x)/eta(x) - ddx_n(x)/n(x))*T(x) + log(eta(x)/n(x))*ddx_T(x))*self.DeltaOver2omega/self.Z
    
        p = Profile(Phi)
        ddx_p = Profile(ddx_Phi)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'PhiFromEtaN'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'PhiFromEtaN'
        ddx_p.generator_dict = vars(self)
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    from profile3Linear import Profile3Linear

    width=0.03
    XPedStop = 1.0
    XPedStart = XPedStop - width
    xPed = 0.97
    XStart = xPed - width
    XStop =  XPedStop #xPed + 2 * width
    
    x=numpy.linspace(XStart,XStop,50)

    nPed=0.4
    ddx_nCore=-0.1*nPed/width
    ddx_nPed=-nPed/width
    ddx_nSOL=-0.05*nPed/width

    etaPed=0.4
    ddx_etaCore=-0.1*nPed/width
    ddx_etaPed=-0.1*nPed/width
    ddx_etaSOL=-0.05*nPed/width

    TPed=0.9
    ddx_TCore=-0.1*TPed/width
    ddx_TPed=-0.1*TPed/width
    ddx_TSOL=-0.05*TPed/width
    
    # bezier curves
    genn = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_nCore,ddx_YSOL=ddx_nSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=nPed,ddx_YPed=ddx_nPed,transWidthToPedWidth=0.2,profile="n")
    geneta = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_etaCore,ddx_YSOL=ddx_etaSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=etaPed,ddx_YPed=ddx_etaPed,transWidthToPedWidth=0.2,profile="eta")
    genT = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_TCore,ddx_YSOL=ddx_TSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=TPed,ddx_YPed=ddx_TPed,transWidthToPedWidth=0.2,profile="T")
    #genn = Profile3Linear(XStart,XStop,XPedStart,XPedStop,nPed,ddx_nCore,ddx_nPed,ddx_nSOL,width,transWidthToPedWidth=0.2,profile="n")
    #geneta = Profile3Linear(XStart,XStop,XPedStart,XPedStop,etaPed,ddx_etaCore,ddx_etaPed,ddx_etaSOL,width,transWidthToPedWidth=0.2,profile="eta")
    #genT = Profile3Linear(XStart,XStop,XPedStart,XPedStop,TPed,ddx_TCore,ddx_TPed,ddx_TSOL,width,transWidthToPedWidth=0.2,profile="T")
    
    (n, ddx_n) = genn.generate()
    (eta, ddx_eta) = geneta.generate()
    (T, ddx_T) = genT.generate()
    profiles = [n,ddx_n,eta,ddx_eta,T,ddx_T]
    
    # potential
    Z = 1
    Delta = 0.0006  #based on He papar
    omega = Delta/2
    DeltaOver2omega=Delta/(2*omega)
    genPhi = ProfilePhiFromEtaN(Z,'',DeltaOver2omega)
    (Phi,ddx_Phi) = genPhi.generate(profiles)
    profiles.append(Phi)
    profiles.append(ddx_Phi)

    plt.plot(x,n(x))
    plt.plot(x,eta(x))
    
    plt.show()



    plt.plot(x,Phi(x))
    plt.show()
