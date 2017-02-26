from profile import Profile
from scipy.interpolate import UnivariateSpline
from findInProfileList import extractProfilesFromList
from numpy import exp

class ProfileEtaFromPhiN(object):
    """Profile generator that takes eta and n profile of some species and gives Phi."""
    
    def __init__(self,Z,species="",DeltaOver2omega=1.0):
        self.Z = Z
        self.DeltaOver2Omega = DeltaOver2omega
        self.profile = "eta"
        self.species = species
        
    def generate(self,profileList=[]):
        #extract needed profiles
        wantedProfiles = ["T","ddx_T","n","ddx_n"]
                
        [T,ddx_T,n,ddx_n] = extractProfilesFromList(profileList,wantedProfiles,self.species)
        [Phi,ddx_Phi] = extractProfilesFromList(profileList,["Phi","ddx_Phi"],"none")
        

        eta=lambda x: n(x)*exp((self.Z/self.DeltaOver2Omega)*Phi(x)/T(x))
        ddx_eta=lambda x: (ddx_n(x) + (self.Z/self.DeltaOver2Omega)*(n(x)/T(x))*(ddx_Phi(x)-Phi(x)*ddx_T(x)/T(x)))*exp((self.Z/self.DeltaOver2Omega)*Phi(x)/T(x))
    
        p = Profile(eta)
        ddx_p = Profile(ddx_eta)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'EtaFromPhiN'
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'EtaFromPhiN'
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

    PhiPed=0.0
    ddx_PhiCore=0
    ddx_PhiPed=0.1*1.0/width
    ddx_PhiSOL=0.0

    TPed=0.9
    ddx_TCore=-0.1*TPed/width
    ddx_TPed=-0.1*TPed/width
    ddx_TSOL=-0.05*TPed/width
    
    # bezier curves
    genn = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_nCore,ddx_YSOL=ddx_nSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=nPed,ddx_YPed=ddx_nPed,transWidthToPedWidth=0.2,profile="n")
    genPhi = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_PhiCore,ddx_YSOL=ddx_PhiSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=PhiPed,ddx_YPed=ddx_PhiPed,transWidthToPedWidth=0.2,profile="Phi")
    genT = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_TCore,ddx_YSOL=ddx_TSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=TPed,ddx_YPed=ddx_TPed,transWidthToPedWidth=0.2,profile="T")
   
    
    (n, ddx_n) = genn.generate()
    (Phi, ddx_Phi) = genPhi.generate()
    (T, ddx_T) = genT.generate()
    profiles = [n,ddx_n,Phi,ddx_Phi,T,ddx_T]
    
    # potential
    Z = 1
    Delta = 0.0006  #based on He papar
    omega = Delta/2
    DeltaOver2omega = Delta/(2*omega)
    genEta = ProfileEtaFromPhiN(Z,'',DeltaOver2omega)
    (Eta,ddx_Eta) = genEta.generate(profiles)
    profiles.append(Eta)
    profiles.append(ddx_Eta)

    plt.plot(x,Eta(x))
    plt.show()
