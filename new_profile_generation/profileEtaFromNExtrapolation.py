from profile import Profile
from generator import Generator
from findInProfileList import extractProfilesFromList
from profileEtaFromPhiN import ProfileEtaFromPhiN
import numpy

class ProfileEtaFromNExtrapolation(Generator):
    """Profile generator that takes a density profile, interval and extrapolation method."""
    
    def __init__(self, interval, species = '', useSpecies=None, ratio=1.0, extrapolation=1):
        self.profile = "eta"
        self.species = species
        self.interval = interval
        if type(extrapolation) is int:
            self.extrapolation = "polyfit"
            self.order = extrapolation
            self.smoothing = 0
        if useSpecies is None:
            self.useSpecies = self.species
        else:
            self.useSpecies = useSpecies
        self.ratio = ratio
        
    def generate(self,profileList=[]):
        wantedProfiles = ["n","ddx_n"]

        [_n,_ddx_n] = extractProfilesFromList(profileList,wantedProfiles,self.useSpecies)
        n = self.ratio * _n
        ddx_n = self.ratio * _ddx_n


        if self.extrapolation == "polyfit":
            sample_points = numpy.linspace(self.interval[0],self.interval[1],self.order+1)
            sample_n = [n(s_p) for s_p in sample_points]
            poly = numpy.polyfit(sample_points,sample_n,self.order)
            eta = lambda x : numpy.polyval(poly,x)
            ddx_eta = lambda x : numpy.polyval(numpy.polyder(poly),x)
            

        p = Profile(eta)
        ddx_p = Profile(ddx_eta)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'EtaFromNExtrapolation'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'EtaFromNExtrapolation'
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
    #genPhi = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_PhiCore,ddx_YSOL=ddx_PhiSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=PhiPed,ddx_YPed=ddx_PhiPed,transWidthToPedWidth=0.2,profile="Phi")
    #genT = Profile3Linear(XStart=XStart,XStop=XStop,ddx_YCore=ddx_TCore,ddx_YSOL=ddx_TSOL,width=width,XPedStart=XPedStart,XPedStop=XPedStop,YPed=TPed,ddx_YPed=ddx_TPed,transWidthToPedWidth=0.2,profile="T")
   
    
    (n, ddx_n) = genn.generate()
    #(Phi, ddx_Phi) = genPhi.generate()
    #(T, ddx_T) = genT.generate()
    profiles = [n,ddx_n]
    
    # potential
    Z = 1
    Delta = 0.0006  #based on He papar
    omega = Delta/2
    DeltaOver2omega = Delta/(2*omega)
    genEta = ProfileEtaFromNExtrapolation([0.9,0.93],'')
    (Eta,ddx_Eta) = genEta.generate(profiles)
    profiles.append(Eta)
    profiles.append(ddx_Eta)

    plt.plot(x,Eta(x))
    plt.plot(x,n(x))
    plt.show()

    plt.plot(x,ddx_Eta(x))
    plt.plot(x,ddx_n(x))
    plt.show()
