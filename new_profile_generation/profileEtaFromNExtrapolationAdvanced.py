from profile import Profile
from generator import Generator
from findInProfileList import extractProfilesFromList
from profileEtaFromPhiN import ProfileEtaFromPhiN
from mtanh import mtanh_transition
from mtanh import ddx_mtanh_transition
import numpy

class ProfileEtaFromNExtrapolationAdvanced(Generator):
    """Profile generator that takes a density profile, interval and extrapolation method."""
    
    def __init__(self, interval1,interval2, point, l=None, species = '', useSpecies=None, ratio=1.0, extrapolation=1):
        self.profile = "eta"
        self.species = species
        self.interval1 = interval1
        self.interval2 = interval2
        self.point = point
        if l is None:
            self.l = (self.interval2[0] - self.interval1[-1])/4.0
        else:
            self.l = l
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
            sample_points1 = numpy.linspace(self.interval1[0],self.interval1[1],self.order+1)
            sample_points2 = numpy.linspace(self.interval2[0],self.interval2[1],self.order+1)
            sample_n1 = [n(s_p) for s_p in sample_points1]
            poly1 = numpy.polyfit(sample_points1,sample_n1,self.order)
            sample_n2 = [n(s_p) for s_p in sample_points2]
            poly2 = numpy.polyfit(sample_points2,sample_n2,self.order)
            
            #add a term to the second polynomial to make them overlap at point
            diff = numpy.polyval(poly1,self.point) - numpy.polyval(poly2,self.point)

            eta1 = lambda x : numpy.polyval(poly1,x)
            ddx_eta1 = lambda x : numpy.polyval(numpy.polyder(poly1),x)
            eta2 = lambda x : diff + numpy.polyval(poly2,x)
            ddx_eta2 = lambda x : numpy.polyval(numpy.polyder(poly2),x)
            eta = mtanh_transition(eta2,eta1,self.l,self.point)
            ddx_eta = ddx_mtanh_transition(eta2,eta1,ddx_eta2,ddx_eta1,self.l,self.point)


        p = Profile(eta)
        ddx_p = Profile(ddx_eta)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'EtaFromNExtrapolationAdvanced'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'EtaFromNExtrapolationAdvanced'
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
    XStop =  1.5 #xPed + 2 * width
    
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
    genEta = ProfileEtaFromNExtrapolationAdvanced([0.94,0.95],[1.03,1.05],0.98,None,'')
    (Eta,ddx_Eta) = genEta.generate(profiles)
    profiles.append(Eta)
    profiles.append(ddx_Eta)

    plt.plot(x,Eta(x))
    plt.plot(x,n(x))
    plt.show()

    plt.plot(x,ddx_Eta(x))
    plt.plot(x,ddx_n(x))
    plt.show()
