from __future__ import division

from profile import Profile
from generator import Generator
from scipy.interpolate import UnivariateSpline
from findInProfileList import extractProfilesFromList
from numpy import log
from setSpeciesParameters import speciesZ
from findUnsharedElements import findUnsharedElements

class ProfileNFromQuasiNeutrality(Generator):
    """Profile generator that takes n profile of all but one species"""
    
    def __init__(self,Zs,species=""):
        self.Zs = Zs
        self.Nspecies = len(self.Zs)
        self.profile = "n"
        self.species = species
        
    def generate(self,profileList=[]):
        #extract the densities for all species but the missing
        ns = extractProfilesFromList(profileList,"n","any")
        ddx_ns = extractProfilesFromList(profileList,"ddx_n","any")
        if len(ns) != self.Nspecies - 1:
            print "ProfileNFromQuasiNeutrality : error: number of provided densities is not equal to Nspecies - 1"
            raise ValueError('Incorrect number of densities')
        if len(ns) != len(ddx_ns):
            print "ProfileNFromQuasiNeutrality : error: number of provided densities is not equal to number of provided ddx_n"
            raise ValueError('Incorrect number of densities or ddx_n')
            
        # get the Zs of the species for which we have densities
        species = [n.species for n in ns]
        species_ddx = [ddx_n.species for ddx_n in ddx_ns]
        if len(findUnsharedElements(species,species_ddx)) > 0:
            print "ProfileNFromQuasiNeutrality : error: n and ddx_n are for different species"
            raise ValueError('Unmatching species between n and ddx_n')
        Zs = speciesZ(species)
        Zs_ddx = speciesZ(species_ddx)
            
        # identify element in self.Zs not in Zs
        Z = findUnsharedElements(self.Zs,Zs)
        if len(Z) > 1:
            print "ProfileNFromQuasiNeutrality : error: was unable to deduce charge of remaining species."
            raise ValueError('Incorrect species charges provided')
        else:
            Z = Z[0]
        
        n = lambda x : sum([(Za*na/(-Z)) for Za,na in zip(Zs,ns)])(x)
        ddx_n = lambda x : sum([(Za*ddx_na/(-Z)) for Za,ddx_na in zip(Zs_ddx,ddx_ns)])(x)
        
    
        p = Profile(n)
        ddx_p = Profile(ddx_n)
        p.profile = self.profile
        p.species = self.species
        p.generator = 'NFromQuasiNeutrality'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'NFromQuasiNeutrality'
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

    n2Ped=0.4
    ddx_n2Core=-0.1*nPed/width
    ddx_n2Ped=-0.1*nPed/width
    ddx_n2SOL=-0.05*nPed/width

    # bezier curves
    genn = Profile3Linear(XStart,XStop,XPedStart,XPedStop,nPed,ddx_nCore,ddx_nPed,ddx_nSOL,width,transWidthToPedWidth=0.2,profile="n",species="H")
    genn2 = Profile3Linear(XStart,XStop,XPedStart,XPedStop,n2Ped,ddx_n2Core,ddx_n2Ped,ddx_n2SOL,width,transWidthToPedWidth=0.2,profile="n",species="He")
    
    (n, ddx_n) = genn.generate()
    (n2, ddx_n2) = genn2.generate()
    profiles = [n,ddx_n,n2,ddx_n2]
    
    # potential
    Zs = [1,2,-1]
    Delta = 0.0006  #based on He papar
    omega = Delta/2
    genn3 = ProfileNFromQuasiNeutrality(Zs,species='e')
    (n3,ddx_n3) = genn3.generate(profiles)
    profiles.append(n3)
    profiles.append(ddx_n3)

    plt.plot(x,n3(x))
    plt.plot(x,(Zs[0]*n(x)+Zs[1]*n2(x))/(-Zs[2]))
    plt.show()
