from __future__ import division

from profile import Profile
from generator import Generator
from findInProfileList import extractProfilesFromList
import scipy.constants as C

class ProfileCollisionalityFromProfiles(Generator):
    """Profile generator that takes n profile of all but one species"""
    
    def __init__(self,Z=None,nu_r=None,species=""):
        self.Z=Z
        self.nu_r = nu_r
        self.profile = "collisionality"
        self.species = species
        
    def generate(self,profileList=[]):
        [n, ddx_n, T, ddx_T] = extractProfilesFromList(profileList,["n", "ddx_n", "T", "ddx_T"],self.species)
        p = Profile(self.nu_r * Z**4 * (n/T**2))
        ddx_p = Profile(self.nu_r * Z**4 * (ddx_n/T**2 - 2 * n*ddx_T/T**3))
        p.profile = self.profile
        p.species = self.species
        p.generator = 'collisionalityFromProfiles'
        p.generator_dict = vars(self)
        ddx_p.profile = "ddx_" + self.profile
        ddx_p.species = self.species
        ddx_p.generator = 'collisionalityFromProfiles'
        ddx_p.generator_dict = vars(self)
        return (p, ddx_p)
        
if __name__ == "__main__":
    import numpy
    import matplotlib.pyplot as plt
    from profileFromData import ProfileFromData

    
    x1 = numpy.linspace(0,1)
    x2 = numpy.linspace(0,1,57)
    y1 = x1
    y2 = x2**2

    RGen=ProfileFromData(x1,y1,kind=3,profile='R',species='none')
    omegatGen=ProfileFromData(x2,y2,kind=3,profile='omegat',species='Z')

    (R,ddx_R) = RGen.generate()
    (omegat,ddx_omegat) = omegatGen.generate()
    
    profiles = [R,ddx_R,omegat,ddx_omegat]
    
    VtGen = ProfileVtFromOmegat(species="Z")

    (Vt,ddx_Vt) = VtGen.generate(profiles)
    
    x = numpy.linspace(0,1,100)
    plt.plot(x,Vt(x))
    plt.plot(x,R(x)*omegat(x))
    plt.show()
    plt.plot(x,ddx_Vt(x))
    plt.plot(x,ddx_R(x)*omegat(x) + R(x)*ddx_omegat(x))
    plt.show()
