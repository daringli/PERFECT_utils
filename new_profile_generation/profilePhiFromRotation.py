from __future__ import division

from profile import Profile
from generator import Generator

class ProfilePhiFromRotation(Generator):
    """Profile generator that takes toroidal rotation and gives Phi."""
    def __init__(self,x,useSpecies,kind=3):
        self.x = x
        self.useSpecies = useSpecies
        self.profile = "Phi"
        self.species = "none"
        if type(kind) is int:
            self.type = "smoothingSpline"
            self.order = kind
            self.smoothing = 0
        
    def generate(self,profileList=[]):
        wantedProfiles = ["Vt","Vp"]
        
        # extract Vt and Vp profile
        [Vt,Vp] = extractProfilesFromList(profileList,wantedProfiles,self.useSpecies)

        if Vtt is not None:
            self.Vt = UnivariateSpline(self.x, Vt, k=self.order,s=self.smoothing)
            self.ddx_Vt = self.Vt.derivative()
        
        if Vp is not None:
            self.Vp = UnivariateSpline(self.x, Vp, k=self.order,s=self.smoothing)
            self.ddx_Vp = self.Vp.derivative()
        else:
            wantedProfiles = ["collisionality","ddx_T"]
            [nu,ddx_T] = extractProfilesFromList(profileList,wantedProfiles,self.useSpecies)
            wantedProfiles = ["BpOutboard","I","FSAB2"]
            [I,BpOutboard,FSAB2]= extractProfilesFromList(profileList,wantedProfiles,"none")
            
            
        
