from profileFromData import ProfileFromData

import h5py


class ProfileFromSimulation(ProfileFromData):
    """Profile generator that takes the profile from a previous simulation. WARNING: will only work if the same grid is used in both simulations."""
    
    def __init__(self,directory,species_filename="species",profiles_filename="input.profiles.h5",psi_filename="psiAHat.h5",kind=3,profile='',species='', useSpecies=None,useProfile=None):
        if useSpecies is None:
            useSpecies = species
        self.useSpecies = useSpecies
            
        if useProfile is None:
            useProfile = profile
        self.useProfile = useProfile
            
        self.species_filename = directory + "/" + species_filename
        self.profiles_filename = directory + "/" + profiles_filename
        self.psi_filename = directory + "/" + psi_filename
        
        
        self.species_file = open(self.species_filename,'r')
        self.species = self.species_file.readline().splitlines()[0].split(',')
        
        self.species_dict={}
        for i,s in enumerate(self.species):
            self.species_dict[s] = i

        print self.species_dict
        
        self.h5file = h5py.File(self.profiles_filename,'r')            
        self.h5file2 = h5py.File(self.psi_filename,'r')            
        self.profile_dict = {"T": "THats",
                             "eta": "etaHats",
                             "n": "nHats",
                             "Phi": "PhiHat"
        }
        self.ddx_profile_dict = {"T": "dTHatdpsis",
                             "eta": "detaHatdpsis",
                             "n": "dnHatdpsis",
                             "Phi": "dPhiHatdpsi"
        }
        
        keys = self.h5file.keys()
        if len(keys) == 1:
            self.groupname = "/"  + keys[0] + "/" 
        else:
            raise ValueError("More than one group in hdf5 file! Cannot decide which one to use.")
        x = self.h5file2[self.groupname+"psiNArray"][()]

        if self.useProfile is not "Phi":
            y = self.h5file[self.groupname+self.profile_dict[self.useProfile]][()][:,self.species_dict[self.useSpecies]]
            ddx_y = self.h5file[self.groupname+self.ddx_profile_dict[self.useProfile]][()][:,self.species_dict[self.useSpecies]]
            
        else:
            y = self.h5file[self.groupname+self.profile_dict[self.useProfile]][()]
            ddx_y = self.h5file[self.groupname+self.ddx_profile_dict[self.useProfile]][()]
        
        ProfileFromData.__init__(self,x,y,kind,profile,species,ddx_y)

        
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    directory = '/home/bstefan/perfectSimuls/JET_82551_smallELMs/__NEWEST/N0'
    gen1 = ProfileFromSimulation(directory,profile='n',species='D')
    x = gen1.x
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    plt.plot(x,ddx_f1(x))

    directory = '/home/bstefan/perfectSimuls/JET_82551_smallELMs/__NEWEST/Nlow_master'
    gen2 = ProfileFromSimulation(directory,profile='n',species='D')
    x = gen2.x
    (f2, ddx_f2) = gen2.generate()
    plt.plot(x,f2(x))
    plt.plot(x,ddx_f2(x))
    
    
    plt.show()
    
