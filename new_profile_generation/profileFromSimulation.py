from profileFromData import ProfileFromData

import h5py


class ProfileFromSimulation(ProfileFromData):
    """Profile generator that takes the profile from a previous simulation"""
    
    def __init__(self,directory,output_filename="perfectOutput.h5",species_filename="species",groupname="/run  1/",kind=3,profile='',species=''):
        self.output_filename = directory + "/" + output_filename
        self.species_filename = directory + "/" + species_filename
        self.groupname = groupname
        self.output_file = h5py.File(self.output_filename,'r')
        self.species_file = open(self.species_filename,'r')
        
        self.species = self.species_file.readline().splitlines()[0].split(',')
        
        self.species_dict={}
        for i,s in enumerate(self.species):
            self.species_dict[s] = i

        print self.species_dict
        
        self.profile_dict = {"T": "THat",
                       "eta": "etaHat",
                       "n": "nHat",
                       "Phi": "PhiHat"
        }

        if profile is not "Phi":
            y = self.output_file[self.groupname+self.profile_dict[profile]][()][:,self.species_dict[species]]
        else:
            y = self.output_file[self.groupname+self.profile_dict[profile]][()]
        x = self.output_file[self.groupname+"psi"][()]

        print y
        ProfileFromData.__init__(self,x,y,kind,profile,species)

        
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    directory = '/home/bstefan/perfectSimuls/JET1_N/allU00'
    gen1 = ProfileFromSimulation(directory,profile='T',species='D')
    x = gen1.x
    (f1, ddx_f1) = gen1.generate()
    plt.plot(x,f1(x))
    plt.plot(x,ddx_f1(x))
    plt.show()
    
