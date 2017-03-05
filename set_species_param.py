import f90nml
import numpy
import os

def calc_species_param(species_list,species_filename,norm_filename):
    norm_file = f90nml.read(norm_filename)
    eBar = norm_file["normalizationParameters"]["eBar"]
    mBar = norm_file["normalizationParameters"]["mBar"]
    
    species_file = f90nml.read(species_filename)
    Z=numpy.array([species_file["speciesCharge"][x] for x in species_list])/eBar
    mHat=numpy.array([species_file["speciesMass"][x] for x in species_list])/mBar
    return [Z,mHat]

def list_to_str(l,delim=' '):
    #Not used, as this is now built into the PERFECT input class
    ret=''
    for entry in l:
        ret=ret+str(entry)+delim
    #remove last delimiter
    return ret.rsplit(delim,1)[0]

def set_species_param(species_list,species_filename,norm_filename,simulation):
    [Z,mHat]=calc_species_param(species_list,species_filename,norm_filename)
    simulation.inputs.charges=Z
    simulation.inputs.masses=mHat

def species_Z(species_list,species_filename="species_database.namelist"):
    species_filename= os.path.join(os.path.dirname(__file__), species_filename)
    species_file = f90nml.read(species_filename)
    Z=numpy.array([species_file["speciesZ"][x] for x in species_list])
    return Z


if __name__=="__main__":
    print list_to_str(['a','b',1])
    print list_to_str([2,'c','dddd'],'-')
