import f90nml
import numpy

def calc_species_param(species_list,species_filename,norm_filename):
    norm_file = f90nml.read(norm_filename)
    eBar = norm_file["normalizationParameters"]["eBar"]
    mBar = norm_file["normalizationParameters"]["mBar"]
    
    species_file = f90nml.read(species_filename)
    Z=numpy.array([species_file["speciesCharge"][x] for x in species_list])/eBar
    mHat=numpy.array([species_file["speciesMass"][x] for x in species_list])/mBar
    return [Z,mHat]

def list_to_str(l,delim=' '):
    ret=''
    for entry in l:
        ret=ret+str(entry)+delim
    #remove last delimiter
    return ret.rsplit(delim,1)[0]

def set_species_param(species_list,species_filename,norm_filename,simulation):
    [Z,mHat]=calc_species_param(species_list,species_filename,norm_filename)
    simulation.inputs.changevar("speciesParameters","charges",list_to_str(Z))
    simulation.inputs.changevar("speciesParameters","masses",list_to_str(mHat))


if __name__=="__main__":
    print list_to_str(['a','b',1])
    print list_to_str([2,'c','dddd'],'-')
