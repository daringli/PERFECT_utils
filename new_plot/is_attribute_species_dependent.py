from array_rank import arraylist_rank
import numpy


def is_attribute_species_dependent(simulList,attrib):
    #quick-n-dirty check to see whether an attribute is species dependent
    #it will be assumed to be so if the lenght along the last axis of the attribute data
    #is equal to the number of species in that simulation for all simulations

    #warning: this will always be false for simulations with one species
    #TODO: species dependence of attributes should be encoded in simulation object
    return all([numpy.shape(getattr(simul,attrib))[-1]  == len(simul.species) for simul in simulList])



if __name__=="__main__":
    class placeholder:
        pass
    
    sim=placeholder()
    sim.data=numpy.array([[1,2],[3,4]])
    sim.species=[1,2]
    simlist=[sim]
    print is_attribute_species_dependent(simlist,"data")
    sim.data=numpy.array([[1,2],[3,4]])
    sim.species=[1]
    simlist=[sim]
    print is_attribute_species_dependent(simlist,"data")
    sim.data=numpy.array([[[1,1],[2,2]],[[3,3],[4,4]]])
    sim.species=[1,2]
    simlist=[sim]
    print is_attribute_species_dependent(simlist,"data")
