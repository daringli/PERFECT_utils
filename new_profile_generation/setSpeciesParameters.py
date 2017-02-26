import f90nml
import numpy
import os

def speciesParams(speciesList,speciesFilename,normFilename):
    normFile = f90nml.read(normFilename)
    eBar = normFile["normalizationParameters"]["eBar"]
    mBar = normFile["normalizationParameters"]["mBar"]    
    speciesFile = f90nml.read(speciesFilename)
    Z=numpy.array([speciesFile["speciesCharge"][x] for x in speciesList])/eBar
    mHat=numpy.array([speciesFile["speciesMass"][x] for x in speciesList])/mBar
    return [Z,mHat]

def speciesZ(speciesList,speciesFilename="species_database.namelist"):
    # same as the Z from speciesParams, but assumes normalization with
    # elementary charge, such that no norms.namelist file is needed.
    # this is useful for prototyping, since it is unecessary to define norms
    # for made-up profiles
    speciesFilename= os.path.join(os.path.dirname(__file__), speciesFilename)
    eBar = 1.602176565e-19
    speciesFile = f90nml.read(speciesFilename)
    Z=numpy.array([speciesFile["speciesCharge"][x] for x in speciesList])/eBar
    return Z
