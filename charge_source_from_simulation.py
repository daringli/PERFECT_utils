#!/usr/bin/python
from perfect_simulations_from_dirs import perfect_simulations_from_dirs
import sys

dirlist=sys.argv[1:]

normlist_filename,species_filename,psiN_to_psi = ("norms.namelist","species","psiAHat")
simulList=perfect_simulations_from_dirs(dirlist,normlist_filename,species_filename,psiN_to_psi)
simul = simulList[0]

simul.export_attribute_to_h5("source_charge","chargeSource","chargeSource.h5")
