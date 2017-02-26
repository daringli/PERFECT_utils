from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simulation
from createProfiles import createProfiles

import numpy
import shutil #copy to copy files

#simulation initialization
input_filename="input.namelist"
norm_filename="norms.namelist"
species_filename="species"
output_filename="perfectOutput.h5"

simul=normalized_perfect_simulation(input_filename,norm_filename,species_filename,None)
species=simul.species_list

width=0.05224
XStart = 0.9
XStop = 1.1
XPedStop = 1.0

# Te
TePed = 0.42
ddx_TeCore = -1.4
TeLCFS = 0.1
ddx_TeSOL = 0.0

#Ti
TiPed = 0.42
ddx_TiCore = -1.2
TiLCFS = 0.35
ddx_TiSOL = 0.0

#n
nPed=0.7 # 10^-20 m^{-3}
nLCFS = 0.06
ddx_nCore = -0.40616
ddx_nSOL = 0



grid={
    "nonuniform":True, "gridType":"IntArctan","transitionLength":0.025,"pedestalGridDensity":0.7,"leftShift":0.025*20,"rightShift":0.025*20,
                }


createProfiles(simul, mode="generateConstPhi", \
               width=width,XPedStart=None,XPedStop=XPedStop, \
               TePed=TePed, ddx_TePed=None, TeLCFS=TeLCFS, ddx_TeCore=ddx_TeCore, ddx_TeSOL=ddx_TeSOL, \
               TiPed=TiPed, ddx_TiPed=None, TiLCFS=TiLCFS, ddx_TiCore=ddx_TiCore, ddx_TiSOL=ddx_TiSOL, \
               nPed=nPed, ddx_nPed=None, nLCFS=nLCFS, ddx_nCore=ddx_nCore, ddx_nSOL=ddx_nSOL, \
               TeKind="bezier", TiKind="bezier", nKind="bezier", \
               **grid
           )

