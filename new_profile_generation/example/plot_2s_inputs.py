from p_1d_plot import perfect_1d_plot
from p_2d_plot import perfect_2d_plot

import sys
from scipy.constants import pi


dirlist=["."]
dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
#vlines=[94.927395957025573, 97.463697978512787]
vlines=[0.94776,1]

xlims=[0.84,1.11]
attribs=["etaHat","nHat","THat","PhiHat"]
ylabels=attribs
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="inputs",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

attribs=["deltaEta","deltaN","deltaT","U"]
ylabels=attribs
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="deltas",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

attribs=["detaHatdpsiN","dnHatdpsiN","dTHatdpsiN","dPhiHatdpsiN"]
ylabels=attribs
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="ddx_inputs",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

attribs=[r"numerical_detaHatdpsiN"]
ylabels=attribs
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="numerical_ddx_inputs",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

attribs=[r"psiAHatArray_normed"]
ylabels=attribs
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="psiAHatArray",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

