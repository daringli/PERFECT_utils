from p_1d_plot import perfect_1d_plot
from p_2d_plot import perfect_2d_plot

import sys
from scipy.constants import pi
import numpy

dirlist=sys.argv[1:]
print dirlist

dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
vlines=[0.94927395957025573, 0.97463697978512787]
vlines2=[94.927395957025573, 97.463697978512787]

#newpsiPedMax = 0.94927395957025573 + (0.97463697978512787-0.94927395957025573)/numpy.sqrt(2)
#vlines=[0.94927395957025573, 0.97463697978512787,newpsiPedMax]

xlims=[0.84,1.11]

attribs=["global_term_multiplier"]
ylabels=[r"$C(\psi_N)$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="globalTermMultiplier",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)


attribs=["etaHat","nHat","THat","PhiHat"]
ylabels=[r"$\eta$",r"$n$",r"$T$",r"$\Phi$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="inputs",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)


attribs=["detaHatdpsiN","dnHatdpsiN","dTHatdpsiN","dPhiHatdpsiN"]
ylabels=[r"$d\eta/d\psi_N$",r"$dn/d\psi_N$",r"$dT/d\psi_N$",r"$d\Phi/d\psi_N$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="ddx_inputs",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)

attribs=["psiAHatArray"]
ylabels=[r"$dh/ds$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="psiAHatArray",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)



attribs=["deltaEta","deltaN","deltaT","U"]
ylabels=[r"$\delta_\eta$",r"$\delta_n$",r"$\delta_T$",r"$U$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="deltas",xattr="actual_psiN",xlims=xlims,hlines=[0],vlines=vlines,generic_labels=False,lg=False)


