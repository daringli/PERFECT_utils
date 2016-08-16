This folder contains profiles generated with various settings with various versions of the profile generation code for the purpose of checking that the output is the same for the different (and any new) version.


These checks are done automatically by running "check.py". This script relies on hdf5 and h5py.

old profiles file:
input.profiles.h5_old
new profiles file:
input.profiles.h5

If old matches the new, the check passes.

For visual inspection, run:
    plot_inputs.py dir1 [dir2...]
and check generated .pdf (inputs.pdf, ddx_inputs.pdf, deltas.pdf)

TODO: eta_i for const Phi with m2tanh generated badly.

_______________________________

subdirectories:

_______________________________

1_const_ne : deuterium-helium with 1:1 core conc, from the constant n_e scan
-----
1_sameflux_const_ne : deuterium-helium with 1:1 core conc, from the constant n_e scan with sameflux option
-----
1_sameflux_const_Phi : deuterium-helium with 1:1 core conc, from the constant Phi scan with sameflux option
-----
1_oldsameflux_const_ne : deuterium-helium with 1:1 core conc, from the constant n_e scan with oldsameflux option, which gives the same ion temperature profile as the the const_Phi case.
-----
1_const_Phi : deuterium-helium with 1:1 core conc, from the constant Phi scan without sameflux option (constant T_i gradients)
-----
imp_baseline : deuterium-nitrogen, baseline from the PPCF 2016 paper (a constant Phi profile).
-----
imp_baseline-nonuniform : same but nonuniform grid
