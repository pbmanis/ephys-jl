# ephys-jl

Julia versions of some of the ephys analysis routines. Faster than the python versions, but
stripped-down also

Problems appear with PyPlot/PyCall. 
You may need to install an arm-64 version of miniconda:
brew install miniconda
But this doesn't work (12/1/2021). 


HDF5 can be tricky to get working on a mac with an M1 processor. Best installed via
homebrew. 

Do the "brew install hdf5" separately, and maybe use "reinstall" to force update.

Check that the correct version of Julia is run. May need to modify .zshrc (or .bash_profile or .bash_rc).
V1.7 is the first version that runs natively on M1 hardware.


