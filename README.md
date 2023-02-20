# ephys-jl

Julia versions of some of the ephys analysis routines. Faster than the python versions, but
stripped-down also. This is beta code and is not fully tested yet.

Problems appear with PyPlot/PyCall. 
You may need to install an arm-64 version of miniconda:
brew install miniconda
But this doesn't work (12/1/2021). 

Start Julia
 use the shell to go to the ephys-jl directory (cd ephys-jl)impo ] pkg
 activate .
 instantiate


just do add HDF5_jll (yes, that is right)
OLD:
Do the "brew install hdf5" separately, and maybe use "reinstall" to force update.
be sure to install from a native M1 terminal, not the rosetta terminal (2023)
Set Julia path in .zshrc:
export JULIA_HDF5_PATH=/opt/homebrew/Cellar/hdf5/1.12.2_2
restart julia to be sure it can be added in the package manager:
] pkg > add HDF5



Check that the correct version of Julia is run. May need to modify .zshrc (or .bash_profile or .bash_rc).
V1.7 is the first version that runs natively on M1 hardware.

using Pkg
Pkg.activate(".")
(Pkg.instantiate(), Pkg.upgrade_manifest())activate

or ] (pkg)
activate
instantiate

